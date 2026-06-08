from __future__ import annotations
from pathlib import Path
import re, random, subprocess
import numpy as np
import pandas as pd
from .utils import run_cmd


def _expected_complete_pair_count(clean_fa):
    """
    Count duplex pairs that have both non-empty A and B arms in the clean FASTA.
    This is used to decide whether an existing output TSV is complete enough to skip.
    """
    pairs = _load_pairs_from_clean_fasta(clean_fa)
    return sum(
        1
        for ab in pairs.values()
        if "A" in ab and "B" in ab and bool(ab["A"]) and bool(ab["B"])
    )


def _valid_tsv(path, required_cols=None, min_rows=1):
    """
    Return True only if a TSV exists, is readable, has the expected columns,
    and has at least min_rows data rows.

    Important: this prevents a failed/partial/empty output from being treated
    as complete. The write functions below write to *.tmp first and then rename
    atomically, which further reduces the risk of partial files on reruns.
    """
    path = Path(path)
    if not path.exists() or path.stat().st_size == 0:
        return False

    try:
        # For large files, do not load the full table if we only need to verify
        # existence/columns/row-count threshold. nrows=0 is allowed for min_rows=0.
        if min_rows and int(min_rows) > 0:
            df = pd.read_csv(path, sep="\t", nrows=int(min_rows))
        else:
            df = pd.read_csv(path, sep="\t", nrows=0)
    except Exception:
        return False

    if required_cols:
        missing = [c for c in required_cols if c not in df.columns]
        if missing:
            return False

    if min_rows and int(min_rows) > 0:
        return len(df) >= int(min_rows)

    return True


def _write_tsv_atomic(df, out):
    """Write TSV through a temporary file, then atomically replace the target."""
    out = Path(out)
    tmp = out.with_suffix(out.suffix + ".tmp")
    df.to_csv(tmp, sep="\t", index=False)
    tmp.replace(out)


def prepare_duplex_inputs(bedtools_exe, fasta, pairs, kept_pair_ids, outdir, tag, condition):
    outdir = Path(outdir)
    pairs = pairs[pairs['pair_id'].isin(kept_pair_ids)]
    rows=[]
    for _, r in pairs.iterrows():
        rows.append([r['A_chrom'],r['A_start'],r['A_end'],f"{r['pair_id']}|A",0,r['A_strand']])
        rows.append([r['B_chrom'],r['B_start'],r['B_end'],f"{r['pair_id']}|B",0,r['B_strand']])
    bed = outdir/f'duplex_arms.{tag}.{condition}.bed'
    pd.DataFrame(rows, columns=['chrom','start','end','name','score','strand']).to_csv(bed, sep='\t', header=False, index=False)
    cp = run_cmd(
        [bedtools_exe, 'getfasta', '-fi', str(fasta), '-bed', str(bed), '-s', '-name'],
        capture=True,
        cwd=str(outdir),
    )
    fa = outdir/f'duplex_arms.{tag}.{condition}.fa'; fa.write_text(cp.stdout)
    clean = outdir/f'duplex_arms.{tag}.{condition}.clean.fa'
    with clean.open('w') as fo, fa.open() as fi:
        for line in fi:
            if line.startswith('>'):
                h=line[1:].strip().split()[0].split('::',1)[0]
                h=re.sub(r'\([+-]\)$','',h)
                if '|' not in h: continue
                pair, arm = h.split('|',1); arm=arm[:1]
                if arm not in ('A','B'): continue
                fo.write(f'>{pair}|{arm}\n')
            else:
                fo.write(line.upper().replace('T','U'))
    return clean


def _load_pairs_from_clean_fasta(clean_fa):
    fa = Path(clean_fa).read_text().splitlines()
    pairs={}; pid=arm=None
    for line in fa:
        if line.startswith('>'):
            pid, arm = line[1:].strip().split('|',1)
            pairs.setdefault(pid,{})[arm] = ''
        else:
            if pid and arm:
                pairs[pid][arm] += line.strip()
    return pairs


def run_rnacofold(args, clean_fa, outdir, tag, condition):
    outdir=Path(outdir)
    pairs=_load_pairs_from_clean_fasta(clean_fa)
    expected_n = _expected_complete_pair_count(clean_fa)

    out_tsv=outdir/f'duplex_pairs.{tag}.{condition}.cofold_mfe.tsv'
    required_cols = ['pair_id','RNAcofold_MFE_kcalmol','lenA','lenB','len_total','MFE_norm_kcalpermkb']
    if _valid_tsv(out_tsv, required_cols=required_cols, min_rows=expected_n):
        print(f"[SKIP] Existing valid RNAcofold MFE file found: {out_tsv}")
        return pairs, out_tsv

    cofold_in=outdir/f'duplex_pairs.{tag}.{condition}.cofold.in'
    lengths=outdir/f'duplex_pairs.{tag}.{condition}.lengths.tsv'
    with cofold_in.open('w') as fo, lengths.open('w') as fl:
        fl.write('pair_id\tlenA\tlenB\tlen_total\n')
        for pid,ab in pairs.items():
            if 'A' in ab and 'B' in ab and ab['A'] and ab['B']:
                fo.write(f'>{pid}\n{ab["A"]}&{ab["B"]}\n')
                fl.write(f'{pid}\t{len(ab["A"])}\t{len(ab["B"])}\t{len(ab["A"])+len(ab["B"])}\n')
    with cofold_in.open("r") as fin:
      proc = subprocess.run(
          [args.rnacofold_exe, "--noPS"],
          stdin=fin,
          text=True,
          capture_output=True,
          check=True,
          cwd=outdir,
      )
    (outdir/f'duplex_pairs.{tag}.{condition}.cofold').write_text(proc.stdout)
    lines=proc.stdout.splitlines(); out=[]; cur=None
    for i,line in enumerate(lines):
        if line.startswith('>'): cur=line[1:].strip()
        elif cur:
            m = re.search(r'\(\s*([-+]?\d+(?:\.\d+)?)\s*\)', line) or (re.search(r'\(\s*([-+]?\d+(?:\.\d+)?)\s*\)', lines[i+1]) if i+1 < len(lines) else None)
            if m:
                out.append((cur, float(m.group(1)))); cur=None
    mfe=pd.DataFrame(out, columns=['pair_id','RNAcofold_MFE_kcalmol'])
    lens=pd.read_csv(lengths, sep='\t')
    X=mfe.merge(lens,on='pair_id', how='left')
    X['MFE_norm_kcalpermkb']=X['RNAcofold_MFE_kcalmol']/(X['len_total']/1000.0)
    _write_tsv_atomic(X, out_tsv)
    return pairs, out_tsv


def run_ddg(args, clean_fa, cofold_tsv, outdir, tag, condition):
    outdir=Path(outdir)
    expected_n = _expected_complete_pair_count(clean_fa)
    out=outdir/f'duplex_pairs.{tag}.{condition}.ddg.tsv'
    required_cols = [
        'pair_id','RNAcofold_MFE_kcalmol','lenA','lenB','len_total',
        'MFE_norm_kcalpermkb','RNAfold_A_MFE_kcalmol','RNAfold_B_MFE_kcalmol',
        'ddG_interaction_kcalmol','ddG_norm_kcalpermkb'
    ]
    if _valid_tsv(out, required_cols=required_cols, min_rows=expected_n):
        print(f"[SKIP] Existing valid ddG file found: {out}")
        return out

    pairs=_load_pairs_from_clean_fasta(clean_fa)
    Ain=outdir/f'duplex_pairs.{tag}.{condition}.A.fold.in'; Bin=outdir/f'duplex_pairs.{tag}.{condition}.B.fold.in'
    with Ain.open('w') as fA, Bin.open('w') as fB:
        for pid,ab in pairs.items():
            if 'A' in ab and 'B' in ab and ab['A'] and ab['B']:
                fA.write(f'>{pid}\n{ab["A"]}\n'); fB.write(f'>{pid}\n{ab["B"]}\n')
    with Ain.open("r") as finA:
        outA = subprocess.run(
            [args.rnafold_exe, "--noPS"],
            stdin=finA,
            text=True,
            capture_output=True,
            check=True,
            cwd=outdir,
        ).stdout

    with Bin.open("r") as finB:
        outB = subprocess.run(
            [args.rnafold_exe, "--noPS"],
            stdin=finB,
            text=True,
            capture_output=True,
            check=True,
            cwd=outdir,
        ).stdout
    def parse_fold(text):
        lines=text.splitlines(); out=[]; cur=None
        for ln in lines:
            if ln.startswith('>'): cur=ln[1:].strip()
            else:
                m=re.search(r'\(\s*([-+]?\d+(?:\.\d+)?)\s*\)', ln)
                if m and cur:
                    out.append((cur, float(m.group(1)))); cur=None
        return pd.DataFrame(out, columns=['pair_id','G_kcalmol'])
    A=parse_fold(outA).rename(columns={'G_kcalmol':'RNAfold_A_MFE_kcalmol'})
    B=parse_fold(outB).rename(columns={'G_kcalmol':'RNAfold_B_MFE_kcalmol'})
    C=pd.read_csv(cofold_tsv, sep='\t')
    X=C.merge(A,on='pair_id',how='left').merge(B,on='pair_id',how='left')
    X['ddG_interaction_kcalmol']=X['RNAcofold_MFE_kcalmol']-(X['RNAfold_A_MFE_kcalmol']+X['RNAfold_B_MFE_kcalmol'])
    X['ddG_norm_kcalpermkb']=X['ddG_interaction_kcalmol']/(X['len_total']/1000.0)
    _write_tsv_atomic(X, out)
    return out


def run_interface_bpp(clean_fa, outdir, tag, condition):
    import RNA
    outdir=Path(outdir)
    expected_n = _expected_complete_pair_count(clean_fa)
    out=outdir/f'duplex_pairs.{tag}.{condition}.interface_bpp.tsv'
    required_cols = ['pair_id','interface_bpp_sum','interface_bpp_max','interface_bpp_n']
    if _valid_tsv(out, required_cols=required_cols, min_rows=expected_n):
        print(f"[SKIP] Existing valid ViennaRNA interface BPP file found: {out}")
        return out

    pairs=_load_pairs_from_clean_fasta(clean_fa); rows=[]
    for pid,ab in pairs.items():
        if 'A' not in ab or 'B' not in ab: continue
        A=ab['A']; B=ab['B'];
        if not A or not B: continue
        seq=A+'&'+B; lenA=len(A); n=lenA+len(B)
        fc=RNA.fold_compound(seq); fc.pf(); bppm=fc.bpp(); one_based=(len(bppm)==n+1)
        def get_p(i,j): return float(bppm[i][j]) if one_based else float(bppm[i-1][j-1])

        # Do not store all probabilities in a Python list. Accumulate directly.
        # This preserves the same output columns while reducing memory overhead.
        p_sum=0.0; max_p=0.0; p_n=0
        for i in range(1,lenA+1):
            for j in range(lenA+1,n+1):
                p=get_p(i,j)
                if p>0:
                    p_sum += p
                    p_n += 1
                    if p > max_p:
                        max_p = p
        rows.append({'pair_id':pid,'interface_bpp_sum':float(p_sum),'interface_bpp_max':float(max_p),'interface_bpp_n':int(p_n)})
    _write_tsv_atomic(pd.DataFrame(rows), out)
    return out


def run_null_z(args, clean_fa, ddg_tsv, outdir, tag, condition):
    outdir=Path(outdir)
    expected_n = _expected_complete_pair_count(clean_fa)
    out=outdir/f'duplex_pairs.{tag}.{condition}.nullZ.tsv'
    required_cols = ['pair_id','ddG_mu_null','ddG_sd_null','ddG_Z']
    if _valid_tsv(out, required_cols=required_cols, min_rows=expected_n):
        print(f"[SKIP] Existing valid null-Z file found: {out}")
        return out

    pairs=_load_pairs_from_clean_fasta(clean_fa); obs=pd.read_csv(ddg_tsv, sep='\t')
    random.seed(int(args.null_seed))
    def dinuc_shuffle(seq):
        from collections import defaultdict
        edges=defaultdict(list)
        for a,b in zip(seq[:-1],seq[1:]): edges[a].append(b)
        for k in edges: random.shuffle(edges[k])
        s=seq[0]; out=[s]; cur=s
        for _ in range(len(seq)-1):
            if not edges[cur]: cur=next((k for k,v in edges.items() if v), cur)
            nxt=edges[cur].pop(); out.append(nxt); cur=nxt
        return ''.join(out)
    rows=[]
    for _, r in obs.iterrows():
        pid=r['pair_id']; A=pairs.get(pid,{}).get('A',''); B=pairs.get(pid,{}).get('B','')
        if not A or not B: continue
        tmpin=outdir/f'.null_{pid}.in'
        with tmpin.open('w') as fo:
            for k in range(int(args.null_n)):
                fo.write(f'>{pid}__null{k}\n{dinuc_shuffle(A)}&{dinuc_shuffle(B)}\n')
        with tmpin.open("r") as fin:
          proc = subprocess.run(
              [args.rnacofold_exe, "--noPS"],
              stdin=fin,
              text=True,
              capture_output=True,
              check=True,
              cwd=outdir,
          )
        mfes=[]; cur=None; lines=proc.stdout.splitlines()
        for i,ln in enumerate(lines):
            if ln.startswith('>'): cur=ln[1:].strip()
            elif cur:
                m = re.search(r'\(\s*([-+]?\d+(?:\.\d+)?)\s*\)', ln) or (re.search(r'\(\s*([-+]?\d+(?:\.\d+)?)\s*\)', lines[i+1]) if i+1 < len(lines) else None)
                if m: mfes.append(float(m.group(1))); cur=None
        ddg_null=np.array(mfes)-(float(r['RNAfold_A_MFE_kcalmol'])+float(r['RNAfold_B_MFE_kcalmol']))
        mu=float(np.nanmean(ddg_null)); sd=float(np.nanstd(ddg_null, ddof=1)) if len(ddg_null)>1 else np.nan
        z=(float(r['ddG_interaction_kcalmol'])-mu)/sd if sd and sd>0 else np.nan
        rows.append({'pair_id':pid,'ddG_mu_null':mu,'ddG_sd_null':sd,'ddG_Z':z})
    _write_tsv_atomic(pd.DataFrame(rows), out)
    return out


def run_intarna(args, clean_fa, outdir, tag, condition):
    outdir=Path(outdir)
    expected_n = _expected_complete_pair_count(clean_fa)
    out=outdir/f'duplex_pairs.{tag}.{condition}.IntaRNA.tsv'
    required_cols = ['pair_id','E']
    if _valid_tsv(out, required_cols=required_cols, min_rows=expected_n):
        print(f"[SKIP] Existing valid IntaRNA file found: {out}")
        return out

    pairs=_load_pairs_from_clean_fasta(clean_fa); rows=[]
    for pid,ab in pairs.items():
        if 'A' not in ab or 'B' not in ab: continue
        try:
            proc=subprocess.run([args.intarna_exe, '--query', ab['A'], '--target', ab['B'], '--outMode', 'C', '--outCsvCols', 'E'], text=True, capture_output=True, check=True)
            val=proc.stdout.strip().splitlines()[-1].strip()
            energy=float(val) if val not in ('','E') else np.nan
        except Exception:
            energy=np.nan
        rows.append({'pair_id':pid,'E':energy})
    _write_tsv_atomic(pd.DataFrame(rows), out)
    return out
