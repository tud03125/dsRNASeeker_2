from __future__ import annotations
from pathlib import Path
import pandas as pd
import numpy as np
from .utils import run_cmd


def condition_samples(samplesheet_df, condition):
    df = samplesheet_df[samplesheet_df['condition'].astype(str) == str(condition)].copy()
    if df.empty:
        raise ValueError(f'No samples found for condition {condition}')
    return df


def _run_multi(exe, bedfile, bw_files, out_npz, out_tsv):
    if not bw_files:
        raise ValueError(f'No bigWig files to summarize for {out_tsv}')
    run_cmd([exe, 'BED-file', '--BED', str(bedfile), '-b', *bw_files, '-out', str(out_npz), '--outRawCounts', str(out_tsv), '--smartLabels'])


def run_coverage(args, samplesheet_df, outdir, condition, tag, bed_use):
    outdir = Path(outdir)
    cond_samples = condition_samples(samplesheet_df, condition)
    fwd_list, rev_list = [], []
    for _, row in cond_samples.iterrows():
        sid = str(row['sample_id'])
        bam = Path(str(row['bam_path']))
        if not bam.exists():
            raise FileNotFoundError(f'Missing BAM: {bam}')
        fwd_bw = outdir / f'{sid}.fwd.bw'
        rev_bw = outdir / f'{sid}.rev.bw'
        if not fwd_bw.exists():
            run_cmd([args.bamcoverage_exe, '-b', str(bam), '-o', str(fwd_bw), '--filterRNAstrand', 'forward', '--binSize', '25', '--numberOfProcessors', '8'], cwd=str(outdir))
        if not rev_bw.exists():
            run_cmd([args.bamcoverage_exe, '-b', str(bam), '-o', str(rev_bw), '--filterRNAstrand', 'reverse', '--binSize', '25', '--numberOfProcessors', '8'], cwd=str(outdir))
        fwd_list.append(str(fwd_bw))
        rev_list.append(str(rev_bw))
    _run_multi(args.multibigwigsummary_exe, bed_use, fwd_list, outdir/f'fwd.{tag}.{condition}.npz', outdir/f'fwd.{tag}.{condition}.tsv')
    _run_multi(args.multibigwigsummary_exe, bed_use, rev_list, outdir/f'rev.{tag}.{condition}.npz', outdir/f'rev.{tag}.{condition}.tsv')
    return cond_samples['sample_id'].astype(str).tolist()

def fuse_pair_level(outdir, tag, condition, bedfile):
    outdir = Path(outdir)

    fwd = pd.read_csv(outdir / f'fwd.{tag}.{condition}.tsv', sep='\t')
    rev = pd.read_csv(outdir / f'rev.{tag}.{condition}.tsv', sep='\t')

    # Clean weird deepTools column names like #'chr', 'start', 'end'
    def clean_cols(df):
        def clean(c):
            c = str(c)
            if c.startswith('#'):
                c = c[1:]
            c = c.strip()
            if len(c) >= 2 and c[0] == c[-1] and c[0] in ("'", '"'):
                c = c[1:-1]
            return c.strip().lower()
        df.columns = [clean(c) for c in df.columns]
        return df

    fwd = clean_cols(fwd)
    rev = clean_cols(rev)

    # Require coordinate columns explicitly
    for df_name, df in [('fwd', fwd), ('rev', rev)]:
        for col in ['chr', 'start', 'end']:
            if col not in df.columns:
                raise ValueError(f"{df_name} summary is missing required column '{col}'. Found: {list(df.columns)}")

    bed = pd.read_csv(
        bedfile,
        sep='\t',
        header=None,
        names=['chrom', 'start', 'end', 'pair_id'],
        dtype={'chrom': str, 'start': int, 'end': int, 'pair_id': str}
    )

    # Normalize chromosome names so '1' and 'chr1' match
    def norm_chr(x):
        s = str(x).strip()
        if s.startswith("chr"):
            return s
        return "chr" + s

    fwd['chr_norm'] = fwd['chr'].map(norm_chr)
    rev['chr_norm'] = rev['chr'].map(norm_chr)
    bed['chrom_norm'] = bed['chrom'].map(norm_chr)

    fwd['start'] = pd.to_numeric(fwd['start'], errors='coerce').astype('Int64')
    fwd['end']   = pd.to_numeric(fwd['end'], errors='coerce').astype('Int64')
    rev['start'] = pd.to_numeric(rev['start'], errors='coerce').astype('Int64')
    rev['end']   = pd.to_numeric(rev['end'], errors='coerce').astype('Int64')

    # Merge by normalized genomic interval, not raw pair_region string
    key_cols_fwd = ['chr_norm', 'start', 'end']
    key_cols_bed = ['chrom_norm', 'start', 'end']

    idx = fwd[key_cols_fwd].merge(
        bed[key_cols_bed + ['pair_id']],
        left_on=key_cols_fwd,
        right_on=key_cols_bed,
        how='left'
    )

    fwd = fwd.join(idx['pair_id']).dropna(subset=['pair_id']).copy()
    rev = rev.loc[fwd.index].copy()

    sample_cols_fwd = [c for c in fwd.columns if c not in ('chr', 'start', 'end', 'chr_norm', 'pair_id')]
    sample_cols_rev = [c for c in rev.columns if c not in ('chr', 'start', 'end', 'chr_norm')]

    for c in sample_cols_fwd:
        fwd[c] = pd.to_numeric(fwd[c], errors='coerce')
    for c in sample_cols_rev:
        rev[c] = pd.to_numeric(rev[c], errors='coerce')

    fwd_mean = fwd[sample_cols_fwd].mean(axis=1, skipna=True)
    rev_mean = rev[sample_cols_rev].mean(axis=1, skipna=True)
    total = fwd_mean + rev_mean

    out = pd.DataFrame({
        'pair_id': fwd['pair_id'].astype(str),
        'chrom': fwd['chr_norm'].astype(str),
        'start': fwd['start'].astype('Int64'),
        'end': fwd['end'].astype('Int64'),
        f'{condition}_fwd_mean': fwd_mean,
        f'{condition}_rev_mean': rev_mean,
        f'{condition}_total': total,
        f'{condition}_fwd_frac': (fwd_mean / total).where(total > 0, np.nan),
        f'{condition}_both_strands': (fwd_mean > 1e-6) & (rev_mean > 1e-6),
    }).drop_duplicates(subset=['pair_id'], keep='first')

    out_path = outdir / f'pair_windows_strand_signal.{tag}.{condition}.tsv'
    out.to_csv(out_path, sep='\t', index=False)
    return out_path

def arm_aware_summaries(args, pairs, kept_pair_ids, outdir, tag, condition):
    outdir = Path(outdir)
    PAD = int(args.arm_pad)
    rows=[]
    pairs = pairs[pairs['pair_id'].isin(kept_pair_ids)]
    for _, r in pairs.iterrows():
        rows.append([r['A_chrom'], max(0,int(r['A_start'])-PAD), int(r['A_end'])+PAD, f"{r['pair_id']}|A"])
        rows.append([r['B_chrom'], max(0,int(r['B_start'])-PAD), int(r['B_end'])+PAD, f"{r['pair_id']}|B"])
    arms_bed = outdir / f'pair_arms.{tag}.{condition}.bed'
    pd.DataFrame(rows, columns=['chrom','start','end','name']).to_csv(arms_bed, sep='\t', header=False, index=False)
    cond_df = pd.read_csv(args.samplesheet, sep='\t')
    cond_df = cond_df[cond_df['condition'].astype(str) == str(condition)]
    fwd_list = [str(outdir/f"{sid}.fwd.bw") for sid in cond_df['sample_id'].astype(str)]
    rev_list = [str(outdir/f"{sid}.rev.bw") for sid in cond_df['sample_id'].astype(str)]
    _run_multi(args.multibigwigsummary_exe, arms_bed, fwd_list, outdir/f'fwd.arms.{tag}.{condition}.npz', outdir/f'fwd.arms.{tag}.{condition}.tsv')
    _run_multi(args.multibigwigsummary_exe, arms_bed, rev_list, outdir/f'rev.arms.{tag}.{condition}.npz', outdir/f'rev.arms.{tag}.{condition}.tsv')
    fwd = pd.read_csv(outdir/f'fwd.arms.{tag}.{condition}.tsv', sep='\t', dtype=str)
    rev = pd.read_csv(outdir/f'rev.arms.{tag}.{condition}.tsv', sep='\t', dtype=str)
    def clean_cols(df):
      def clean(c):
          c = str(c)
          if c.startswith('#'):
              c = c[1:]
          c = c.strip()
          if len(c) >= 2 and c[0] == c[-1] and c[0] in ("'", '"'):
              c = c[1:-1]
          return c.strip().lower()
      df.columns = [clean(c) for c in df.columns]
      return df
    fwd = clean_cols(fwd); rev = clean_cols(rev)
    if 'name' not in fwd.columns:
      bed = pd.read_csv(
          arms_bed,
          sep='\t',
          header=None,
          names=['chr', 'start', 'end', 'name'],
          dtype={'chr': str, 'start': int, 'end': int, 'name': str}
      )
  
      for df in (fwd, rev):
          if 'start' not in df.columns or 'end' not in df.columns:
              raise ValueError(
                  f"Arm summary TSV is missing expected coordinate columns after cleaning. "
                  f"Found columns: {list(df.columns)}"
              )
          df['start'] = pd.to_numeric(df['start'], errors='coerce').astype('Int64')
          df['end'] = pd.to_numeric(df['end'], errors='coerce').astype('Int64')
  
      fwd = fwd.merge(bed, on=['chr', 'start', 'end'], how='left')
      rev = rev.merge(bed, on=['chr', 'start', 'end'], how='left')
    fwd = fwd.dropna(subset=['name']).copy(); rev = rev.loc[fwd.index].copy()
    for c in fwd.columns:
        if c not in ('chr','start','end','name'): fwd[c]=pd.to_numeric(fwd[c], errors='coerce')
    for c in rev.columns:
        if c not in ('chr','start','end','name'): rev[c]=pd.to_numeric(rev[c], errors='coerce')
    sf=[c for c in fwd.columns if c not in ('chr','start','end','name')]
    sr=[c for c in rev.columns if c not in ('chr','start','end','name')]
    fwd_m=fwd[sf].mean(axis=1, skipna=True); rev_m=rev[sr].mean(axis=1, skipna=True); tot=fwd_m+rev_m
    def splitname(s):
        s='' if pd.isna(s) else str(s)
        if '|' not in s: return (None,None)
        pid, arm = s.split('|',1); arm=(arm or '').strip()[:1]
        if arm not in ('A','B'): return (None,None)
        return (pid,arm)
    pid_arm=fwd['name'].map(splitname)
    ok=pid_arm.map(lambda x:isinstance(x,tuple) and x[0] is not None)
    pair_ids=[p for p,a in pid_arm.loc[ok]]
    arms=[a for p,a in pid_arm.loc[ok]]
    A=pd.DataFrame({'pair_id':pair_ids,'arm':arms,'f':fwd_m.loc[ok].values,'r':rev_m.loc[ok].values,'t':tot.loc[ok].values})
    A_A=A.query('arm=="A"').set_index('pair_id')
    A_B=A.query('arm=="B"').set_index('pair_id')
    J=A_A.join(A_B, lsuffix='_A', rsuffix='_B', how='inner')
    def sgn(x):
        if not np.isfinite(x) or abs(x)<1e-12: return 0.0
        return 1.0 if x>0 else -1.0
    J['A_bias']=J['f_A']-J['r_A']; J['B_bias']=J['f_B']-J['r_B']
    J['A_sgn']=J['A_bias'].map(sgn); J['B_sgn']=J['B_bias'].map(sgn)
    J['arms_both_cov']=(J['t_A']>float(args.arm_min_cov))&(J['t_B']>float(args.arm_min_cov))
    J['arm_opposite']=(J['A_sgn']*J['B_sgn']==-1)
    out = J[['arm_opposite','arms_both_cov']].rename(columns={'arm_opposite':f'{condition}_arm_opposite','arms_both_cov':f'{condition}_arms_both_cov'}).reset_index()
    out_path = outdir / f'pair_arms_signal.{tag}.{condition}.tsv'
    out.to_csv(out_path, sep='\t', index=False)
    return out_path
