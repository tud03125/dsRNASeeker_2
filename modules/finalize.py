from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd


def finalize_condition(outdir, pairs, txmap_path, tag, condition, args):
    outdir=Path(outdir)
    COFOLD_STRONG=float(args.cofold_strong); COFOLD_MODERATE=float(args.cofold_moderate)
    strand=pd.read_csv(outdir/f'pair_windows_strand_signal.{tag}.{condition}.tsv', sep='\t')
    cofold_path=outdir/f'duplex_pairs.{tag}.{condition}.cofold_mfe.tsv'
    ddg_path=outdir/f'duplex_pairs.{tag}.{condition}.ddg.tsv'
    mfe=pd.read_csv(ddg_path if ddg_path.exists() else cofold_path, sep='\t')
    AtoI_fn=outdir/f'AtoI_counts_window.{tag}.{condition}.tsv'; REDI_fn=outdir/f'REDI_counts_window.{tag}.{condition}.tsv'
    AtoI=pd.read_csv(AtoI_fn, sep='\t') if AtoI_fn.exists() else pd.DataFrame(columns=['pair_id','AtoI_hits_window'])
    REDI=pd.read_csv(REDI_fn, sep='\t') if REDI_fn.exists() else pd.DataFrame(columns=['pair_id','REDI_hits_window'])

    # Orientation is already computed in modules.pairs.classify_orientations().
    # Do not recompute through transcript-map merges here; that can expand rows
    # if TE_name appears more than once in txmap.
    meta = pairs.copy().reset_index(drop=True)

    if 'genomic_orientation' not in meta.columns:
        meta['genomic_orientation'] = (meta['A_strand'].astype(str) != meta['B_strand'].astype(str)).map({True:'inverted',False:'direct'})

    if 'transcript_orientation' not in meta.columns:
        meta['transcript_orientation'] = meta['genomic_orientation']
    X=strand.merge(mfe, on='pair_id', how='left').merge(meta[['pair_id','A_repFamily','A_repName','B_repFamily','B_repName','A_SYMBOL','B_SYMBOL','A_annotation','B_annotation','genomic_orientation','transcript_orientation']], on='pair_id', how='left').merge(AtoI, on='pair_id', how='left').merge(REDI, on='pair_id', how='left')
    for opt in ['nullZ','interface_bpp','IntaRNA']:
        p=outdir/f'duplex_pairs.{tag}.{condition}.{opt}.tsv'
        if p.exists(): X=X.merge(pd.read_csv(p, sep='\t'), on='pair_id', how='left')
    arms_path=outdir/f'pair_arms_signal.{tag}.{condition}.tsv'
    if arms_path.exists():
        X=X.merge(pd.read_csv(arms_path, sep='\t'), on='pair_id', how='left')
    else:
        X[f'{condition}_arm_opposite']=np.nan; X[f'{condition}_arms_both_cov']=np.nan
    def energy_bin(e):
        try: e=float(e)
        except Exception: return 'NA'
        if e<=-60: return 'very_strong'
        if e<=COFOLD_STRONG: return 'strong'
        if e<=COFOLD_MODERATE: return 'moderate'
        return 'weak'
    X['cofold_energy_bin']=X['RNAcofold_MFE_kcalmol'].apply(energy_bin)
    def conf_row(r):
        orient_ok=(r.get('genomic_orientation')=='inverted') or (r.get('transcript_orientation')=='inverted')
        both_ok=bool(r.get(f'{condition}_both_strands', False))
        arm_ok=bool(r.get(f'{condition}_arm_opposite', False)) and bool(r.get(f'{condition}_arms_both_cov', False))
        mfe=r.get('RNAcofold_MFE_kcalmol', np.nan)
        if orient_ok and (both_ok or arm_ok) and pd.notna(mfe):
            if mfe<=COFOLD_STRONG: return 'high'
            elif mfe<=COFOLD_MODERATE: return 'probable'
        if (both_ok or arm_ok) and pd.notna(mfe) and mfe<-5: return 'possible'
        return 'uncertain'
    X['dsRNA_confidence']=X.apply(conf_row, axis=1)
    energy_points={'very_strong':3,'strong':2,'moderate':1,'weak':0,'NA':0}
    X['expr_points']=X[f'{condition}_both_strands'].fillna(False).astype(int)*2 + X.get(f'{condition}_arm_opposite', pd.Series(False,index=X.index)).fillna(False).astype(int)*2 + X.get(f'{condition}_arms_both_cov', pd.Series(False,index=X.index)).fillna(False).astype(int)
    X['energy_points']=X['cofold_energy_bin'].map(energy_points).fillna(0).astype(int)
    def editing_points(row):
        total=0.0
        for col in ['AtoI_hits_window','REDI_hits_window']:
            try:
                v=float(row.get(col, np.nan))
                if not np.isnan(v): total+=v
            except Exception: pass
        if total>=20: return 3
        elif total>=5: return 2
        elif total>=1: return 1
        return 0
    X['editing_points']=X.apply(editing_points, axis=1)
    def bias(frac):
        try: x=float(frac); return abs(x-0.5)
        except Exception: return None
    bias_avg=X.get(f'{condition}_fwd_frac', pd.Series(np.nan,index=X.index)).apply(bias)
    X['bias_penalty']=pd.Series(bias_avg).clip(0,0.5).fillna(0)*2
    X['rank_score']=X['expr_points']+X['energy_points']+X['editing_points']-X['bias_penalty']
    cols_order=['pair_id','A_SYMBOL','B_SYMBOL','A_annotation','B_annotation','A_repFamily','A_repName','B_repFamily','B_repName','genomic_orientation','transcript_orientation','RNAcofold_MFE_kcalmol','MFE_norm_kcalpermkb','RNAfold_A_MFE_kcalmol','RNAfold_B_MFE_kcalmol','ddG_interaction_kcalmol','ddG_norm_kcalpermkb','ddG_Z','interface_bpp_sum','interface_bpp_max','interface_bpp_n','cofold_energy_bin',f'{condition}_total',f'{condition}_fwd_frac',f'{condition}_both_strands',f'{condition}_arm_opposite',f'{condition}_arms_both_cov','AtoI_hits_window','REDI_hits_window','bias_penalty','expr_points','energy_points','editing_points','rank_score','dsRNA_confidence']
    for cand in ['E','E_total','Eall','E_hybrid','E_init','E_open','seedStart1','seedEnd1','seedStart2','seedEnd2']:
        if cand in X.columns and cand not in cols_order: cols_order.append(cand)
    for c in cols_order:
        if c not in X.columns: X[c]=np.nan
    out_sub=outdir/tag/condition; out_sub.mkdir(parents=True, exist_ok=True)
    X[cols_order].to_csv(out_sub/f'TEpair_dsRNA_master.{condition}.tsv', sep='\t', index=False)
    X_dedup=X[cols_order].drop_duplicates(subset='pair_id', keep='first').sort_values('rank_score', ascending=False)
    X_dedup.to_csv(out_sub/f'TEpair_dsRNA_master.{condition}.dedup.tsv', sep='\t', index=False)
    shortlist_cols=[c for c in ['pair_id','RNAcofold_MFE_kcalmol','MFE_norm_kcalpermkb','ddG_interaction_kcalmol','rank_score'] if c in X.columns]
    X.query("dsRNA_confidence == 'high'").sort_values('rank_score', ascending=False).head(50)[shortlist_cols].to_csv(out_sub/f'shortlist_high.{condition}.tsv', sep='\t', index=False)
    X.query("dsRNA_confidence == 'probable'").sort_values('rank_score', ascending=False).head(50)[shortlist_cols].to_csv(out_sub/f'shortlist_probable.{condition}.tsv', sep='\t', index=False)
    X.sort_values('rank_score', ascending=False).head(20)[[c for c in ['pair_id','RNAcofold_MFE_kcalmol','MFE_norm_kcalpermkb','cofold_energy_bin','dsRNA_confidence','rank_score'] if c in X.columns]].to_csv(out_sub/f'shortlist_top20_overall.{condition}.tsv', sep='\t', index=False)
