import os
import pandas as pd

from pathlib import Path

def get_file_list(main_path, fc, out_dir):
    """flowcell result, sample alignment and sample result directories contain directories for each flowcell/sample
    appended by "v##". The latest output does not necessarily have the greatest ## value. This function puts together a
    list of the latest run directories based on creation time stamp.

    Args:
        main_path (Path): Path to directory to search
        fc (string): Flowcell ID

    Returns: List

    """
    files = pd.DataFrame()
    files['BFX_ID'] = [s for s in os.listdir(main_path) if s.split('_')[0] in fc]
    files['SAMPLE_ID'] = files['BFX_ID'].str[:-4] # remove version tag

    files['TIME'] = [os.stat(main_path / s).st_ctime for s in files['BFX_ID'].tolist()]

    files.sort_values(by='TIME', inplace=True)
    files.drop_duplicates(subset='SAMPLE_ID', keep='last', inplace=True)

    files.to_csv(out_dir / f'{fc}_dirs.txt', sep='\t', index=None)

    return files['BFX_ID'].tolist()


def aggregate_files(main_path, file_list, out_dir, fc):

    agg_cnv = []
    agg_fp = []
    agg_gaps = []
    agg_general = []
    agg_refmat = []
    agg_rt = []
    for s in file_list:
        agg_cnv.append(pd.read_csv(main_path / s / 'results' / 'cnv.tsv', sep='\t', header=0))
        agg_fp.append(pd.read_csv(main_path / s / 'results' / 'fingerprint.tsv', sep='\t', header=0))
        agg_gaps.append(pd.read_csv(main_path / s / 'results' / 'gaps.tsv', sep='\t', header=0))
        agg_general.append(pd.read_csv(main_path / s / 'results' / 'general.tsv', sep='\t', header=0))
        agg_refmat.append(pd.read_csv(main_path / s / 'results' / 'refmat.tsv', sep='\t', header=0))
        agg_rt.append(pd.read_csv(main_path / s / 'results' / 'rt.tsv', sep='\t', header=0))

    all_cnv = pd.concat(agg_cnv, axis=0,)
    all_fp = pd.concat(agg_fp, axis=0)
    all_gaps = pd.concat(agg_gaps, axis=0)
    all_general = pd.concat(agg_general, axis=0)
    all_refmat = pd.concat(agg_refmat, axis=0)
    all_rt = pd.concat(agg_rt, axis=0)

    all_cnv.to_csv(out_dir /f'{fc}_cnv.tsv', sep='\t', index=None)
    all_fp.to_csv(out_dir /f'{fc}_fp.tsv', sep='\t', index=None)
    all_gaps.to_csv(out_dir /f'{fc}_gaps.tsv', sep='\t', index=None)
    all_general.to_csv(out_dir /f'{fc}_general.tsv', sep='\t', index=None)
    all_refmat.to_csv(out_dir /f'{fc}_refmat.tsv', sep='\t', index=None)
    all_rt.to_csv(out_dir /f'{fc}_rt.tsv', sep='\t', index=None)

def compare_cnv(out_dir, fc):
    cnv_old = pd.read_csv(out_dir / f'{fc}_v1.0.6_cnv.tsv', sep='\t', header=0)
    cnv_new = pd.read_csv(out_dir / f'{fc}_v1.1.0_cnv.tsv', sep='\t', header=0)

    cnv_old = cnv_old.sort_values(by=['BFX_NEBULA_BFXID', 'BFX_NEBULA_ALLELEID']).reset_index(drop=True)
    cnv_new = cnv_new.sort_values(by=['BFX_NEBULA_BFXID', 'BFX_NEBULA_ALLELEID']).reset_index(drop=True)
    cnv_old = cnv_old.reindex(sorted(cnv_old.columns), axis=1)
    cnv_new = cnv_new.reindex(sorted(cnv_new.columns), axis=1)

    cnv_new.to_csv('/mnt/ruo_rw/rnd/staff/nilanthy.balendra/agg_test/comp.cnv_new.tsv', sep='\t', index=None)
    cnv_old.to_csv('/mnt/ruo_rw/rnd/staff/nilanthy.balendra/agg_test/comp.cnv_old.tsv', sep='\t', index=None)

    cols = ['BFX_CALL_OVERLAPGAP', 'BFX_CALL_REQUESTED', 'BFX_NEBULA_ALLELEID', 'BFX_NEBULA_BFXID', 'BFX_NEBULA_CALL', 
    'BFX_NEBULA_CALLPLOIDY', 'BFX_NEBULA_CALLSTATUS', 'BFX_NEBULA_CALLTYPE', 'BFX_NEBULA_DISORDER', 'BFX_NEBULA_GENENAME','BFX_NEBULA_HGVSC', 'BFX_NEBULA_REGION', 'BFX_NEBULA_ZYGOSITY', 
    'BFX_RESULT_DEFAULTCONFIRM', 'BFX_RESULT_HOLDFORREVIEW', 'BFX_RESULT_MULTIPLEALLELES', 'BFX_RESULT_REVIEWNEEDED']

    diff_col = [s for s in cnv_new.columns if s not in cols]
    print('cnv exclude')
    print(diff_col)
    pd.testing.assert_frame_equal(cnv_new[cols], cnv_old[cols])

def compare_rt(out_dir, fc):
    rt_old = pd.read_csv(out_dir / f'{fc}_v1.0.6_rt.tsv', sep='\t', header=0)
    rt_new = pd.read_csv(out_dir / f'{fc}_v1.1.0_rt.tsv', sep='\t', header=0)

    changed_vars = ['VSIUIJMY2S', 'VSI86YLIFK', 'VSCJ214FZP', 'VS7STJZ1F2', 'VSDJIGLZPQ',  'VSDAACPGML', 'VSSDKF4IAX', 'VS19NMNYM5', 'VSXTKJSYJK', 'VSEMRG78Y9', 'VSKMF3WCYB', 'VS22MP5ZGX']

    rt_old = rt_old.loc[~rt_old['BFX_RT_ALLELEID'].isin(changed_vars + ['CYP21A2_CYP21A1P_P453S', 'CYP21A2_CYP21A1P_V281L', 'CYP21A2_CYP21A1P_V282L'])]
    rt_new = rt_new.loc[~rt_new['BFX_RT_ALLELEID'].isin(changed_vars + ['CYP21A2_CYP21A1P_P453S', 'CYP21A2_CYP21A1P_V281L', 'CYP21A2_CYP21A1P_V282L'])]
    rt_old = rt_old.sort_values(by=['BFX_RT_BFXID', 'BFX_RT_ALLELEID']).reset_index(drop=True)
    rt_new = rt_new.sort_values(by=['BFX_RT_BFXID', 'BFX_RT_ALLELEID']).reset_index(drop=True)

    # drop_cols = ['BFX_RT_VARIANTQUAL', 'BFX_RT_CALLQUAL', 'BFX_RT_ALLELEREADS', 'BFX_RT_REFERENCEREADS', 'BFX_RT_OTHERREADS', 'BFX_COVATRON_REGION', 'BFX_COVATRON_REFALLELE','BFX_COVATRON_REFALLELECOUNT','BFX_COVATRON_ALTALLELE','BFX_COVATRON_ALTALLELECOUNT','BFX_COVATRON_OTHERALLELECOUNT' ]
    # pd.testing.assert_frame_equal(rt_new.drop(columns=drop_cols), rt_old.drop(columns=drop_cols))


    cols = ['BFX_RT_BFXID', 'BFX_RT_ALLELEID', 'BFX_RT_CALLSTATUS', 'BFX_RT_ALLELEPLOIDY', 'BFX_RT_REFERENCEPLOIDY', 'BFX_RT_OTHERPLOIDY',
    'BFX_RT_PLOIDY',
    'BFX_RT_GENENAME',
    'BFX_RT_DISORDER',
    'BFX_RT_ZYGOSITY',
    'BFX_RT_XVARCALL',
    'BFX_CALL_OVERLAPGAP',
    'BFX_RESULT_MULTIPLEALLELES',
    'BFX_RESULT_DEFAULTCONFIRM',
    'BFX_RESULT_REVIEWNEEDED',
    'BFX_RESULT_HOLDFORREVIEW', 'BFX_RT_ALLELEREADS', 'BFX_RT_REFERENCEREADS', 'BFX_RT_OTHERREADS', 'BFX_RT_COVERAGE', 'BFX_COVATRON_REGION', 'BFX_COVATRON_REFALLELE', 'BFX_COVATRON_REFALLELECOUNT', 'BFX_COVATRON_ALTALLELE', 'BFX_COVATRON_ALTALLELECOUNT', 'BFX_COVATRON_OTHERALLELECOUNT', 'BFX_COVATRON_UNINFORMEDCOUNT', 'BFX_COVATRON_DISCORDANTCOUNT', 'BFX_COVATRON_TOTALCOUNT', 'BFX_RT_VARIANTQUAL', 'BFX_RT_CALLQUAL', 'BFX_RT_CALLTYPE', 'BFX_RT_HGVSC', 'BFX_RT_ALLELERATIO']

    exclude = [c for c in rt_new.columns if c not in cols]
    print('these RT cols are excluded')
    print(exclude)

    pd.testing.assert_frame_equal(rt_new[cols], rt_old[cols])

def compare_general(out_dir, fc):
    gen_old = pd.read_csv(out_dir / f'{fc}_v1.0.6_general.tsv', sep='\t', header=0)
    gen_new = pd.read_csv(out_dir / f'{fc}_v1.1.0_general.tsv', sep='\t', header=0)

    gen_old = gen_old.sort_values(by=['LIS_REQUEST_BFXID']).reset_index(drop=True)
    gen_new = gen_new.sort_values(by=['LIS_REQUEST_BFXID']).reset_index(drop=True)

    drop_cols = ['BFX_ALIGNMENT_PAIRPFALIGNEDBASES', 'BFX_ALIGNMENT_PAIRPFINDELRATE', 'BFX_ALIGNMENT_PAIRFRACTIONCHIMERAS', 'BFX_ALIGNMENT_PAIRSTRANDBALANCE', 
    'BFX_ALIGNMENT_PAIRPFHQALIGNEDREADS', 'BFX_ALIGNMENT_PAIRPFALIGNEDREADS', 'BFX_ALIGNMENT_PAIRPFHQALIGNEDBASES', 'BFX_ALIGNMENT_PAIRALIGNEDIMPROPERPAIRS', 'BFX_ALIGNMENT_PAIRALIGNEDPROPERPAIRS', 'BFX_ALIGNMENT_PAIRPFHQALIGNEDQ20BASES',
    'BFX_GAPSRT_CALLABLE', 'BFX_GAPSRT_CALLABLECDS', 'BFX_RESULT_CREATEDDATETIME', 'LIS_RESULT_REQUESTGUID', 'LIS_REQUEST_CREATEDDATETIME', 'LIS_PIPELINE_VERSION', 'LIS_RESULT_RESULTGUID', 'BFX_GAPS_CALLABLECDS', 'BFX_GAPS_CALLABLE', 'BFX_RESULT_GAPSQCPASS',
    'BFX_RESULT_SAMPLEALIGNMENTRESULTURI', 'BFX_NEBULA_CHRXPLOIDY', 'BFX_NEBULA_CHRXPLOIDYTVALUE', 'BFX_NEBULA_CHRYPLOIDY', 'BFX_NEBULA_CHRYPLOIDYSE', 'BFX_NEBULA_CHRYPLOIDYTVALUE', 'BFX_RESULT_SAMPLEQCFAILUREREASON', 'BFX_NEBULA_OUTLIERSITECOUNT']
    pd.testing.assert_frame_equal(gen_new[gen_old.columns].drop(columns=drop_cols), gen_old.drop(columns=drop_cols))


def compare_fp(out_dir, fc):
    fp_old = pd.read_csv(out_dir / f'{fc}_v1.0.6_fp.tsv', sep='\t', header=0)
    fp_new = pd.read_csv(out_dir / f'{fc}_v1.1.0_fp.tsv', sep='\t', header=0)

    fp_old = fp_old.sort_values(by=['LIS_REQUEST_BFXID', 'BFX_IDENTITY_VARIANTID']).reset_index(drop=True)
    fp_new = fp_new.sort_values(by=['LIS_REQUEST_BFXID', 'BFX_IDENTITY_VARIANTID']).reset_index(drop=True)
    pd.testing.assert_frame_equal(fp_old, fp_new)

def compare_metrics(out_dir, fc):
    fc_metrics_old = pd.read_csv(out_dir / fc / 'v1.0.6' / 'flowcell_metrics.tsv', sep='\t', header=0)
    fc_metrics_new = pd.read_csv(out_dir / fc / 'v1.1.0' / 'flowcell_metrics.tsv', sep='\t', header=0)

    drop_cols = ['BFX_RESULT_CREATEDDATETIME', 'LIS_RESULT_REQUESTGUID', 'LIS_REQUEST_CREATEDDATETIME', 'LIS_PIPELINE_VERSION', 'LIS_RESULT_RESULTGUID', 'BFX_RESULT_FLOWCELLRESULTURI']
    pd.testing.assert_frame_equal(fc_metrics_old.drop(columns=drop_cols), fc_metrics_new[fc_metrics_old.columns].drop(columns=drop_cols))

    b_metrics_old = pd.read_csv(out_dir / fc / 'v1.0.6' / 'batch_metrics.tsv', sep='\t', header=0)
    b_metrics_new = pd.read_csv(out_dir / fc / 'v1.1.0' / 'batch_metrics.tsv', sep='\t', header=0)

    b_metrics_old = b_metrics_old.sort_values(by=['LIS_RESULT_NAME']).reset_index(drop=True)
    b_metrics_new = b_metrics_new.sort_values(by=['LIS_RESULT_NAME']).reset_index(drop=True)

    drop_cols = ['BFX_RESULT_CREATEDDATETIME', 'LIS_RESULT_REQUESTGUID', 'LIS_REQUEST_CREATEDDATETIME', 'LIS_PIPELINE_VERSION', 'LIS_RESULT_RESULTGUID', 'BFX_RESULT_SAMPLESFAILEDGAPSQCPASS', 'BFX_RESULT_FRACTIONSAMPLESFAILEDGAPSQCPASS']
    pd.testing.assert_frame_equal(b_metrics_old.drop(columns=drop_cols), b_metrics_new.drop(columns=drop_cols))

    s_metrics_old = pd.read_csv(out_dir / fc / 'v1.0.6' / 'sample_metrics.tsv', sep='\t', header=0)
    s_metrics_new = pd.read_csv(out_dir / fc / 'v1.1.0' / 'sample_metrics.tsv', sep='\t', header=0)

    s_metrics_old = s_metrics_old.sort_values(by=['LIS_RESULT_NAME']).reset_index(drop=True)
    s_metrics_new = s_metrics_new.sort_values(by=['LIS_RESULT_NAME']).reset_index(drop=True)

    drop_cols = ['BFX_RESULT_CREATEDDATETIME', 'LIS_RESULT_REQUESTGUID', 'LIS_REQUEST_CREATEDDATETIME', 'LIS_PIPELINE_VERSION', 'LIS_RESULT_RESULTGUID', 'BFX_ALIGNMENT_PAIRPFALIGNEDBASES', 'BFX_ALIGNMENT_PAIRPFINDELRATE', 'BFX_ALIGNMENT_PAIRFRACTIONCHIMERAS', 'BFX_ALIGNMENT_PAIRSTRANDBALANCE', 
    'BFX_ALIGNMENT_PAIRPFHQALIGNEDREADS', 'BFX_ALIGNMENT_PAIRPFALIGNEDREADS', 'BFX_ALIGNMENT_PAIRPFHQALIGNEDBASES', 'BFX_ALIGNMENT_PAIRALIGNEDIMPROPERPAIRS', 'BFX_ALIGNMENT_PAIRALIGNEDPROPERPAIRS', 'BFX_ALIGNMENT_PAIRPFHQALIGNEDQ20BASES', 'BFX_GAPSRT_CALLABLE', 'BFX_GAPSRT_CALLABLECDS', 'BFX_RESULT_CREATEDDATETIME', 'LIS_RESULT_REQUESTGUID', 'LIS_REQUEST_CREATEDDATETIME', 'LIS_PIPELINE_VERSION', 'LIS_RESULT_RESULTGUID', 'BFX_GAPS_CALLABLECDS', 'BFX_GAPS_CALLABLE', 'BFX_RESULT_GAPSQCPASS',
    'BFX_RESULT_SAMPLEALIGNMENTRESULTURI']
    pd.testing.assert_frame_equal(s_metrics_old.drop(columns=drop_cols), s_metrics_new[s_metrics_old.columns].drop(columns=drop_cols))
    print('wooooo')



def main():
    #fc = 'HCNJLDRXX'
    #fc = 'HCWM2DRXX'
    fc = 'HCNLGDRXX'
    main_dir = Path('/mnt/bfx_analysis_ro/analysis_data/progenity_workflow_levitate/sample_result')
    out_dir = Path('/mnt/ruo_rw/rnd/staff/nilanthy.balendra/agg_test')

    # # aggregate the old samples using the old file liist
    # with open(f'/mnt/ruo_rw/rnd/SQA_Data/output/progenity_workflow_levitate/vv_report_v1.0.6/objective_4/dirs/{fc}_dirs.txt') as f:
    #     content = f.readlines()
    
    # samples_old = [x.strip() for x in content] 
    # aggregate_files(main_dir, samples_old, out_dir, f'{fc}_v1.0.6')

    # #aggregate the latest stuff (candidate release)
    # samples_new = get_file_list(main_dir, fc, out_dir / 'dirs')
    # aggregate_files(main_dir, samples_new, out_dir, f'{fc}_v1.1.0')

    compare_cnv(out_dir, fc)
    compare_rt(out_dir, fc)
    compare_fp(out_dir, fc)
    compare_general(out_dir, fc)

    # compare_metrics(Path('/mnt/ruo_rw/rnd/SQA_Data/output/progenity_workflow_levitate/vv_report_v1.1.0/objective_8'), fc)




if __name__ == '__main__':
    main()