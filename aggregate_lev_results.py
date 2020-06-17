"""This script pulls together data from flowcells of interest for review

Arguments:
    - string: flowcells to pull data from, one value or comma separated list
    - path: output directory
    - string: output file prefix

Optional:
    --data, path to workflow results directory
    --f_result, path to flowcell_result directory
    --s_aligns, path to sample_alignment directory
    --s_result, path to sample_result directory

    Either provide --data or --f_result,--s_aligns, and --s_result. They are mutually exclusive options.

    -tsv, default value is False, if true, the following individual tsv files are outputted
        - {prefix}_sample_metrics.tsv
        - {prefix}_hs.tsv
        - {prefix}_other_cnv_pp2.tsv
        - {prefix}_smn_alpha.tsv
        - {prefix}_rt_pp1.tsv

The following file is always created and contains all contents in the above tsv files on separate sheets:
    - {prefix}_output.xlsx

Example commands:
python aggregate_lev_results.py H737HDRXX,HCYTJDRXX  /mnt/ruo_rw/rnd/staff/nilanthy.balendra/ my_prefix
--data /mnt/bfx_analysis_ro/analysis_data/progenity_workflow_levitate/ -tsv

python aggregate_results.py H737HDRXX,HCYTJDRXX /mnt/ruo_rw/rnd/staff/nilanthy.balendra/ my_prefix
--f_result /mnt/bfx_analysis_ro/analysis_data/progenity_workflow_levitate/flowcell_result/
--s_align /mnt/bfx_analysis_ro/analysis_data/progenity_workflow_levitate/sample_alignment/
--s_result /mnt/bfx_analysis_ro/analysis_data/progenity_workflow_levitate/sample_result/ -tsv

"""

import argparse
import os
import pandas as pd

from pathlib import Path

fp = ['SAMPLE_ID','VSE2LNKMRW', 'VSHKUT54EL', 'VS2CJ12KCY', 'VS5FN9VIQG', 'VSLJYQ6EVS', 'VSFMENYNMW', 'VSLV8ZHH33',
      'VSZ46X4R7A', 'VSLU4IJRML', 'VS93LGTXXY', 'VSS8GIC16K', 'VSU99R3RB5', 'VSX1RSIIP6']


def get_file_list(main_path, fc):
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

    return files['BFX_ID'].tolist()


def reformat_calls(calls, header):
    """Reformat the CNV output files according to variant region (SMN or ALPHA).

    Args:
        calls (dataframe): CNV calls aggregated together
        header: header to be added to reformated dataframe

    Returns: dataframe

    """
    new_format = pd.DataFrame()
    new_format['SAMPLE_ID'] = calls['SAMPLE_ID']
    new_format[header[0]] = calls['CALL']

    if 'SMN' in header[0]:
        new_format[header[1]] = calls['CALL_PLOIDY']
    else:
        new_format[header[1]] = calls['VARIANT_ID']

    new_format[header[2]] = calls['CALL_QUAL']
    new_format.set_index('SAMPLE_ID', inplace=True)

    return new_format


def aggregate_qc(flowcell_result_path, aligns_path, fc_id):
    """Given filepaths and a list of flowcells, sample_metrics.tsv files for each flowcell and hs_metrics.tsv files for
    each sample are aggregated and returned as two dataframes.

    Args:
        flowcell_result_path (Path): path to flowcell directory
        aligns_path (Path): path to sample_alignment directory
        fc_id (list): list of flowcell ids

    Returns: dataframe, dataframe

    """
    sample_dirs = get_file_list(aligns_path, fc_id)
    f_dirs = get_file_list(flowcell_result_path, fc_id)
    hs = []
    sample_metrics = pd.DataFrame()

    for d in sample_dirs:
        hs_metrics_path = aligns_path / d / 'metrics' / 'hs_metrics.tsv'

        if os.path.isfile(hs_metrics_path):
            hs_lines = open(hs_metrics_path).readlines()
            rest = hs_lines[7].split('\t')
            headings = hs_lines[6].split('\t')
            hs.append([d[:-4]] + rest)

    header = ['SAMPLE_ID'] + headings
    hs_metrics = pd.DataFrame(data=hs, columns=header)

    for d in f_dirs:
        sample_metrics_path = flowcell_result_path / d / 'sample_metrics.tsv'

        if os.path.isfile(sample_metrics_path):
            temp= pd.read_csv(sample_metrics_path, sep='\t', header=0)
            cols = list(temp)
            cols = sorted(cols)
            cols.insert(0, cols.pop(cols.index('LIS_REQUEST_BFXID'))) #move LIS_REQUEST_BFXID to first column position
            sample_metrics = pd.concat([sample_metrics, temp.loc[:, cols]], axis=0, sort=False)

    return sample_metrics, hs_metrics


def aggregate_calls(main_path, fc_id):
    """Given path to rt and cnv results directories, flowcell ids and version, aggregates calls and returns two
    dataframes.

    Args:
        main_path: Path to sample_result directory
        fc_id: flowcell ids of interest

    Returns: dataframe, dataframe

    """
    meta_data = pd.read_csv('Classic_Variant_Metadata_v6.txt', sep='\t', header=0)
    meta_data.set_index('VSID', inplace=True)

    dirs = get_file_list(main_path, fc_id)

    rt_dfs = []
    rt_dfs_nf = []
    agg_cnv_dfs = []

    for d in dirs:
        rt = main_path / d / 'results' / 'rt.tsv'
        cnv = main_path / d / 'cnv' / 'calls.tsv'

        if os.path.isfile(rt) and os.path.isfile(cnv):

            rt_output = pd.read_csv(rt, sep='\t', header=0)
            no_fp = rt_output.loc[~rt_output['BFX_RT_ALLELEID'].isin(fp)] #drop fingerprint variants
            no_fp = no_fp.loc[~no_fp['BFX_RT_ALLELEID'].str.contains('CYP21A2')]

            no_fp = no_fp[['BFX_RT_BFXID', 'BFX_RT_ALLELEID', 'BFX_RT_CALLSTATUS', 'BFX_RT_CALLTYPE', 'BFX_RT_VARIANTQUAL',
                             'BFX_RT_CALLQUAL', 'BFX_RT_ALLELEPLOIDY', 'BFX_RT_REFERENCEPLOIDY', 'BFX_RT_OTHERPLOIDY',
                             'BFX_RT_ALLELEREADS', 'BFX_RT_REFERENCEREADS', 'BFX_RT_OTHERREADS', 'BFX_RT_COVERAGE', 'BFX_RT_ALLELERATIO',
                             'BFX_RT_ZYGOSITY']]

            rt_calls = no_fp.loc[no_fp['BFX_RT_CALLSTATUS'].isin(['FOUND', 'NOCALL'])]  # drop NOTFOUND calls
            rt_calls_nf = no_fp.loc[no_fp['BFX_RT_CALLSTATUS'].isin(['NOTFOUND'])]  # drop NOTFOUND calls

            cnv_calls = pd.read_csv(cnv, sep='\t', header=0)

            rt_dfs.append(rt_calls)
            rt_dfs_nf.append(rt_calls_nf)
            agg_cnv_dfs.append(cnv_calls)

    all_rt_calls = pd.concat(rt_dfs, axis=0, sort=False)
    all_nf_rt_calls = pd.concat(rt_dfs_nf, axis=0, sort=False)
    agg_cnv = pd.concat(agg_cnv_dfs, axis=0, sort=False)

    gene = []
    market_name = []
    disease = []
    # add variant metadata
    for ind, row in all_rt_calls.iterrows():
        vsid = row[1]
        gene.append(meta_data.loc[vsid, 'GENE'])
        market_name.append(meta_data.loc[vsid, 'MARKET NAME'])
        disease.append(meta_data.loc[vsid, 'DISEASES'])

    all_rt_calls.rename(columns={'BFX_RT_BFXID': 'SAMPLE_ID',
                           'BFX_RT_ALLELEID': 'ALLELE_ID',
                           'BFX_RT_CALLSTATUS': 'STATUS',
                           'BFX_RT_CALLTYPE': 'TYPE',
                           'BFX_RT_VARIANTQUAL': 'VARIANT_QUALITY',
                           'BFX_RT_CALLQUAL': 'CALL QUALITY',
                           'BFX_RT_ALLELEPLOIDY': 'ALLELE_PLOIDY',
                           'BFX_RT_REFERENCEPLOIDY': 'REFERENCE_PLOIDY',
                           'BFX_RT_OTHERPLOIDY': 'OTHER_PLOIDY',
                           'BFX_RT_ALLELEREADS': 'ALLELE_READS',
                           'BFX_RT_REFERENCEREADS': 'REFERENCE_READS',
                           'BFX_RT_OTHERREADS': 'OTHER_READS'}, inplace=True)

    all_rt_calls.insert(4, 'GENE', gene)
    all_rt_calls.insert(5, 'MARKET_NAME', market_name)
    all_rt_calls.insert(6, 'DISEASE', disease)

    return all_rt_calls, all_nf_rt_calls, agg_cnv


def create_cnv_results(calls):
    """Reformats aggregated CNV calls into desired format. Returns two dataframes of SMN and CNV calls

    Args:
        calls (dataframe): aggregated CNV calls

    Returns: dataframe, dataframe

    """
    # pull all SMN and alpha calls
    smn1_calls = reformat_calls(calls[calls['REGION']=='SMN1'], ['SMN1 Call', 'SMN1 Ploidy', 'SMN1 Quality'])
    smn2_calls = reformat_calls(calls[calls['REGION']=='SMN2'], ['SMN2 Call', 'SMN2 Ploidy', 'SMN2 Quality'])
    hba_calls = reformat_calls(calls[calls['REGION']=='HBA'], ['HBA Call', 'HBA Variant', 'HBA Quality'])

    smn_alpha_only = pd.concat([smn1_calls, smn2_calls, hba_calls], axis=1, join='outer')
    smn_alpha_only.reset_index(inplace=True)

    # other CNVs and PP2, but only positives
    other = calls[~calls['REGION'].isin(['SMN1','SMN2','HBA'])]
    only_positives = other.dropna(subset=['VARIANT_ID'])
    only_positives.drop(labels='EXPECTED_CALL', inplace=True, axis=1)

    if not only_positives.empty: #check if there are other positive CNV variants
        only_positives.loc[:,'STATUS'] = 'FOUND'

    other_to_output = pd.concat([only_positives])

    return smn_alpha_only, other_to_output

def arg_parser():
    parser = argparse.ArgumentParser()

    parser.add_argument('fcid', type=str, help='flowcell ID')
    parser.add_argument('output_path', type=Path, help='output path')
    parser.add_argument('prefix', type=str, help='output filename prefix')

    parser.add_argument('--data', type=Path, help='path to Levitate analysis output directory')
    parser.add_argument('--f_result', type=Path, help='path to flowcell result directory', default=None)
    parser.add_argument('--s_align', type=Path, help='path to sample alignment directory', default=None)
    parser.add_argument('--s_result', type=Path, help='path to sample result directory', default=None)

    parser.add_argument('-tsv', action="store_true", default=False)

    return parser

def cli():
    parser      = arg_parser()
    args        = parser.parse_args()
    output_path = args.output_path
    prefix      = args.prefix

    tsv = args.tsv

    if args.data and all(p is None for p in [args.f_result, args.s_align, args.s_result]):
        flowcell_result_path = args.data / 'flowcell_result'
        sample_alignment_path = args.data / 'sample_alignment'
        sample_result_path = args.data / 'sample_result'
    elif all(p is not None for p in [args.f_result, args.s_align, args.s_result]) and (args.data is None):
        flowcell_result_path = args.f_result
        sample_alignment_path = args.s_align
        sample_result_path = args.s_result
    else:
        raise ValueError('Either provide --f_result, --s_align, and --s_result args or provide --data only')

    fcs = args.fcid.split(',')

    print('Aggregating QC Metrics')
    sample_output, hs_output = aggregate_qc(flowcell_result_path, sample_alignment_path, fcs)

    print('Aggregating RT and PP1')
    rt_pp1, rt_nf, agg_cnv = aggregate_calls(sample_result_path, fcs)

    print('Aggregating SMA, Alphathal, & other positive CNVs')
    smn_alpha, positive_CNV_PP2 = create_cnv_results(agg_cnv)

    print('Outputting to excel')
    with pd.ExcelWriter(f'{output_path}/{prefix}_output.xlsx') as writer:

        sample_output.to_excel(writer, sheet_name='Sample Metrics', index=None)
        hs_output.to_excel(writer, sheet_name='HS Metrics', index=None)

        rt_pp1.to_excel(writer, sheet_name='RT & PP1 Variants', index=None)
        smn_alpha.to_excel(writer, sheet_name='SMN & Alpha', index=None)
        positive_CNV_PP2[positive_CNV_PP2['CALL_TYPE'] == 'GENO'].to_excel(writer, sheet_name='Other CNV', index=None)
        positive_CNV_PP2[positive_CNV_PP2['CALL_TYPE'] == 'PLOIDY'].to_excel(writer, sheet_name='PP2', index=None)

    if tsv:
        print('Outputting TSVs')

        sample_output.to_csv(args.output_path / f'{prefix}_sample_metrics.tsv', sep='\t', index=None)
        hs_output.to_csv(args.output_path / f'{prefix}_hs.tsv', sep='\t', index=None)
        rt_pp1.to_csv(args.output_path / f'{prefix}_rt_pp1.tsv', sep='\t', index=None)
        smn_alpha.to_csv(args.output_path / f'{prefix}_smn_alpha.tsv', sep='\t', index=None)
        positive_CNV_PP2.to_csv(args.output_path / f'{prefix}_other_cnv_pp2.tsv', sep='\t', index=None)
        rt_nf.to_csv(args.output_path / f'{prefix}_rt_pp1_notfound.tsv', sep='\t', index=None)

if __name__ == '__main__':
    cli()
