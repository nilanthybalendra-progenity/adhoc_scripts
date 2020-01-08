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
python aggregate_results.py H737HDRXX,HCYTJDRXX  /mnt/ruo_rw/rnd/staff/nilanthy.balendra/ my_prefix
--data /mnt/bfx_analysis_ro/analysis_data/progenity_workflow_levitate/

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

def reformat_calls(calls, header):
    """reformat the CNV output files according to variant region (SMN or ALPHA)"""

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


def get_file_list(main_path, fc):
    """flowcell result, sample alignment and sample result directories contain directories for each flowcell/sample
    appended by "v##". The latest output has the greatest ## value. This function puts together a list of the latest run
    directories
    """

    files = pd.DataFrame()
    files['BFX_ID'] = [s for s in os.listdir(main_path) if s.split('_')[0] in fc]
    files['SAMPLE_ID'] = files['BFX_ID'].str[:-4] # remove version tag

    files.sort_values(by='BFX_ID', inplace=True)
    files.drop_duplicates(subset='SAMPLE_ID', keep='last', inplace=True)

    return files['BFX_ID'].tolist()


def aggregate_rt(main_path, fc_id):
    """Given path to rt and cnv results directories, flowcell ids and version, aggregates calls"""

    meta_data = pd.read_csv('Classic_Variant_Metadata_v6.txt', sep='\t', header=0)
    meta_data.set_index('VSID', inplace=True)

    dirs = get_file_list(main_path, fc_id)

    all_vgraph = pd.DataFrame()
    agg_cnv = pd.DataFrame()

    for d in dirs:
        rt = main_path / d / 'readthrough' / 'calls_rt.tsv'
        pp1 = main_path / d / 'readthrough' / 'calls_pp1.tsv'
        cnv = main_path / d / 'cnv' / 'calls.tsv'

        if os.path.isfile(rt) and os.path.isfile(pp1)and os.path.isfile(cnv):

            rt_output = pd.read_csv(rt, sep='\t', header=0)
            no_fp = rt_output[~rt_output['ALLELE_ID'].isin(fp)]
            all_rt = no_fp[no_fp['STATUS'].isin(['FOUND','NOCALL'])]
            all_rt.insert(3, 'TYPE', 'READTHROUGH')

            pp1_output = pd.read_csv(pp1, sep='\t', header=0)
            found_pp1  = pp1_output[pp1_output['STATUS'].isin(['FOUND','NOCALL'])]
            all_pp1 = found_pp1[~found_pp1['ALLELE_ID'].str.contains('CYP21A2')]
            all_pp1.insert(3,'TYPE', 'PP1')

            cnv_calls = pd.read_csv(cnv, sep='\t', header=0)

            all_vgraph = pd.concat([all_vgraph, all_rt, all_pp1], axis=0, sort=False)
            agg_cnv = pd.concat([agg_cnv, cnv_calls], axis=0, sort=False)

    gene = []
    market_name = []
    disease = []

    for ind, row in all_vgraph.iterrows():
        vsid = row[1]
        gene.append(meta_data.loc[vsid, 'GENE'])
        market_name.append(meta_data.loc[vsid, 'MARKET NAME'])
        disease.append(meta_data.loc[vsid, 'DISEASES'])

    all_vgraph.insert(4, 'GENE', gene)
    all_vgraph.insert(5, 'MARKET_NAME', market_name)
    all_vgraph.insert(6, 'DISEASE', disease)

    return all_vgraph, agg_cnv


def agg_cnv_results(calls):
    """Reformats aggregated CNV calls into desired format"""

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
    only_positives.loc[:,'STATUS'] = 'FOUND'

    other_to_output = pd.concat([only_positives])

    return smn_alpha_only, other_to_output


def aggregate_qc(flowcell_result_path, aligns_path, fc_id):
    ''' Given a flowcell ID or list of them, aggregate the following files together required qc metrics
    - hs-metrics.txt
    - gaps_classic.txt: CDS, NEAR CDS, DEEP, FREEMIX gaps for PP1 and PP2
    - calls.tsv (from CNV): get GOF
    '''

    sample_dirs = get_file_list(aligns_path, fc_id)
    f_dir = get_file_list(flowcell_result_path, fc_id)
    hs = []

    for d in sample_dirs:
        hs_metrics_path = aligns_path / d / 'metrics' / 'hs_metrics.tsv'

        if os.path.isfile(hs_metrics_path):
            hs_lines = open(hs_metrics_path).readlines()
            rest = hs_lines[7].split('\t')
            headings = hs_lines[6].split('\t')
            hs.append([d[:-4]] + rest)

    header = ['SAMPLE_ID'] + headings
    hs_metrics = pd.DataFrame(data=hs, columns=header)

    sample_metrics_path = flowcell_result_path / f_dir[0] / 'sample_metrics.tsv'

    if os.path.isfile(sample_metrics_path):
        sample_metrics = pd.read_csv(sample_metrics_path, sep='\t', header=0)
        cols = list(sample_metrics)
        cols = sorted(cols)
        cols.insert(0, cols.pop(cols.index('LIS_REQUEST_BFXID')))
        sample_metrics = sample_metrics.loc[:, cols]

    return sample_metrics, hs_metrics

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

    rt_pp1           = pd.DataFrame()
    smn_alpha        = pd.DataFrame()
    positive_CNV_PP2 = pd.DataFrame()
    sample_output    = pd.DataFrame()
    hs_output        = pd.DataFrame()

    for fc in fcs:
        print(fc)

        if sample_alignment_path:
            print('Aggregating QC Metrics')
            f_sample_output, f_hs_output = aggregate_qc(flowcell_result_path, sample_alignment_path, fc)
            sample_output = pd.concat([sample_output, f_sample_output])
            hs_output = pd.concat([hs_output, f_hs_output])

        print('Aggregating RT and PP1')
        f_rt_pp1, agg_cnv = aggregate_rt(sample_result_path, fc)

        print('Aggregating SMA, Alphathal, & other positive CNVs')
        f_smn_alpha, f_positive_CNV_PP2 = agg_cnv_results(agg_cnv)

        rt_pp1 = pd.concat([rt_pp1, f_rt_pp1])
        smn_alpha = pd.concat([smn_alpha, f_smn_alpha])
        positive_CNV_PP2 = pd.concat([positive_CNV_PP2, f_positive_CNV_PP2])

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


if __name__ == '__main__':
    cli()
