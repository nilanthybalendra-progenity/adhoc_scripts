"""This script pulls together data from flowcells of interest for review

Arguments:
    - string: flowcells to pull data from, one value or comma separated list
    - path: sample_result directory
    - path: output directory
    - string: output file prefix

Optional:
    --aligns, default value is None, if a path is given to sample alignment directory, qc values are aggregated
    -tsv, default value is False, if true, the following individual tsv files are outputted
        - {prefix}_hs.tsv
        - {prefix}_other_cnv_pp2.tsv
        - {prefix}_qc.tsv (if -aligns is not None)
        - {prefix}_smna_alpha.tsv
        - {prefix}_rt_pp1.tsv

Example command:
python aggregate_results.py H737HDRXX,HCYTJDRXX /mnt/bfx_analysis_ro/analysis_data/progenity_workflow_levitate/sample_result/
/mnt/ruo_rw/rnd/staff/nilanthy.balendra/ BLAH --aligns /mnt/bfx_analysis_ro/analysis_data/progenity_workflow_levitate/sample_alignment/ -tsv

The following file is always created and contains all contents in the above tsv files on separate sheets:
    - {prefix}_output.xlsx

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
    """sample alignment and sample result directories contain directories for each sample appended by "v##". The latest
    output for that sample has the greatest ## value. This function puts together a list of the latest run directories"""

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


def aggregate_qc(main_path, fc_id):
    ''' Given a flowcell ID or list of them, aggregate the following files together required qc metrics
    - hs-metrics.txt
    - gaps_classic.txt: CDS, NEAR CDS, DEEP, FREEMIX gaps for PP1 and PP2
    - calls.tsv (from CNV): get GOF
    '''

    dirs = get_file_list(main_path, fc_id)

    qc = []
    hs = []

    for d in dirs:
        hs_metrics = main_path / d / 'metrics' / 'hs_metrics.tsv'
        gaps_rt = main_path / d / 'metrics' / 'gaps_rt_classic_summary.tsv'
        gaps_pp = main_path / d / 'metrics' / 'gaps_pp_classic_summary.tsv'
        crud = main_path/ d / 'crud' / 'contamination_metrics.tsv'
        fit = main_path/ d / 'fit' / 'fit.tsv'

        if os.path.isfile(hs_metrics) and os.path.isfile(gaps_rt) and os.path.isfile(crud) and os.path.isfile(fit):

            # hs_metrics
            hs_lines = open(hs_metrics).readlines()
            rest = hs_lines[7].split('\t')
            headings = hs_lines[6].split('\t')

            fold_80 = rest[43]
            pf_unique_reads = rest[25]

            # classic gaps
            gaps_rt = open(gaps_rt).readlines()
            gaps_pp = open(gaps_pp).readlines()

            # contamination
            contamination = pd.read_csv(crud, sep='\t', header=0)
            free_mix = contamination['FREEMIX'][0]

            #gof
            fit = pd.read_csv(fit, sep='\t', header=0)
            gof = fit['GOF'][0]

            # aggregate
            qc.append([d[:-4], gof, pf_unique_reads, fold_80] + gaps_rt[1][:-1].split('\t')[1:-1] + gaps_pp[1][:-1].split('\t')[1:-1] + [free_mix])
            hs.append([d[:-4]] + rest)


    qc_output = pd.DataFrame(data=qc, columns=['SAMPLE_ID', 'GOF', 'UNIQUE_ALIGNED_READS', 'FOLD80_BASE_PENALTY', 'RT_CDS_GAPS', 'RT_NEAR_CDS_GAPS', 'RT_DEEP_GAPS', 'PP1_CDS_GAPS', 'PP1_NEAR_CDS_GAPS', 'PP1_DEEP_GAPS', 'FREEMIX'])
    header = ['SAMPLE_ID'] + headings
    hs_metrics = pd.DataFrame(data=hs, columns=header)

    return qc_output, hs_metrics

def arg_parser():
    parser = argparse.ArgumentParser()

    parser.add_argument('fcid', type=str, help='flowcell ID')
    parser.add_argument('sample_result_path', type=Path, help='path to sample result directory')
    parser.add_argument('output_path', type=Path, help='output path')
    parser.add_argument('prefix', type=Path, help='output filename prefix')
    parser.add_argument('-tsv', action="store_true", default=False)
    parser.add_argument('--aligns', type=Path, help='path to sample alignment directory', default=None)

    return parser

def cli():
    parser      = arg_parser()
    args        = parser.parse_args()
    output_path = args.output_path
    prefix      = args.prefix

    tsv = args.tsv
    sample_alignment_path = args.aligns

    fcs = args.fcid.split(',')

    rt_pp1           = pd.DataFrame()
    smn_alpha        = pd.DataFrame()
    positive_CNV_PP2 = pd.DataFrame()
    qc_output        = pd.DataFrame()
    hs_output        = pd.DataFrame()

    for fc in fcs:
        print(fc)

        if sample_alignment_path:
            print('Aggregating QC Metrics')
            f_qc_output, f_hs_output = aggregate_qc(Path(sample_alignment_path), fc)
            qc_output = pd.concat([qc_output, f_qc_output])
            hs_output = pd.concat([hs_output, f_hs_output])

        print('Aggregating RT and PP1')
        f_rt_pp1, agg_cnv = aggregate_rt(args.sample_result_path, fc)

        print('Aggregating SMA, Alphathal, & other positive CNVs')
        f_smn_alpha, f_positive_CNV_PP2 = agg_cnv_results(agg_cnv)

        rt_pp1 = pd.concat([rt_pp1, f_rt_pp1])
        smn_alpha = pd.concat([smn_alpha, f_smn_alpha])
        positive_CNV_PP2 = pd.concat([positive_CNV_PP2, f_positive_CNV_PP2])

    print('Outputting to excel')
    with pd.ExcelWriter(f'{output_path}/{prefix}_output.xlsx') as writer:

        if sample_alignment_path:
            qc_output.to_excel(writer, sheet_name='Sample QC Summary', index=None)
            hs_output.to_excel(writer, sheet_name='HS Metrics', index=None)

        rt_pp1.to_excel(writer, sheet_name='RT & PP1 Variants', index=None)
        smn_alpha.to_excel(writer, sheet_name='SMN & Alpha', index=None)
        positive_CNV_PP2[positive_CNV_PP2['CALL_TYPE'] == 'GENO'].to_excel(writer, sheet_name='Other CNV', index=None)
        positive_CNV_PP2[positive_CNV_PP2['CALL_TYPE'] == 'PLOIDY'].to_excel(writer, sheet_name='PP2', index=None)

    if tsv:
        print('Outputting TSVs')

        if sample_alignment_path:
            qc_output.to_csv(args.output_path / f'{prefix}_qc.tsv', sep='\t', index=None)

        hs_output.to_csv(args.output_path / f'{prefix}_hs.tsv', sep='\t', index=None)
        rt_pp1.to_csv(args.output_path / f'{prefix}_rt_pp1.tsv', sep='\t', index=None)
        smn_alpha.to_csv(args.output_path / f'{prefix}_smn_alpha.tsv', sep='\t', index=None)
        positive_CNV_PP2.to_csv(args.output_path / f'{prefix}_other_cnv_pp2.tsv', sep='\t', index=None)


if __name__ == '__main__':
    cli()
