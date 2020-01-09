"""This script compares variant calls between two Levitate flowcells. Specifically it compares the following:

    * SMN1 and SMN2 calls
    * Readthough/PP1 calls
    * PP2/CNV calls

Arguments:
    - path: tsv file containing Levitate RT/PP1 calls  - flowcell 1
    - path: tsv file containing Levitate RT/PP1 calls  - flowcell 2
    - path: tsv file containing Levitate PP2/CNV calls - flowcell 1
    - path: tsv file containing Levitate PP2/CNV calls - flowcell 2
    - path: tsv file containing Levitate SMN calls     - flowcell 1
    - path: tsv file containing Levitate SMN calls     - flowcell 2
"""

import argparse
import numpy as np

from pathlib import Path

import pandas as pd

def format_compare(compare, fc_id1, fc_id2, sorter, dup_logic=False, cnv=False):

    var_col = 'ALLELE_ID'
    if cnv:
        var_col = 'VARIANT_ID'

    if dup_logic:
        compare['unique'] = compare['PROPS_ID'] + compare[var_col]
        one_fc_call = compare.drop_duplicates(subset='unique', keep=False)['unique'].to_list()
        compare = compare.astype({'Exist' :'str'})
        compare.loc[(compare['unique'].isin(one_fc_call)) & (compare['Exist'] == 'left_only'), 'Exist'] = f'{fc_id1}_only'
        compare.loc[(compare['unique'].isin(one_fc_call)) & (compare['Exist'] == 'right_only'), 'Exist'] = f'{fc_id2}_only'
        compare.drop('unique', axis=1, inplace=True)

    compare.replace('left_only', f'{fc_id1}', inplace=True)
    compare.replace('right_only', f'{fc_id2}', inplace=True)
    compare.replace('both', 'none', inplace=True)
    compare = compare.rename(columns={'Exist': 'DIFFERENCE'})
    compare.sort_values(by=sorter, inplace=True)

    # move PROPS_ID to the front of the dataframe
    cols = list(compare)
    cols.insert(0, cols.pop(cols.index('PROPS_ID')))

    return compare.loc[:, cols]


def compare_rt(f1_rt, f2_rt, fc_id1, fc_id2):

    f1_rt['PROPS_ID'] = f1_rt['SAMPLE_ID'].str.split('_').str[-1]
    f2_rt['PROPS_ID'] = f2_rt['SAMPLE_ID'].str.split('_').str[-1]

    f1_rt = f1_rt.astype({'ALLELE_PLOIDY': 'str', 'REFERENCE_PLOIDY': 'str', 'OTHER_PLOIDY': 'str'})
    f2_rt = f2_rt.astype({'ALLELE_PLOIDY': 'str', 'REFERENCE_PLOIDY': 'str', 'OTHER_PLOIDY': 'str'})

    #drop columns that probably won't match. Not sure about including variant quality here.
    f1_rt.drop(['SAMPLE_ID', 'VARIANT_QUALITY', 'ALLELE_READS',	'REFERENCE_READS', 'OTHER_READS'], axis=1, inplace=True) #drop columns that definitely aren't the same between the two
    f2_rt.drop(['SAMPLE_ID', 'VARIANT_QUALITY', 'ALLELE_READS', 'REFERENCE_READS', 'OTHER_READS'], axis=1, inplace=True)

    compare = pd.merge(f1_rt, f2_rt,
                       on=['PROPS_ID', 'ALLELE_ID', 'STATUS', 'TYPE', 'GENE', 'MARKET_NAME', 'DISEASE', 'CALL_QUALITY',
                           'ALLELE_PLOIDY', 'REFERENCE_PLOIDY', 'OTHER_PLOIDY'], how='outer', indicator='Exist')

    return format_compare(compare, fc_id1, fc_id2,['PROPS_ID', 'ALLELE_ID'], dup_logic=True)

def compare_cnv(f1_cnv, f2_cnv, fc_id1, fc_id2):

    f1_cnv['PROPS_ID'] = f1_cnv['SAMPLE_ID'].str.split('_').str[-1]
    f2_cnv['PROPS_ID'] = f2_cnv['SAMPLE_ID'].str.split('_').str[-1]

    f1_cnv = f1_cnv.astype({'CALL': 'str'})
    f2_cnv = f2_cnv.astype({'CALL': 'str'})

    f1_cnv.drop(['SAMPLE_ID', 'GOF', 'CALL_PLOIDY', 'CALL_QUAL'], axis=1, inplace=True) #drop columns that definitely aren't the same between the two
    f2_cnv.drop(['SAMPLE_ID', 'GOF', 'CALL_PLOIDY', 'CALL_QUAL'], axis=1, inplace=True)

    compare = pd.merge(f1_cnv, f2_cnv,
                       on=['PROPS_ID', 'CALL_TYPE', 'REGION', 'VARIANT_ID', 'CALL', 'STATUS'],
                       how='outer', indicator='Exist')

    return format_compare(compare, fc_id1, fc_id2,['PROPS_ID', 'VARIANT_ID'], dup_logic=True, cnv=True)


def compare_smn_alpha(f1_smn_alpha, f2_smn_alpha, fc_id1, fc_id2):

    f1_smn_alpha['PROPS_ID'] = f1_smn_alpha['SAMPLE_ID'].str.split('_').str[-1]
    f2_smn_alpha['PROPS_ID'] = f2_smn_alpha['SAMPLE_ID'].str.split('_').str[-1]

    f1_smn_alpha = f1_smn_alpha.astype({'SMN1 Call': 'str', 'SMN2 Call': 'str'})
    f2_smn_alpha = f2_smn_alpha.astype({'SMN1 Call': 'str', 'SMN2 Call': 'str'})

    f1_smn_alpha.drop(['SAMPLE_ID', 'SMN1 Ploidy', 'SMN1 Quality', 'SMN2 Ploidy', 'SMN2 Quality', 'HBA Quality'],
                      axis=1, inplace=True)  # drop columns that definitely aren't the same between the two
    f2_smn_alpha.drop(['SAMPLE_ID', 'SMN1 Ploidy', 'SMN1 Quality', 'SMN2 Ploidy', 'SMN2 Quality', 'HBA Quality'],
                      axis=1, inplace=True)

    compare = pd.merge(f1_smn_alpha, f2_smn_alpha,
                       on=['PROPS_ID', 'SMN1 Call', 'SMN2 Call', 'HBA Call', 'HBA Variant'],
                       how='outer', indicator='Exist')

    return format_compare(compare, fc_id1, fc_id2,['PROPS_ID'])


def arg_parser():
    parser = argparse.ArgumentParser()

    parser.add_argument('f1_prefix', type=str, help='flowcell 1 prefix')
    parser.add_argument('f2_prefix', type=str, help='flowcell 2 prefix')
    parser.add_argument('f1_calls', type=Path, help='path to flowcell 1 call tsvs')
    parser.add_argument('f2_calls', type=Path, help='path to flowcell 2 call tsv')

    parser.add_argument('output_path', type=Path, help='output path')

    parser.add_argument('-tsv', action="store_true", default=False)

    return parser

def cli():
    parser = arg_parser()
    args = parser.parse_args()

    output_path = args.output_path

    fc_prefix_1 = args.f1_prefix
    fc_prefix_2 = args.f2_prefix

    rt_file_1 = args.f1_calls / f'{fc_prefix_1}_rt_pp1.tsv'
    rt_file_2 = args.f2_calls / f'{fc_prefix_2}_rt_pp1.tsv'

    cnv_file_1 = args.f1_calls / f'{fc_prefix_1}_other_cnv_pp2.tsv'
    cnv_file_2 = args.f2_calls / f'{fc_prefix_2}_other_cnv_pp2.tsv'

    smn_file_1 = args.f1_calls / f'{fc_prefix_1}_smn_alpha.tsv'
    smn_file_2 = args.f2_calls / f'{fc_prefix_2}_smn_alpha.tsv'

    # fc_id1 = str(args.f1_rt).split('/')[-1].split('_')[0]
    f1_rt = pd.read_csv(rt_file_1, sep='\t', header=0)
    f1_cnv = pd.read_csv(cnv_file_1, sep='\t', header=0)
    f1_smn_alpha = pd.read_csv(smn_file_1, sep='\t', header=0)

    # fc_id2 = str(args.f2_rt).split('/')[-1].split('_')[0]
    f2_rt = pd.read_csv(rt_file_2, sep='\t', header=0)
    f2_cnv = pd.read_csv(cnv_file_2, sep='\t', header=0)
    f2_smn_alpha = pd.read_csv(smn_file_2, sep='\t', header=0)

    rt_compare = compare_rt(f1_rt, f2_rt, fc_prefix_1, fc_prefix_2)
    cnv_compare = compare_cnv(f1_cnv, f2_cnv, fc_prefix_1, fc_prefix_2)
    smn_alpha_compare = compare_smn_alpha(f1_smn_alpha, f2_smn_alpha, fc_prefix_1, fc_prefix_2)

    with pd.ExcelWriter( output_path / 'levitate_comparison_output.xlsx') as writer:

        rt_compare.to_excel(writer, sheet_name='RT & PP1 Comparison', index=None)
        cnv_compare.to_excel(writer, sheet_name='CNV & PP2 Comparison', index=None)
        smn_alpha_compare.to_excel(writer, sheet_name='SMN & Alpha Comparison', index=None)

    if args.tsv:
        rt_compare.to_csv(output_path / 'lev_lev_rt_compare.tsv',sep='\t', index=None)
        cnv_compare.to_csv(output_path / 'lev_lev_cnv_pp2_compare.tsv',sep='\t', index=None)
        smn_alpha_compare.to_csv(output_path / 'lev_lev_smn_alpha.tsv', sep='\t', index=None)

if __name__ == '__main__':
    cli()