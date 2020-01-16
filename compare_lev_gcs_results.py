"""This script compares GCS calls aggregated using M.Hess's script and levitate calls aggregated using
aggregate_results.py. Specifically it compares the following:

    * SMN1 and SMN2 calls
    * Readthrough and CNV calls

An excel workbook named call_comparison.xlsx is generated.

The Classic_Variant_Metadata_v6.txt file is used to determine which variants to compare. This list contains all variants
called by Levitate. Note that alphathal calls are not compared as they were not called using GCS.

Arguments:
    - path: tsv/txt file containing GCS SMN calls
    - path: tsv/txt file containing Levitate SMN calls
    - path: tsv/txt file containing GCS variant calls
    - path: tsv/txt file containing Levitate RT and PP1 variant calls
    - path: tsv/txt file containing Levitate CNV calls
    - path: location to save output file

Example command:

python compare_lev_gcs_results.py /mnt/ruo_rw/rnd/staff/nilanthy.balendra/Tickets/BFX-802_BFXSD-432/gcs_smn_calls.txt
/mnt/ruo_rw/rnd/staff/nilanthy.balendra/Tickets/BFX-802_BFXSD-432/H72T3DRXX_H73JWDRXX_smn_alpha.tsv
/mnt/ruo_rw/rnd/staff/nilanthy.balendra/Tickets/BFX-802_BFXSD-432/gcs_all_calls.txt
/mnt/ruo_rw/rnd/staff/nilanthy.balendra/Tickets/BFX-802_BFXSD-432/H72T3DRXX_H73JWDRXX_rt_pp1.tsv
/mnt/ruo_rw/rnd/staff/nilanthy.balendra/Tickets/BFX-802_BFXSD-432/H72T3DRXX_H73JWDRXX_other_cnv_pp2.tsv
/mnt/ruo_rw/rnd/staff/nilanthy.balendra/Tickets/BFX-802_BFXSD-432
"""

import argparse
import numpy as np

from pathlib import Path

import pandas as pd


def compare_smn(gcs_smn, lev_smn):
    """Compares SMN calls from GCS and Levitate and outputs their status as EQUAL or NOT EQUAL
    Args:
        gcs_smn(dataframe): gcs SMN calls
        lev_smn(dataframe): levitate SMN calls

    Returns:
        dataframe with GCS and Levitate SMN calls and columns indicating their equivalence

    """
    gcs_smn.replace('>= 4', '4', inplace=True)

    lev_smn['PROPS_ID'] = lev_smn['SAMPLE_ID'].str.split('_').str[-1]

    gcs_smn_1 = gcs_smn.loc[gcs_smn['Variant'] == 'SMN1'].set_index('SampleID')
    gcs_smn_1['GCS_SMN1_CALL'] = gcs_smn_1['SMN FINAL']
    gcs_smn_1['GCS_SMN1_RAW'] = gcs_smn_1['smn_raw']

    gcs_smn_2 = gcs_smn.loc[gcs_smn['Variant'] == 'SMN2'].set_index('SampleID')
    gcs_smn_2['GCS_SMN2_CALL'] = gcs_smn_2['SMN FINAL']
    gcs_smn_2['GCS_SMN2_RAW'] = gcs_smn_2['smn_raw']

    all_gcs_smn = gcs_smn_1[['GCS_SMN1_CALL', 'GCS_SMN1_RAW']].join(gcs_smn_2[['GCS_SMN2_CALL', 'GCS_SMN2_RAW']], sort=False, how='left')

    full = lev_smn.join(all_gcs_smn, sort='False', how='left', on='PROPS_ID')

    full['SMN1_COMPARE'] = np.where(full['GCS_SMN1_CALL'] == full['SMN1 Call'].astype(str), 'EQUAL', 'NOT_EQUAL')
    full['SMN2_COMPARE'] = np.where(full['GCS_SMN2_CALL'] == full['SMN2 Call'].astype(str), 'EQUAL', 'NOT_EQUAL')

    return full


def compare_variants(gcs_vars, lev_rt, lev_cnv):
    """Compares calls between GCS and Levitate for each of the given dataframes. For the GCS calls, variants that are
    not called by the Levitate Classic product are not compared or outputted.

    Args:
        gcs_vars: GCS variant calls
        lev_rt: Levitate readthrough and PP1 calls
        lev_cnv: Levitate CNV variant calls

    Returns:
        The above dataframes with additional comparison information.

    """

    lev_rt['PROPS_ID'] = lev_rt['SAMPLE_ID'].str.split('_').str[-1]
    lev_cnv['PROPS_ID'] = lev_cnv['SAMPLE_ID'].str.split('_').str[-1]

    # drop variants that aren't called in levitate
    lev_meta_data = pd.read_csv('Classic_Variant_Metadata_v6.txt', sep='\t', header=0)
    gcs_vars = gcs_vars.loc[(gcs_vars['VariantID'].isin(lev_meta_data['VSID'])) & (gcs_vars['CallType'] == 'Positive')]

    # check if gcs calls are called by levitate
    check = []
    ploidy = []
    for i, row in gcs_vars.iterrows():
        # try to find matching levitate call
        lev_rt_call = lev_rt.loc[(lev_rt['PROPS_ID'].astype(str) == str(row['SampleID'])) & (lev_rt['ALLELE_ID'] == row['VariantID'])]
        lev_cnv_call = lev_cnv.loc[(lev_rt['PROPS_ID'].astype(str) == str(row['SampleID'])) & (lev_cnv['VARIANT_ID'] == row['VariantID'])]

        if lev_rt_call.empty & lev_cnv_call.empty:
            ploidy.append(np.nan)
            check.append('NOTFOUND in Levitate')
        elif lev_rt_call.empty: #then it is a CNV
            ploidy.append(lev_cnv_call['CALL'].item())
            if str(row['AlleleCount']) == str(lev_cnv_call['CALL'].item()):
                check.append('FOUND')
            else:
                check.append('FOUND, ploidy mismatch')
        else: #it is an rt
            ploidy.append(lev_rt_call['ALLELE_PLOIDY'].item())
            if str(row['AlleleCount']) == str(lev_rt_call['ALLELE_PLOIDY'].item()):
                check.append('FOUND')
            else:
                check.append('FOUND, ploidy mismatch')
    gcs_vars['Levitate_Ploidy'] = ploidy
    gcs_vars['COMPARISON'] = check

    # check if Levitate RT/PP1 variants are
    ploidy = []
    check = []
    for i, row in lev_rt.iterrows():
        gcs_call = gcs_vars.loc[(gcs_vars['SampleID'].astype(str) == row['PROPS_ID']) & (gcs_vars['VariantID'] == row['ALLELE_ID'])]

        if gcs_call.empty:
            ploidy.append(np.nan)
            check.append('NOTFOUND in GCS')
        else:
            ploidy.append(gcs_call['AlleleCount'].item())
            if str(row[9]) == str(gcs_call['AlleleCount'].item()):
                check.append('FOUND')
            else:
                check.append('FOUND, ploidy mismatch')

    lev_rt['GCS_Ploidy'] = ploidy
    lev_rt['COMPARISON'] = check

    ploidy = []
    check = []
    for i, row in lev_cnv.iterrows():
        gcs_call = gcs_vars.loc[(gcs_vars['SampleID'].astype(str) == row['PROPS_ID']) & (gcs_vars['VariantID'] == row['VARIANT_ID'])]

        if gcs_call.empty:
            ploidy.append(np.nan)
            check.append('NOTFOUND in GCS')
        else:
            ploidy.append(gcs_call['AlleleCount'].item())
            if row['CALL'] == str(gcs_call['AlleleCount'].item()):
                check.append('FOUND')
            else:
                check.append('FOUND, ploidy mismatch')

    lev_cnv['GCS_Ploidy'] = ploidy
    lev_cnv['COMPARISON'] = check

    return gcs_vars, lev_rt, lev_cnv


def arg_parser():
    parser = argparse.ArgumentParser()

    parser.add_argument('gcs_smn', type=Path, help='path to gcs SMN calls')
    parser.add_argument('lev_smn', type=Path, help='path to lev SMN calls')

    parser.add_argument('gcs_var_calls', type=Path, help='path to gcs variant calls')
    parser.add_argument('lev_rt_calls', type=Path, help='path to lev readthrough calls')
    parser.add_argument('lev_cnv_calls', type=Path, help='path to lev CNV/PP2 calls')

    parser.add_argument('output_path', type=Path, help='output path')

    return parser


def cli():
    parser = arg_parser()
    args = parser.parse_args()

    output_path = args.output_path

    gcs_smn = pd.read_csv(args.gcs_smn, sep='\t', header=0)
    lev_smn = pd.read_csv(args.lev_smn, sep='\t', header=0)

    gcs_vars = pd.read_csv(args.gcs_var_calls, sep='\t', header=0)
    lev_rt   = pd.read_csv(args.lev_rt_calls, sep='\t', header=0)
    lev_cnv  = pd.read_csv(args.lev_cnv_calls, sep='\t', header=0)

    comparison_smn = compare_smn(gcs_smn, lev_smn[['SAMPLE_ID', 'SMN1 Call', 'SMN1 Ploidy', 'SMN1 Quality', 'SMN2 Call', 'SMN2 Ploidy', 'SMN2 Quality']])

    comparison_gcs, comparison_lev, comparison_lev_cnv = compare_variants(gcs_vars, lev_rt, lev_cnv)

    with pd.ExcelWriter(f'{output_path}/call_comparison.xlsx') as writer:
        comparison_gcs.to_excel(writer, sheet_name='GCS Call Comparison', index=None)
        comparison_lev.to_excel(writer, sheet_name='Levitate Readthrough Comparison', index=None)
        comparison_lev_cnv.to_excel(writer, sheet_name='Levitate CNV_PP2 Comparison', index=None)
        comparison_smn.to_excel(writer, sheet_name='SMN Comparison', index=None)

    # comparison_smn.to_csv(output_path / 'comparison_smn.tsv', sep='\t')
    # comparison_gcs.to_csv(output_path / 'tmp_gcs.tsv', sep='\t', index=None)
    # comparison_lev.to_csv(output_path / 'tmp_lev.tsv', sep='\t', index=None)
    # comparison_lev_cnv.to_csv(output_path / 'tmp_lev_cnv.tsv', sep='\t', index=None)

if __name__ == '__main__':
    cli()
