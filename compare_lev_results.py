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


def compare_rt(f1_rt, f2_rt, fc1, fc2):

    f1_rt['PROPS_ID'] = f1_rt['SAMPLE_ID'].str.split('_').str[-1]
    f2_rt['PROPS_ID'] = f2_rt['SAMPLE_ID'].str.split('_').str[-1]

    #drop columns that probably won't match. Not sure about including variant quality here.
    f1_rt.drop(['SAMPLE_ID', 'VARIANT_QUALITY', 'ALLELE_READS',	'REFERENCE_READS', 'OTHER_READS'], axis=1, inplace=True)
    f2_rt.drop(['SAMPLE_ID', 'VARIANT_QUALITY', 'ALLELE_READS', 'REFERENCE_READS', 'OTHER_READS'], axis=1, inplace=True)

    compare = pd.merge(f1_rt, f2_rt, on=['PROPS_ID', 'ALLELE_ID', 'STATUS', 'TYPE', 'GENE', 'MARKET_NAME', 'DISEASE', 'CALL_QUALITY', 'ALLELE_PLOIDY', 'REFERENCE_PLOIDY', 'OTHER_PLOIDY'], how='outer', indicator='Exist')

    return compare


def arg_parser():
    parser = argparse.ArgumentParser()

    parser.add_argument('f1_rt', type=Path, help='path to flowcell 1 readthrough/PP1 calls')
    parser.add_argument('f2_rt', type=Path, help='path to flowcell 2 readthough/PP1 calls')
    # parser.add_argument('f1_cnv', type=Path, help='path to flowcell 1 CNV/PP2 calls')
    # parser.add_argument('f2_cnv', type=Path, help='path to flowcell 2 CNV/PP2 calls')
    # parser.add_argument('f1_smn', type=Path, help='path to flowcell 1 SMN calls')
    # parser.add_argument('f2_smn', type=Path, help='path to flowcell 2 SMN calls')

    #parser.add_argument('output_path', type=Path, help='output path')

    return parser

def cli():
    parser = arg_parser()
    args = parser.parse_args()

    f1_rt = pd.read_csv(args.f1_rt, sep='\t', header=0)
    #f1_cnv = pd.read_csv(args.f1_cnv, sep='\t', header=0)
    #f1_smn = pd.read_csv(args.f1_smn, sep='\t', header=0)

    f2_rt = pd.read_csv(args.f2_rt, sep='\t', header=0)
    #f2_cnv = pd.read_csv(args.f2_cnv, sep='\t', header=0)
    #f2_smn = pd.read_csv(args.f2_smn, sep='\t', header=0)

    rt_compare = compare_rt(f1_rt, f2_rt, 'flowcell_1', 'flowcell_2')

    rt_compare.to_csv('blarf.tsv',sep='\t', index=None)


if __name__ == '__main__':
    cli()