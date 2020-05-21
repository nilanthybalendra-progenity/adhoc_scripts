import argparse
from pathlib import Path
from subprocess                   import PIPE, Popen


def cluster(args):

    fc_path = args.fc_path
    out_dir = args.out_dir

    fc_list = [line.rstrip('\n') for line in open(fc_path)]

    sbatch_cmd = ['sbatch']

    for i, fc in enumerate(fc_list):
        print(i)
        cmd = f'python3 agg_counts.py {fc} {out_dir}'
        proc = Popen(sbatch_cmd, stdin=PIPE, stdout=PIPE, encoding='utf-8')
        proc.communicate('#!/bin/bash\n\n' + cmd + '\n')


def arg_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('fc_path', type=Path, help='path to fc list')
    parser.add_argument('out_dir', type=Path, help='path to output directory')

    return parser

def cli():
    parser = arg_parser()
    args = parser.parse_args()

    cluster(args)


if __name__ == '__main__':
    cli()
