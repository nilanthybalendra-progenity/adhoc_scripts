import argparse
import os
from pathlib import Path

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import plotnine as p9
from plotnine import *

from scipy import stats


def standard_scatter(data_set,title, xaxis, yaxis, xlab, ylab, out_file):

    # or with a regression line
    slope, intercept, r_value, p_value, std_err = stats.linregress(data_set[xaxis], data_set[yaxis])

    fig = plt.figure()

    if intercept < 0:
        ax = sns.regplot(x=xaxis, y=yaxis, data=data_set, color='b', ci=None,
                         scatter_kws={'s': 8}, line_kws={'color': 'red',
                                                         'label': f'y={slope:.3f}x - {abs(intercept):.3f} \n R-Sq= {r_value**2:.3f}'})
    else:
        ax = sns.regplot(x=xaxis, y=yaxis, data=data_set, color='b', ci=None,
                         scatter_kws={'s': 8}, line_kws={'color': 'red',
                                                         'label': f'y={slope:.3f}x + {intercept:.3f} \n R-Sq= {r_value**2:.3f}'})

    ax.legend()
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    plt.title(title)

    fig.savefig(out_file, bbox_inches='tight', dpi=250)


def standard_scatter2(data_set1, title, xaxis, yaxis, xlab, ylab, out_file):

    p = ggplot(data_set1, aes(xaxis, yaxis, color='Model')) + geom_point(alpha=0.1) \
        + labs(x=xlab, y=ylab, title=title) + theme_bw() + stat_smooth(method='lm', se=False)

    p.save(out_file, format='png', dpi=500)



def make_scatter(call_data_path, call_data_path2, output_path):

    data1 = pd.read_csv(call_data_path, sep='\t', header=0)
    data2 = pd.read_csv(call_data_path2, sep='\t', header=0)

    data1['Model'] = 'Model 3.1'
    data2['Model'] = 'Model 4o'

    data = pd.concat([data1, data2], join='inner', axis=0)
    print(data.columns)

    m = '3.1'
    chr_num = ['13', '18', '21']

    for num in chr_num:
        standard_scatter2(data[data[f'CHR{num}_CALL'] != 'FETAL EUPLOIDY'],
                         f'Model{m} Trisomy {num}: T-Value vs. SNP Fetal Fraction',
                         'SNP_FETAL_PCT',
                         f'CHR{num}_TVALUE',
                         'SNP Fetal Fraction (%)',
                         f'chr{num} T-Value',
                         output_path / f'T{num}_Model{m}_t_event.png')


def nipt_null_histogram(call_file_path, title, output_file, plot_x=False):

    data = pd.read_csv(call_file_path, sep='\t', header=0)

    fig = plt.figure()
    x = np.linspace(-5, 5, 1000)
    norm = stats.norm().pdf(x)

    sns.lineplot(x, norm, color='k', label='Standard Normal')

    sns.kdeplot(data.loc[data['CHR13_CALL'] == 'FETAL EUPLOIDY', 'CHR13_TVALUE'], color='r', label='CHR13 T-Value')
    sns.kdeplot(data.loc[data['CHR18_CALL'] == 'FETAL EUPLOIDY', 'CHR18_TVALUE'], color='g', label='CHR18 T-Value')
    sns.kdeplot(data.loc[data['CHR21_CALL'] == 'FETAL EUPLOIDY', 'CHR21_TVALUE'], color='b', label='CHR21 T-Value')

    if plot_x:
        sns.kdeplot(data.loc[data['CHRXY_CALL'] == 'FETAL EUPLOIDY, FEMALE', 'CHRX_TVALUE'], color='orange', label='CHRX T-Value')

    plt.xlabel('T-Value')
    plt.title(title)
    plt.legend(frameon=False)
    fig.show()
    fig.savefig(output_file, bbox_inches='tight', dpi=250)


def arg_parser():
    parser = argparse.ArgumentParser()
    commands = parser.add_subparsers(title='Plots', dest='plot', metavar='')

    null_dist = commands.add_parser('null_dist', help='plot null distribution of t values')
    null_dist.add_argument('call_file', type=Path, help='Path to aggregated call results')
    null_dist.add_argument('plot_title', type=str, help='Plot title')
    null_dist.add_argument('output_path', type=Path, help='Plot output path')
    null_dist.add_argument('-x', action='store_true', default=False, help='Plots chrx t-value distribution of female fetuses')

    scatter = commands.add_parser('scatter', help='plot scatter plots')
    scatter.add_argument('call_file', type=Path, help='Path to aggregated call results')
    scatter.add_argument('call_file2', type=Path, help='Path to aggregated call results')
    scatter.add_argument('output_path', type=Path, help='Plot output path')

    return parser


def cli():

    parser = arg_parser()

    args = parser.parse_args()

    if args.plot == 'null_dist':
        nipt_null_histogram(os.path.abspath(args.call_file), args.plot_title, os.path.abspath(args.output_path), args.x)

    elif args.plot == 'scatter':
        make_scatter(os.path.abspath(args.call_file), os.path.abspath(args.call_file2), Path(os.path.abspath(args.output_path)))


if __name__ == '__main__':
    cli()










    #
    # for num in chr_num:
    #     standard_scatter(data,
    #                      f'Model{m}: chr{num} T-Value vs. SNP Fetal Fraction',
    #                      'SNP_FETAL_PCT',
    #                      f'CHR{num}_TVALUE',
    #                      'SNP Fetal Fraction (%)',
    #                      f'chr{num} T-Value',
    #                      output_path / f'chr{num}_Model{m}_t_snp.png')
    #
    #     standard_scatter(data[data[f'CHR{num}_CALL'] == 'FETAL TRISOMY'],
    #                      f'Model{m} Trisomy {num}: T-Value vs. chr{num} Fetal Fraction',
    #                      f'CHR{num}_FETAL_PCT',
    #                      f'CHR{num}_TVALUE',
    #                      f'chr{num} Fetal Fraction (%)',
    #                      f'chr{num} T-Value',
    #                      output_path / f'T{num}_Model{m}_t_event.png')
    #
    #     standard_scatter(data[data[f'CHR{num}_CALL'] == 'FETAL TRISOMY'],
    #                      f'Model{m} Trisomy {num}: chr{num} Fetal Fraction vs. SNP Fetal Fraction',
    #                      'SNP_FETAL_PCT',
    #                      f'CHR{num}_FETAL_PCT',
    #                      'SNP Fetal Fraction (%)',
    #                      f'chr{num} Fetal Fraction (%)',
    #                      output_path / f'T{num}_Model{m}_event_snp.png')