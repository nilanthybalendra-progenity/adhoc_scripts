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


def standard_scatter2(data_set1, title, xaxis, yaxis, xlab, ylab, out_file, xlim=None, ylim=None):

    p = ggplot(data_set1, aes(xaxis, yaxis, color='Model')) + geom_point(alpha=0.1) \
        + labs(x=xlab, y=ylab, title=title) + theme_bw() + stat_smooth(method='lm', se=False)

    if xlim:
        p += p9.scales.xlim(xlim[0], xlim[1])

    if ylim:
        p += p9.scales.ylim(ylim[0], ylim[1])

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

    # for num in chr_num:
    #     standard_scatter2(data[data[f'CHR{num}_CALL'] == 'FETAL TRISOMY'], #what about monosomy?
    #                      f'Model{m} Trisomy {num}: T-Value vs. SNP Fetal Fraction',
    #                      'SNP_FETAL_PCT',
    #                      f'CHR{num}_TVALUE',
    #                      'SNP Fetal Fraction (%)',
    #                      f'chr{num} T-Value',
    #                      output_path / f'T{num}_Model{m}_t_event.png', xlim=(0,10), ylim=(0,15))

    standard_scatter2(data[~data[f'CHRXY_CALL'].isin(['FETAL EUPLOIDY, FEMALE', 'FETAL EUPLOIDY, MALE'])], #what about monosomy?
                         f'chrXY Aneuploidy: T-Value vs. SNP Fetal Fraction',
                         'SNP_FETAL_PCT',
                         f'CHRX_TVALUE',
                         'SNP Fetal Fraction (%)',
                         f'chrX T-Value',
                         output_path / f'X_Model{m}_t_event.png', xlim=(0,10), ylim=(0,15))


def nipt_null_x(call_1_path, call_2_path, title, output_file, male=False):
    data_1 = pd.read_csv(call_1_path, sep='\t', header=0)
    data_2 = pd.read_csv(call_2_path, sep='\t', header=0)

    # if call_3_path:
    #      data_3 = pd.read_csv(call_3_path, sep='\t', header=0)
    #      mod_3 = "Model 5, 1mm"

    mod_1 = "Model 3.1"
    mod_2 = "Model 5"

    if male:
        sort_crit = 'FETAL EUPLOIDY, MALE'
        x = np.linspace(-50, 50, 1000)
    else:
        sort_crit = 'FETAL EUPLOIDY, FEMALE'
        x = np.linspace(-10, 10, 1000)

    norm = stats.norm().pdf(x)

    fig = plt.figure()
    sns.lineplot(x, norm, color='k', label='Standard Normal')
    sns.kdeplot(data_1.loc[data_1['CHRXY_CALL'] == sort_crit, 'CHRX_TVALUE'], color='red',
                label=f'{mod_1}')
    print(f'{mod_1} CHRX_TVALUE: {data_1.loc[data_1["CHRXY_CALL"] == sort_crit, "CHRX_TVALUE"].mean()}')

    sns.kdeplot(data_2.loc[data_2['CHRXY_CALL'] == sort_crit, 'CHRX_TVALUE'], color='blue',
                label=f'{mod_2}')
    print(f'{mod_2} CHRX_TVALUE: {data_2.loc[data_2["CHRXY_CALL"] == sort_crit, "CHRX_TVALUE"].mean()}')
    
    # if call_3_path:
    #      sns.kdeplot(data_3.loc[data_3['CHRXY_CALL'] == sort_crit, 'CHRX_TVALUE'], color='blue',
    #             label=f'{mod_3}')
    # print(f'{mod_3} CHRX_TVALUE: {data_3.loc[data_3["CHRXY_CALL"] == sort_crit, "CHRX_TVALUE"].mean()}')
    plt.xlabel('T-Value')
    plt.title(title)
    plt.vlines(0, 0, 0.42, colors='k', linestyles='dotted')
    #plt.ylim(0, 0.2)
    plt.ylim(min(norm), 0.42)

    if male:
        plt.xlim(-50, 5)

    plt.legend(frameon=False)

    fig.savefig(output_file, bbox_inches='tight', dpi=250)


def chrY_scatter(call_file_path, manifest_file_path, out_file):


    data = pd.read_csv(call_file_path, sep='\t', header=0)
    manifest = pd.read_csv(manifest_file_path, sep='\t', header=0)

    new_data = data.join(manifest[['WELL', 'SAMPLE_ID', 'SAMPLE_TYPE', 'COMPANY']].set_index('SAMPLE_ID'), sort=False, how='left', on='SAMPLE_ID')

    new_data.loc[((new_data['WELL'] == 'E12') & (new_data['SAMPLE_TYPE'] == 'Control')), 'FETAL_SEX'] = 'NON PREGNANT MALE'

    new_data.dropna(subset=['FETAL_SEX'], inplace=True)
  # new_data = new_data.loc[new_data["COMPANY"] == 'Avero']

    # p = ggplot(new_data, aes('SNP_FETAL_PCT', 'CHRY_TVALUE', color='FETAL_SEX')) + geom_point(alpha=0.1) \
    #     + labs(x='SNP FF', y='chrY T-Value', title='chrY T value vs SNP FF') + theme_bw()
    #
    # p.save(out_file, format='png', dpi=500)

    model31 = manifest.loc[manifest['SAMPLE_ID'].isin(data['SAMPLE_ID'])] #subset to only the samples in the calls
    model31.loc[((model31['WELL'] == 'E12') & (model31['SAMPLE_TYPE'] == 'Control')), 'FETAL_SEX'] = 'NON PREGNANT MALE'
    #model31.dropna(subset=['FETAL_SEX'], inplace=True)
    model31.to_csv('check.tsv', sep='\t', index=None)
    #model31 = model31.loc[model31["COMPANY"] == 'Avero']

    new_data = model31

    x = np.linspace(-50, 50, 1000)
    norm = stats.norm().pdf(x)

    fig = plt.figure()
    #sns.lineplot(x, norm, color='k', label='Standard Normal')
    male_crit = (new_data["CHRXY_CALL"] == "FETAL EUPLOIDY, MALE") & (new_data["SAMPLE_TYPE"] == "Test")
    sns.kdeplot(new_data.loc[male_crit, 'CHRY_TVALUE'], color='blue', label='MALE')
    print(f'MALE CHRY_TVALUE: {new_data.loc[male_crit, "CHRY_TVALUE"].mean():.2f}')

    female_crit = (new_data["CHRXY_CALL"] == "FETAL EUPLOIDY, FEMALE") & (new_data["SAMPLE_TYPE"] == "Test")
    sns.kdeplot(new_data.loc[female_crit, 'CHRY_TVALUE'], color='red', label='FEMALE')
    print(f'FEMALE CHRY_TVALUE: {new_data.loc[female_crit, "CHRY_TVALUE"].mean():.2f}')

    sns.kdeplot(new_data.loc[new_data['FETAL_SEX'] == 'NON PREGNANT MALE', 'CHRY_TVALUE'], color='black', label='NON PREGNANT MALE')
    print(f'NON PREGNANT MALE CHRY_TVALUE: {new_data.loc[new_data["FETAL_SEX"] == "NON PREGNANT MALE", "CHRY_TVALUE"].mean():.2f}')

    # sns.kdeplot(model31.loc[model31['FETAL_SEX'] == 'MALE', 'CHRY_TVALUE'], color='blue', label='MALE', linestyle='--')
    # print(f'MALE CHRY_TVALUE: {model31.loc[model31["FETAL_SEX"] == "MALE", "CHRY_TVALUE"].mean():.2f}')
    #
    # sns.kdeplot(model31.loc[model31['FETAL_SEX'] == 'FEMALE', 'CHRY_TVALUE'], color='red', label='FEMALE', linestyle='--')
    # print(f'FEMALE CHRY_TVALUE: {model31.loc[model31["FETAL_SEX"] == "FEMALE", "CHRY_TVALUE"].mean():.2f}')
    #
    # sns.kdeplot(model31.loc[model31['FETAL_SEX'] == 'NON PREGNANT MALE', 'CHRY_TVALUE'], color='black', label='NON PREGNANT MALE', linestyle='--')
    # print(f'NON PREGNANT MALE CHRY_TVALUE: {model31.loc[model31["FETAL_SEX"] == "NON PREGNANT MALE", "CHRY_TVALUE"].mean():.2f}')

    plt.xlabel('T Value')
    plt.title('Model 3.1: chrY Null Distributions')

    #plt.vlines(0, 0, 0.42, colors='k', linestyles='dotted')
    # plt.ylim(0, 0.2)
    #plt.ylim(min(norm), 0.42)

    # if male:
    #     plt.xlim(-50, 5)

    plt.legend(frameon=False)

    fig.savefig(out_file, bbox_inches='tight', dpi=250)




def nipt_null_histogram(call_file_path, title, output_file, plot_x=False, by_sex=False):

    data = pd.read_csv(call_file_path, sep='\t', header=0)

    x = np.linspace(-5, 5, 1000)
    norm = stats.norm().pdf(x)

    if not by_sex:
        fig = plt.figure()
        sns.lineplot(x, norm, color='k', label='Standard Normal')
        sns.kdeplot(data.loc[data['CHR13_CALL'] == 'FETAL EUPLOIDY', 'CHR13_TVALUE'], color='r', label='CHR13 T-Value')
        sns.kdeplot(data.loc[data['CHR18_CALL'] == 'FETAL EUPLOIDY', 'CHR18_TVALUE'], color='g', label='CHR18 T-Value')
        sns.kdeplot(data.loc[data['CHR21_CALL'] == 'FETAL EUPLOIDY', 'CHR21_TVALUE'], color='b', label='CHR21 T-Value')

        if plot_x:
            sns.kdeplot(data.loc[data['CHRXY_CALL'] == 'FETAL EUPLOIDY, FEMALE', 'CHRX_TVALUE'], color='orange',
                        label='CHRX T-Value')
            print(f'CHRX_TVALUE: {data.loc[data["CHRXY_CALL"] == "FETAL EUPLOIDY, FEMALE", "CHRX_TVALUE"].mean()}')

        plt.xlabel('T-Value')
        plt.title(title)
        plt.vlines(0,min(norm),0.42, colors='k', linestyles='dotted')
        plt.ylim(min(norm), 0.42)
        plt.legend(frameon=False)

        print(f'CHR13_TVALUE: {data.loc[data["CHR13_CALL"] == "FETAL EUPLOIDY", "CHR13_TVALUE"].mean()}')
        print(f'CHR18_TVALUE: {data.loc[data["CHR18_CALL"] == "FETAL EUPLOIDY", "CHR18_TVALUE"].mean()}')
        print(f'CHR21_TVALUE: {data.loc[data["CHR21_CALL"] == "FETAL EUPLOIDY", "CHR21_TVALUE"].mean()}')


    else:
        fig, [ax1, ax2] = plt.subplots(2, 1, sharex=True)
        sns.lineplot(x, norm, color='k', label='Standard Normal',ax=ax1)
        sns.kdeplot(data.loc[(data['CHR13_CALL'] == 'FETAL EUPLOIDY') & (data['CHRXY_CALL'] == 'FETAL EUPLOIDY, MALE'), 'CHR13_TVALUE'], color='r', label='chr13', ax=ax1)
        sns.kdeplot(data.loc[(data['CHR18_CALL'] == 'FETAL EUPLOIDY') & (data['CHRXY_CALL'] == 'FETAL EUPLOIDY, MALE'), 'CHR18_TVALUE'], color='g', label='chr18', ax=ax1)
        sns.kdeplot(data.loc[(data['CHR21_CALL'] == 'FETAL EUPLOIDY') & (data['CHRXY_CALL'] == 'FETAL EUPLOIDY, MALE'), 'CHR21_TVALUE'], color='b', label='chr21', ax=ax1)
        ax1.set_title(f'{title} - Male Fetuses')
        ax1.get_legend().remove()
        sns.lineplot(x, norm, color='k', label='Standard Normal',ax=ax2)
        sns.kdeplot(data.loc[(data['CHR13_CALL'] == 'FETAL EUPLOIDY') & (data['CHRXY_CALL'] == 'FETAL EUPLOIDY, FEMALE'), 'CHR13_TVALUE'], color='r', label='chr13', ax=ax2)
        sns.kdeplot(data.loc[(data['CHR18_CALL'] == 'FETAL EUPLOIDY') & (data['CHRXY_CALL'] == 'FETAL EUPLOIDY, FEMALE'), 'CHR18_TVALUE'], color='g', label='chr18', ax=ax2)
        sns.kdeplot(data.loc[(data['CHR21_CALL'] == 'FETAL EUPLOIDY') & (data['CHRXY_CALL'] == 'FETAL EUPLOIDY, FEMALE'), 'CHR21_TVALUE'], color='b', label='chr21', ax=ax2)

        if plot_x:
            sns.kdeplot(data.loc[data['CHRXY_CALL'] == 'FETAL EUPLOIDY, FEMALE', 'CHRX_TVALUE'], color='orange', label='chrX', ax=ax2)

        ax2.set_title(f'{title} - Female Fetuses')
        ax2.set_xlabel('T-Value')
        fig.subplots_adjust(bottom=0.3, wspace=0.43)
        ax2.legend(frameon=False, loc='upper center', bbox_to_anchor=(0.5, -0.5), ncol=5)
        fig.tight_layout(pad=1.0)

    fig.savefig(output_file, bbox_inches='tight', dpi=250)


def nipt_null_dist_gof_model(call_file_path, call_file_path2, output_file, by_sex=False):

    data = pd.read_csv(call_file_path, sep='\t', header=0)
    data['CONTROL_SAMPLE'] = data['SAMPLE_ID'].str.split('_').str[-1].str[0]
    data.to_csv('test.tsv', sep='\t', index=None)
    data2 = pd.read_csv(call_file_path2, sep='\t', header=0)

    if by_sex:
        fig, ax = plt.subplots(2, 2, sharex=True)
        sns.kdeplot(data.loc[(data['KNOWN_PLOIDY'] == 'FETAL EUPLOIDY') &
                             (data['CONTROL_SAMPLE'] != 'C') &
                             (data['FETAL_SEX'] == 'MALE'), 'GOF'], color='b', label='Known Negative Samples', ax=ax[0,0])
        sns.kdeplot(data.loc[(data['KNOWN_PLOIDY'] != 'FETAL EUPLOIDY') & (data['KNOWN_PLOIDY'].notnull()) &
                             (data['CONTROL_SAMPLE'] != 'C') &
                             (data['FETAL_SEX'] == 'MALE'), 'GOF'], color='r', label='Known Positive Samples', ax=ax[0,0])
        ax[0, 0].set_title('Model4o: GOF Male')
        ax[0, 0].grid(True)
        ax[0, 0].get_legend().remove()
        ax[0, 0].set_xlim([0.75, 2.0])

        sns.kdeplot(data.loc[(data['KNOWN_PLOIDY'] == 'FETAL EUPLOIDY') &
                             (data['CONTROL_SAMPLE'] != 'C') &
                             (data['FETAL_SEX'] == 'FEMALE'), 'GOF'], color='b', label='Known Negative Samples', ax=ax[0,1])
        sns.kdeplot(data.loc[(data['KNOWN_PLOIDY'] != 'FETAL EUPLOIDY') & (data['KNOWN_PLOIDY'].notnull()) &
                             (data['CONTROL_SAMPLE'] != 'C') &
                             (data['FETAL_SEX'] == 'FEMALE'), 'GOF'], color='r', label='Known Positive Samples', ax=ax[0,1])
        ax[0, 1].set_title('Model4o: GOF Female')
        ax[0, 1].grid(True)
        ax[0, 1].get_legend().remove()

        sns.kdeplot(data2.loc[(data2['KNOWN_PLOIDY'] == 'FETAL EUPLOIDY') &
                              (data2['CONTROL_SAMPLE'] == 'Test')  &
                              (data2['FETAL_SEX'] == 'MALE'), 'GOF'], color='b', label='Known Negative Samples', ax=ax[1,0])
        sns.kdeplot(data2.loc[(data2['KNOWN_PLOIDY'] != 'FETAL EUPLOIDY') & (data2['KNOWN_PLOIDY'].notnull()) &
                              (data2['CONTROL_SAMPLE'] == 'Test') &
                              (data2['FETAL_SEX'] == 'MALE'), 'GOF'], color='r', label='Known Positive Samples', ax=ax[1,0])
        ax[1, 0].set_title('Model3.1: GOF Male')
        ax[1, 0].grid(True)
        ax[1, 0].get_legend().remove()

        sns.kdeplot(data2.loc[(data2['KNOWN_PLOIDY'] == 'FETAL EUPLOIDY') &
                              (data2['CONTROL_SAMPLE'] == 'Test') &
                              (data2['FETAL_SEX'] == 'MALE'), 'GOF'], color='b', label='Known Negative Samples', ax=ax[1,1])
        sns.kdeplot(data2.loc[(data2['KNOWN_PLOIDY'] != 'FETAL EUPLOIDY') & (data2['KNOWN_PLOIDY'].notnull()) &
                              (data2['CONTROL_SAMPLE'] == 'Test') &
                              (data2['FETAL_SEX'] == 'FEMALE'), 'GOF'], color='r', label='Known Positive Samples', ax=ax[1,1])
        ax[1, 1].set_title('Model3.1: GOF Female')
        ax[1, 1].grid(True)
        fig.subplots_adjust(bottom=0.3, wspace=0.43)
        fig.tight_layout(pad=1.0)

        ax[1, 1].legend(frameon=False, loc='upper center', bbox_to_anchor=(0.5, -0.5), ncol=5)

    else:
        fig, [ax1, ax2] = plt.subplots(2, 1, sharex=True)
        sns.kdeplot(data.loc[(data['KNOWN_PLOIDY'] == 'FETAL EUPLOIDY') & (data['CONTROL_SAMPLE'] != 'C'), 'GOF'], color='b', label='Known Negative Samples', ax=ax1)
        sns.kdeplot(data.loc[(data['KNOWN_PLOIDY'] != 'FETAL EUPLOIDY') & (data['KNOWN_PLOIDY'].notnull()) & (data['CONTROL_SAMPLE'] != 'C'), 'GOF'], color='r', label='Known Positive Samples', ax=ax1)
        ax1.set_title('Model4o: GOF Distribution')
        ax1.get_legend().remove()
        ax1.set_xlim([0.75, 2.0])

        sns.kdeplot(data2.loc[(data2['KNOWN_PLOIDY'] == 'FETAL EUPLOIDY') & (data2['CONTROL_SAMPLE'] == 'Test'), 'GOF'], color='b', label='Known Negative Samples', ax=ax2)
        sns.kdeplot(data2.loc[(data2['KNOWN_PLOIDY'] != 'FETAL EUPLOIDY') & (data2['KNOWN_PLOIDY'].notnull())  & (data2['CONTROL_SAMPLE'] == 'Test'), 'GOF'], color='r', label='Known Positive Samples', ax=ax2)
        ax2.set_title('Model3.1: GOF Distribution')
        plt.xlabel('GOF')
        fig.subplots_adjust(bottom=0.3, wspace=0.43)
        ax2.legend(frameon=False, loc='upper center', bbox_to_anchor=(0.5, -0.5), ncol=5)
        fig.tight_layout(pad=1.0)

    fig.savefig(output_file, bbox_inches='tight', dpi=250)


def nipt_null_dist_gof_control(call_file_path, call_file_path2, output_file):
    data = pd.read_csv(call_file_path, sep='\t', header=0)
    data['CONTROL_SAMPLE'] = data['SAMPLE_ID'].str.split('_').str[-1].str[0]
    data['PROPS_ID'] = data['SAMPLE_ID'].str.split('_').str[-1]

    data2 = pd.read_csv(call_file_path2, sep='\t', header=0)


    non_pregnant_male = ['C00259', 'C00260', 'C00261', 'C00262', 'C00263',]

    fig = plt.figure()

    sns.kdeplot(data.loc[data['PROPS_ID'].isin(non_pregnant_male), 'GOF'],
                color='r', label='Model 4o')
    sns.kdeplot(data2.loc[data2['PROPS_ID'].isin(non_pregnant_male), 'GOF'],

    # sns.kdeplot(data.loc[(data['KNOWN_PLOIDY'] == 'FETAL EUPLOIDY') & (data['PROPS_ID'].isin(non_pregnant_male)), 'GOF'],
    #             color='r', label='Model 4o')
    # sns.kdeplot(data2.loc[(data2['KNOWN_PLOIDY'] == 'FETAL EUPLOIDY') & (data2['CONTROL_SAMPLE'] != 'Test' & (data2['PROPS_ID'].isin(non_pregnant_male))), 'GOF'],
                color='b', label='Model 3.1')

    plt.xlabel('GOF')
    plt.title('Male Controls: GOF Distribution')
    plt.legend(frameon=False)
    plt.grid(True)
    fig.savefig(output_file, bbox_inches='tight', dpi=250)

    # #investigate :/
    #
    # control_manifest =  data.loc[data['CONTROL_SAMPLE'] == 'Control']
    # p = ggplot(control_manifest, aes('CHRX_TVALUE', yax, color='Model')) + geom_point(alpha=0.1) \
    #     + labs(x=xlab, y=ylab, title=title) + theme_bw() + stat_smooth(method='lm', se=False)


def nipt_null_dist_gof(call_file_path, title, output_file):

    data = pd.read_csv(call_file_path, sep='\t', header=0)

    fig = plt.figure()
    sns.kdeplot(data.loc[data['KNOWN_PLOIDY'] == 'FETAL EUPLOIDY', 'GOF'], color='b', label='Known Negative Samples')
    sns.kdeplot(data.loc[(data['KNOWN_PLOIDY'] != 'FETAL EUPLOIDY') & (data['KNOWN_PLOIDY'].notnull()), 'GOF'], color='r', label='Known Positive Samples')

    data.loc[(data['KNOWN_PLOIDY'] != 'FETAL EUPLOIDY') & (data['KNOWN_PLOIDY'].notnull())].to_csv('blrop.tsv', sep='\t', index=None)

    plt.xlabel('GOF')
    plt.title(title)
    plt.legend(frameon=False)

    fig.savefig(output_file, bbox_inches='tight', dpi=250)


def arg_parser():
    parser = argparse.ArgumentParser()
    commands = parser.add_subparsers(title='Plots', dest='plot', metavar='')

    null_dist = commands.add_parser('null_dist', help='plot null distribution of t values')
    null_dist.add_argument('call_file', type=Path, help='Path to aggregated call results')
    null_dist.add_argument('plot_title', type=str, help='Plot title')
    null_dist.add_argument('output_path', type=Path, help='Plot output path')
    null_dist.add_argument('-x', action='store_true', default=False, help='Plots chrx t-value distribution of female fetuses')
    null_dist.add_argument('-by_sex', action='store_true', default=False,
                           help='Plots chrx t-value distribution of female fetuses')

    null_dist_x = commands.add_parser('null_dist_x', help='plot null distribution of t values')
    null_dist_x.add_argument('call_file', type=Path, help='Path to aggregated call results')
    null_dist_x.add_argument('call_file2', type=Path, help='Path to aggregated call results')
    null_dist_x.add_argument('plot_title', type=str, help='Plot title')
    null_dist_x.add_argument('output_path', type=Path, help='Plot output path')
    null_dist_x.add_argument('-male', action='store_true', default=False,
                           help='Plot male euploid fetuses')
    # null_dist_x.add_argument('--mod3', type=Path, default=None,
    #                        help='Possible third model to plot')

    null_dist_gof = commands.add_parser('null_dist_gof', help='plot null distribution of t values')
    null_dist_gof.add_argument('call_file', type=Path, help='Path to aggregated call results')
    null_dist_gof.add_argument('plot_title', type=str, help='Plot title')
    null_dist_gof.add_argument('output_path', type=Path, help='Plot output path')
    null_dist_gof.add_argument('-x', action='store_true', default=False, help='Plots chrx t-value distribution of female fetuses')
    null_dist_gof.add_argument('-by_sex', action='store_true', default=False,
                           help='Plots chrx t-value distribution of female fetuses')

    null_dist_gof_mod = commands.add_parser('null_dist_gof_mod', help='plot null distribution of t values')
    null_dist_gof_mod.add_argument('call_file', type=Path, help='Path to aggregated call results')
    null_dist_gof_mod.add_argument('call_file2', type=Path, help='Path to aggregated call results')
    null_dist_gof_mod.add_argument('output_path', type=Path, help='Plot output path')
    null_dist_gof_mod.add_argument('-by_sex', action='store_true', default=False,)

    null_dist_control = commands.add_parser('null_dist_control', help='plot null distribution of t values')
    null_dist_control.add_argument('call_file', type=Path, help='Path to aggregated call results')
    null_dist_control.add_argument('call_file2', type=Path, help='Path to aggregated call results')
    null_dist_control.add_argument('output_path', type=Path, help='Plot output path')

    chry_scatter = commands.add_parser('chry', help='plot null distribution of t values')
    chry_scatter.add_argument('call_file', type=Path, help='Path to aggregated call results')
    chry_scatter.add_argument('manifest', type=Path, help='Path to aggregated call results')
    chry_scatter.add_argument('output_path', type=Path, help='Plot output path')

    scatter = commands.add_parser('scatter', help='plot scatter plots')
    scatter.add_argument('call_file', type=Path, help='Path to aggregated call results')
    scatter.add_argument('call_file2', type=Path, help='Path to aggregated call results')
    scatter.add_argument('output_path', type=Path, help='Plot output path')

    return parser


def cli():

    parser = arg_parser()

    args = parser.parse_args()

    if args.plot == 'null_dist':
        nipt_null_histogram(os.path.abspath(args.call_file), args.plot_title, os.path.abspath(args.output_path), args.x, args.by_sex)

    elif args.plot == 'null_dist_gof':
        nipt_null_dist_gof(os.path.abspath(args.call_file), args.plot_title, os.path.abspath(args.output_path))

    elif args.plot == 'null_dist_gof_mod':
        nipt_null_dist_gof_model(os.path.abspath(args.call_file), os.path.abspath(args.call_file2), os.path.abspath(args.output_path), args.by_sex)

    elif args.plot == 'scatter':
        make_scatter(os.path.abspath(args.call_file), os.path.abspath(args.call_file2), Path(os.path.abspath(args.output_path)))

    elif args.plot == 'null_dist_control':
        nipt_null_dist_gof_control(os.path.abspath(args.call_file), os.path.abspath(args.call_file2), Path(os.path.abspath(args.output_path)))

    elif args.plot == 'null_dist_x':
        nipt_null_x(os.path.abspath(args.call_file), os.path.abspath(args.call_file2),args.plot_title, os.path.abspath(args.output_path), args.male)

    elif args.plot == 'chry':
        chrY_scatter(os.path.abspath(args.call_file), os.path.abspath(args.manifest), os.path.abspath(args.output_path))


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
