import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import plotnine

from datetime import datetime as dt

from plotnine          import *
from pathlib import Path

def format_df(data, prod=False):
    euploid_manifest = data.loc[(data['KNOWN_PLOIDY'] == 'FETAL EUPLOIDY') &(data['SAMPLE_TYPE'] == 'Test')]
    #euploid_manifest = data.loc[data['KNOWN_PLOIDY'] == 'FETAL EUPLOIDY']

    if not prod:
        euploid_manifest['ANALYSIS_DATETIME'] = pd.to_datetime(euploid_manifest['ANALYSIS_DATETIME'])
        euploid_manifest['ANALYSIS_WEEK'] = euploid_manifest['ANALYSIS_DATETIME'].apply(
            lambda dt: dt.strftime('%Y-%W') if dt else None)

        avg = euploid_manifest.groupby(['COMPANY', 'ANALYSIS_WEEK'])[['CHR13_TVALUE', 'CHR18_TVALUE', 'CHR21_TVALUE']].mean().reset_index()
        avg = avg.rename(columns={'CHR13_TVALUE': 'AVG_CHR13_TVALUE', 'CHR18_TVALUE': 'AVG_CHR18_TVALUE',
                                  'CHR21_TVALUE': 'AVG_CHR21_TVALUE'})

    else:
        print('prooooooooooooooooooooooood')
        euploid_manifest['ANALYSIS_DATETIME'] = pd.to_datetime(euploid_manifest['ANALYSIS_DATETIME'])
        euploid_manifest['ANALYSIS_WEEK'] = euploid_manifest['ANALYSIS_DATETIME'].apply(
            lambda dt: dt.strftime('%Y-%W') if dt else None)

        avg = euploid_manifest.groupby(['COMPANY', 'ANALYSIS_WEEK'])[['CHR13_TVALUE_PROD', 'CHR18_TVALUE_PROD', 'CHR21_TVALUE_PROD']].mean().reset_index()
        avg = avg.rename(columns={'CHR13_TVALUE_PROD': 'AVG_CHR13_TVALUE', 'CHR18_TVALUE_PROD': 'AVG_CHR18_TVALUE',
                                  'CHR21_TVALUE_PROD': 'AVG_CHR21_TVALUE'})

    return euploid_manifest, avg


def plot_time_series(metric, progenity, avero, outputfile, model, low_y, high_y):
    plt.figure(figsize=[14, 8])
    plt.plot(progenity['ANALYSIS_WEEK'], progenity[metric], 'b.-')
    plt.plot(avero['ANALYSIS_WEEK'], avero[metric], 'r.-')

    plt.legend(('Progenity', 'Avero'), shadow=True, loc='upper right')
    dates = np.sort(progenity['ANALYSIS_WEEK'].unique())
    plt.xticks(range(0, 35,1), dates, fontsize=9.5)
    plt.xticks(rotation=90)

    xmin, xmax, ymin, ymax = plt.axis()
    update_fs = 12

    title = {'AVG_CHR13_TVALUE': f'{model}: Chromosome 13 T-value',
             'AVG_CHR18_TVALUE': f'{model}: Chromosome 18 T-value',
             'AVG_CHR21_TVALUE': f'{model}: Chromosome 21 T-value'}

    midpt = (max(progenity[metric]) + min(progenity[metric])) / 2

    plt.title(title[metric], fontsize=20)
    plt.xlabel('ANALYSIS_WEEK', fontsize=14)

    plt.ylabel(metric, fontsize=14)
    plt.yticks(fontsize=10)
    plt.ylim(low_y, high_y)

    '''
    DATES of Polyphemus Software Version Updates
        v2.0.0 2019-10-08 is Week 41, which is x=7
        v2.1.1 2019-10-14 is Week 42, Day 1 of 2019 which is x=8
        v2.1.3 2019-11-12 is Week 46, Day 4 of 2019 which is x=12
        v2.2.0 2019-12-08 is Week 50 of 2019 which is x=16
        v2.3.0 2019-12-18 is Week 51, Day 3 of 2019 which is x=17
        v2.3.1 2020-02-03 is Week 6 of 2020 which is x=25
    '''

    x = [7, 8, 12, 16, 17, 25]
    y = [0.09, 0.09, 0.09, 0.09, 0.09, 0.09]
    s = ['v2.0.0', 'v2.1.1', 'v2.1.3', 'v2.2.0', 'v2.3.0', 'v2.3.1']

    for i in range(len(x)):
        plt.text(x[i]-1.5, low_y + y[i], s[i], rotation=-45, fontsize=update_fs)
        plt.axvline(x[i], ymax=0.1 - 0.02, color='orange', linestyle='--')

    plt.savefig(outputfile)


#for production
main_dir = Path('/mnt/ruo_rw/rnd/SCRUM_Outputs/NIPT_9002/BFX-1130_NB_model4o/')
file_path = main_dir / 'manifest_sept2019.tsv'
model='Production'
is_prod = True
data = pd.read_csv(file_path, sep='\t', header=0)
euploid_manifest, avg = format_df(data, prod=is_prod)
col_values = ['CHR13_TVALUE_PROD','CHR18_TVALUE_PROD', 'CHR21_TVALUE_PROD']
map_dict = {'CHR13_TVALUE_PROD':'CHR13_TVALUE_PROD_STD','CHR18_TVALUE_PROD':'CHR18_TVALUE_PROD_STD', 'CHR21_TVALUE_PROD':'CHR21_TVALUE_PROD_STD', 'INDIVIDUAL_ID': 'COUNTS'}

# # for new model
# main_dir = Path('/mnt/ruo_rw/rnd/SCRUM_Outputs/NIPT_9002/BFX-1130_NB_model4o/')
# file_path = main_dir / 'calls_model4o.tsv'
# prod_path = main_dir / 'manifest_sept2019.tsv'
# model='Model4o'
# is_prod = False
# col_values = ['CHR13_TVALUE','CHR18_TVALUE', 'CHR21_TVALUE']
#
# model4o = pd.read_csv(file_path, sep='\t', header=0)
# prod = pd.read_csv(prod_path, sep='\t', header=0)
# model4o.drop(labels=['KNOWN_PLOIDY'], axis=1, inplace=True)
# model4o_done = model4o.join(prod[['SAMPLE_ID', 'INDIVIDUAL_ID', 'SAMPLE_TYPE', 'ANALYSIS_DATETIME', 'KNOWN_PLOIDY', 'COMPANY']].set_index(keys='SAMPLE_ID'), sort=False, how='left', on='SAMPLE_ID')
# model4o_done.to_csv(main_dir / 'model40_test.tsv', sep='\t')
# map_dict = {'CHR13_TVALUE':'CHR13_TVALUE_PROD_STD','CHR18_TVALUE':'CHR18_TVALUE_PROD_STD', 'CHR21_TVALUE':'CHR21_TVALUE_PROD_STD', 'INDIVIDUAL_ID': 'COUNTS'}
#
# euploid_manifest, avg = format_df(model4o_done, prod=is_prod)


# calculations

relevant_data = pd.merge(euploid_manifest,avg, on=['COMPANY', 'ANALYSIS_WEEK'])[['INDIVIDUAL_ID', 'AVG_CHR13_TVALUE', 'AVG_CHR18_TVALUE', 'AVG_CHR21_TVALUE', 'ANALYSIS_WEEK', 'COMPANY']]

std_deviations = euploid_manifest.groupby(['COMPANY', 'ANALYSIS_WEEK']).std()[col_values].reset_index(drop=False)

counts = euploid_manifest.groupby(['COMPANY','ANALYSIS_WEEK']).count().reset_index(drop=False)[['COMPANY','ANALYSIS_WEEK', 'INDIVIDUAL_ID']]
calcs = pd.merge(std_deviations, counts, on=['COMPANY', 'ANALYSIS_WEEK']).rename(columns=map_dict)


calcs['CHR13_TVALUE_SE'] = calcs['CHR13_TVALUE_PROD_STD'] / np.sqrt(calcs['COUNTS'])
calcs['CHR18_TVALUE_SE'] = calcs['CHR18_TVALUE_PROD_STD'] / np.sqrt(calcs['COUNTS'])
calcs['CHR21_TVALUE_SE'] = calcs['CHR21_TVALUE_PROD_STD'] / np.sqrt(calcs['COUNTS'])

calculations = pd.merge(calcs,avg, on=['COMPANY', 'ANALYSIS_WEEK'])

sept_cutoff = pd.to_datetime('2019-09-01').strftime('%Y-%W')

dates = np.sort(relevant_data['ANALYSIS_WEEK'].unique())[0::5]
calculations['ANALYSIS_WEEK'].max()


progenity = relevant_data[(relevant_data.COMPANY == 'Progenity')][['AVG_CHR13_TVALUE', 'AVG_CHR18_TVALUE', 'AVG_CHR21_TVALUE', 'ANALYSIS_WEEK', 'COMPANY']]
progenity = progenity.drop_duplicates().sort_values(by=['COMPANY', 'ANALYSIS_WEEK']).reset_index(drop=True)

avero = relevant_data[(relevant_data.COMPANY == 'Avero')][['AVG_CHR13_TVALUE', 'AVG_CHR18_TVALUE', 'AVG_CHR21_TVALUE', 'ANALYSIS_WEEK', 'COMPANY']]
avero = avero.drop_duplicates().sort_values(by=['COMPANY', 'ANALYSIS_WEEK']).reset_index(drop=True)

#plooooooooooooooooooooooooooooooooooooooo


plot_time_series('AVG_CHR13_TVALUE', progenity, avero, main_dir/'time_series' / f'{model}_chr13.png', model,-0.95, 0.05)
plot_time_series('AVG_CHR18_TVALUE', progenity, avero, main_dir/'time_series' / f'{model}_chr18.png', model,-0.6, 0.5)
plot_time_series('AVG_CHR21_TVALUE', progenity, avero, main_dir/'time_series' / f'{model}_chr21.png', model,-0.65, 0.05)