import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import plotnine as p9

from datetime import datetime as dt
from mizani.breaks import date_breaks
from mizani.formatters import date_format

from plotnine          import *
from pathlib import Path

def format_df(data, prod=False):
    euploid_manifest = data.loc[(data['KNOWN_PLOIDY'] == 'FETAL EUPLOIDY') &(data['SAMPLE_TYPE'] == 'Test')]

    suf=''
    if prod:
        suf='_PROD'

    col_values = [f'CHR13_TVALUE{suf}', f'CHR18_TVALUE{suf}', f'CHR21_TVALUE{suf}']
    euploid_manifest['ANALYSIS_DATETIME'] = pd.to_datetime(euploid_manifest['ANALYSIS_DATETIME'])
    euploid_manifest['ANALYSIS_WEEK'] = euploid_manifest['ANALYSIS_DATETIME'].apply(lambda dt: dt.strftime('%Y-%W') if dt else None)

    avg = euploid_manifest.groupby(['COMPANY', 'ANALYSIS_WEEK'])[[f'CHR13_TVALUE{suf}', f'CHR18_TVALUE{suf}', f'CHR21_TVALUE{suf}']].mean().reset_index()
    avg = avg.rename(columns={f'CHR13_TVALUE{suf}': f'AVG_CHR13_TVALUE', f'CHR18_TVALUE{suf}': f'AVG_CHR18_TVALUE',
                                  f'CHR21_TVALUE{suf}': f'AVG_CHR21_TVALUE'})

    relevant_data = pd.merge(euploid_manifest, avg, on=['COMPANY', 'ANALYSIS_WEEK'])[
        ['INDIVIDUAL_ID', 'AVG_CHR13_TVALUE', 'AVG_CHR18_TVALUE', 'AVG_CHR21_TVALUE', 'ANALYSIS_WEEK', 'COMPANY']]

    std_deviations = euploid_manifest.groupby(['COMPANY', 'ANALYSIS_WEEK']).std()[col_values].reset_index(drop=False)

    counts = euploid_manifest.groupby(['COMPANY', 'ANALYSIS_WEEK']).count().reset_index(drop=False)[
        ['COMPANY', 'ANALYSIS_WEEK', 'INDIVIDUAL_ID']]

    map_dict = {f'CHR13_TVALUE{suf}': 'CHR13_TVALUE_STD', f'CHR18_TVALUE{suf}': 'CHR18_TVALUE_STD',
                f'CHR21_TVALUE{suf}': 'CHR21_TVALUE_STD', 'INDIVIDUAL_ID': 'COUNTS'}

    calcs = pd.merge(std_deviations, counts, on=['COMPANY', 'ANALYSIS_WEEK']).rename(columns=map_dict)

    calcs['CHR13_TVALUE_SE'] = calcs['CHR13_TVALUE_STD'] / np.sqrt(calcs['COUNTS'])
    calcs['CHR18_TVALUE_SE'] = calcs['CHR18_TVALUE_STD'] / np.sqrt(calcs['COUNTS'])
    calcs['CHR21_TVALUE_SE'] = calcs['CHR21_TVALUE_STD'] / np.sqrt(calcs['COUNTS'])

    calculations = pd.merge(calcs, avg, on=['COMPANY', 'ANALYSIS_WEEK'])

    return calculations


def t_ggplot(calculations, chrom,ymin,ymax, title, poly_model_changes, output_file):
    text_pos =((ymax - ymin) / 2) -1
    k = ggplot(calculations, aes(x='ANALYSIS_WEEK', y=f'AVG_CHR{chrom}_TVALUE', color='COMPANY'))\
        + geom_point() + geom_line(aes(group='COMPANY')) \
        + geom_errorbar(aes(ymin=f'AVG_CHR{chrom}_TVALUE-1.96*CHR{chrom}_TVALUE_SE', ymax=f'AVG_CHR{chrom}_TVALUE+1.96*CHR{chrom}_TVALUE_SE'))\
        + theme_bw() \
        + theme(axis_text=element_text(color='black'),
            axis_text_x=element_text(angle=90, hjust=1, size=6),
            panel_grid_major = element_blank(),
            panel_grid_minor = element_blank(), figure_size=(15,2.5), ) \
        + p9.labels.xlab('Analysis Week') + p9.labels.ylab(f'Mean chr{chrom} T-value')\
        + scale_y_continuous(breaks=np.arange(ymin, ymax, 0.2), limits=(ymin, ymax)) \
        + geom_segment(poly_model_changes, aes(x='ANALYSIS_WEEK', xend='ANALYSIS_WEEK', y=ymin, yend=(ymax)), inherit_aes=False, linetype='--', color='purple') \
        + geom_text(poly_model_changes, aes(x='ANALYSIS_WEEK', y=text_pos, label='version', angle=(90)), inherit_aes=False, ha='right', size=8) \
        + p9.geom_hline(yintercept=0, linetype='dotted', color='grey')

    k.save(output_file, format='png', dpi=500)

def t_ggplot_together(calculations, chrom,ymin,ymax, title, poly_model_changes, output_file, line=False):
    text_pos =((ymax - ymin) / 2) -1

    if line:
        k = ggplot(calculations, aes(x='ANALYSIS_WEEK', y=f'AVG_CHR{chrom}_TVALUE', color='MODEL'))\
            + geom_point() + geom_line(aes(group='MODEL')) \
            + geom_errorbar(aes(ymin=f'AVG_CHR{chrom}_TVALUE-1.96*CHR{chrom}_TVALUE_SE', ymax=f'AVG_CHR{chrom}_TVALUE+1.96*CHR{chrom}_TVALUE_SE'))\
            + theme_bw() \
            + theme(axis_text=element_text(color='black'),
                axis_text_x=element_text(angle=45, hjust=1),
                panel_grid_major = element_blank(),
                panel_grid_minor = element_blank(), figure_size=(15,4)) \
            + p9.labels.xlab('Analysis Week') + p9.labels.ylab(f'Mean chr{chrom} T-value')+p9.ggtitle(title) + p9.labs(color='') \
            + scale_y_continuous(breaks=np.arange(ymin, ymax, 0.1), limits=(ymin, ymax)) \
            + geom_segment(poly_model_changes, aes(x='ANALYSIS_WEEK', xend='ANALYSIS_WEEK', y=ymin, yend=(ymax)), inherit_aes=False, linetype='--', color='purple') \
            + geom_text(poly_model_changes, aes(x='ANALYSIS_WEEK', y=text_pos, label='version', angle=(90)), inherit_aes=False, ha='right')\
            + p9.geom_hline(yintercept=0, linetype='dotted', color='grey')

    else:
        k = ggplot(calculations, aes(x='ANALYSIS_WEEK', y=f'AVG_CHR{chrom}_TVALUE', color='MODEL'))\
            + geom_point() + geom_line(aes(group='MODEL')) \
            + geom_errorbar(aes(ymin=f'AVG_CHR{chrom}_TVALUE-1.96*CHR{chrom}_TVALUE_SE', ymax=f'AVG_CHR{chrom}_TVALUE+1.96*CHR{chrom}_TVALUE_SE'))\
            + theme_bw() \
            + theme(axis_text=element_text(color='black'),
                axis_text_x=element_text(angle=45, hjust=1),
                panel_grid_major = element_blank(),
                panel_grid_minor = element_blank(), figure_size=(15,4)) \
            + p9.labels.xlab('Analysis Week') + p9.labels.ylab(f'Mean chr{chrom} T-value')+p9.ggtitle(title) + p9.labs(color='') \
            + scale_y_continuous(breaks=np.arange(ymin, ymax, 0.1), limits=(ymin, ymax)) \
            + p9.geom_hline(yintercept=0, linetype='dotted', color='grey')

    k.save(output_file, format='png', dpi=500)


def ratio_ggplot(calculations, chrom,title, output_file):
    k = ggplot(calculations, aes(x='ANALYSIS_WEEK', y=f'chr{chrom}_ratio'))\
        + geom_point() + geom_line(aes(group=1)) \
        + theme_bw() \
        + theme(axis_text=element_text(color='black'),
            axis_text_x=element_text(angle=45, hjust=1),
            panel_grid_major = element_blank(),
            panel_grid_minor = element_blank(), figure_size=(15,4)) \
        + p9.labels.xlab('Analysis Week') + p9.labels.ylab(f'Ratio')+p9.ggtitle(title)\
        + p9.geom_hline(yintercept=1, linetype='dotted', color='grey')

    k.save(output_file, format='png', dpi=500)


def main():

    sept_cutoff = pd.to_datetime('2019-09-01').strftime('%Y-%W')

    poly_version_changes = pd.DataFrame({
        'date': ['2018-10-22', '2019-01-10', '2019-02-14', '2019-03-21', '2019-10-14',
                 '2019-11-12', '2019-12-08', '2019-12-18', '2020-02-03'],
        'version': ['v1.0.0', 'v1.0.4', 'v1.1.0', 'v1.2.0', 'v2.1.1',
                    'v2.1.2', 'v2.1.3', 'v2.2.0', 'v2.3.0'],
    })
    poly_version_changes['date'] = pd.to_datetime(poly_version_changes['date'])
    poly_version_changes['ANALYSIS_WEEK'] = poly_version_changes['date'].apply(lambda dt: dt.strftime('%Y-%W') if dt else None)
    #poly_version_changes = poly_version_changes[poly_version_changes['ANALYSIS_WEEK'] >= sept_cutoff]

    poly_model_changes = pd.DataFrame({
        'date': ['2019-05-02', '2019-10-08'],
        'version': ['model v2.0.0', 'model v3.1.0'],
    })

    poly_model_changes['date'] = pd.to_datetime(poly_model_changes['date'])
    poly_model_changes['ANALYSIS_WEEK'] = poly_model_changes['date'].apply(lambda dt: dt.strftime('%Y-%W') if dt else None)
    #poly_model_changes = poly_model_changes[poly_model_changes['ANALYSIS_WEEK'] >= sept_cutoff]

    poly_version_changes['date'] = pd.to_datetime(poly_version_changes['date'])
    poly_version_changes['ANALYSIS_WEEK'] = poly_version_changes['date'].apply(lambda dt: dt.strftime('%Y-%W') if dt else None)
    #poly_version_changes = poly_version_changes[poly_version_changes['ANALYSIS_WEEK'] >= sept_cutoff]

    #main_dir = Path('/mnt/ruo_rw/rnd/SCRUM_Outputs/NIPT_9002/BFX-1130_NB_model4o/')
    main_dir = Path('/mnt/ruo_rw/rnd/SCRUM_Outputs/NIPT_9002/BFX-1296_NB_test_model_4/time_series')
    #prod = pd.read_csv(main_dir / 'manifest_sept2019.tsv', sep='\t', header=0)
    prod = pd.read_csv('/mnt/bfx_projects/nipt_lifecycle/data/manifest.tsv', sep='\t', header=0)
    # model40 = pd.read_csv(main_dir / 'calls_model4o.tsv', sep='\t', header=0)
    # model40.drop(labels=['KNOWN_PLOIDY'], axis=1, inplace=True)
    # model40 = model40.join(prod[['SAMPLE_ID', 'INDIVIDUAL_ID', 'SAMPLE_TYPE', 'ANALYSIS_DATETIME', 'KNOWN_PLOIDY', 'COMPANY']].set_index(keys='SAMPLE_ID'), sort=False, how='left', on='SAMPLE_ID')


    # calculations
    prod_calculations = format_df(prod, prod=True)
    #mod4_calculations = format_df(model40, prod=False)

    chr = ['13', '18', '21']
    # ymin = [-1.3, -1, -1]
    # ymax = [0.2, 0.5, 0.2]

    ymin = [-1.2, -1.2, -1.2]
    ymax = [1.2, 1.2, 1.2]

    for i, c in enumerate(chr):

        t_ggplot(prod_calculations, c, ymin[i], ymax[i], f'Production: chr{c} Mean T-Value by Week', poly_model_changes, main_dir / f'prod_chr{c}.png')
        #t_ggplot(mod4_calculations, c, ymin[i], ymax[i], f'Model 4.0: chr{c} Mean T-Value by Week', poly_model_changes, main_dir / 'time_series'/  f'mod4_chr{c}.png')


    # prod_calculations.loc[prod_calculations['COMPANY'] == 'Avero', 'MODEL']     = 'Avero - Production'
    # prod_calculations.loc[prod_calculations['COMPANY'] == 'Progenity', 'MODEL'] = 'Progenity - Production'
    # mod4_calculations.loc[mod4_calculations['COMPANY'] == 'Avero', 'MODEL']     = 'Avero - Model 4.0'
    # mod4_calculations.loc[mod4_calculations['COMPANY'] == 'Progenity', 'MODEL']     = 'Progenity - Model 4.0'

    # prod_calculations['MODEL'] = 'Production'
    # mod4_calculations['MODEL'] = 'Model 4.0'
    #
    # all_calc = pd.concat([prod_calculations, mod4_calculations], axis=0)
    # all_calc.to_csv('all_cal.tsv', sep='\t', index=None)
    #
    # for i, c in enumerate(chr):
    #     t_ggplot_together(all_calc.loc[all_calc['COMPANY'] == 'Progenity'], c, ymin[i], ymax[i], f'Progenity: chr{c} Mean T-Value by Analysis Week', poly_model_changes, main_dir / 'time_series'/ f'chr{c}_progenity.png', line=True)
    #     t_ggplot_together(all_calc.loc[all_calc['COMPANY'] == 'Avero'], c, ymin[i], ymax[i], f'Avero: chr{c} Mean T-Value by Analysis Week', poly_model_changes, main_dir / 'time_series' / f'chr{c}_avero.png')

    # avero_data     = pd.read_csv(main_dir / 'avero_ratio_4_prod.txt', sep='\t', header=0)
    # progenity_data = pd.read_csv(main_dir / 'progenity_ratio_4_prod.txt', sep='\t', header=0)
    #
    # for c in chr:
    #     ratio_ggplot(avero_data, c,     f'chr{c} Model 4.0:Prod Progenity T-value Ratio', main_dir / 'ratios' / f'chr{c}_avero_ratio.png')
    #     ratio_ggplot(progenity_data, c, f'chr{c} Model 4.0:Prod Avero T-value Ratio',     main_dir / 'ratios' / f'chr{c}_progenity_ratio.png')

if __name__ == '__main__':
    main()