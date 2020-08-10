import argparse
import pandas as pd
import plotnine
from plotnine import *
import statsmodels.api as sm
import numpy as np
from pathlib import Path
#
# import warnings
#
from plotnine.stats.smoothers import wls_prediction_std

def ols_smoother(data, xseq, **params):
    x, y = data['x'], data['y']

    model   = sm.OLS(y, x)
    results = model.fit()

    print(results.params)
    data = pd.DataFrame(
        {
            'x': xseq,
            'y': results.predict(xseq)
        },
    )

    if params['se']:
        prstd, iv_l, iv_u = wls_prediction_std(results, xseq[:,None], alpha=1 - params['level'])
        data['se']        = prstd
        data['ymin']      = iv_l
        data['ymax']      = iv_u

    return data


def extract(data, contig):
    return pd.DataFrame(
        {
            'CONTIG':  contig.replace('CHR', 'chr'),
            'TVALUE':  abs(data[f'{contig}_TVALUE']),
            'SNP_FETAL_PCT'    :  data['SNP_FETAL_PCT'],
            'MODEL'            :  data['MODEL'],
        },
        index=data.index,
    )


def plot_stuff(joint_calls, joint_calls2, outpath, xlim=None):
    t13 = joint_calls.loc[joint_calls['KNOWN_PLOIDY'].str.contains('TRISOMY 13', na=False)]
    t18 = joint_calls.loc[joint_calls['KNOWN_PLOIDY'].str.contains('TRISOMY 18', na=False)]
    t21 = joint_calls.loc[joint_calls['KNOWN_PLOIDY'].str.contains('TRISOMY 21', na=False)]
    tX = joint_calls.loc[(joint_calls['KNOWN_PLOIDY'].str.contains('XXX', na=False)) | (joint_calls['KNOWN_PLOIDY'].str.contains('X0', na=False))]
    #tX0 = joint_calls.loc[joint_calls['KNOWN_PLOIDY'].str.contains('X0', na=False)]

    aneuploid = pd.concat([extract(t13, 'CHR13'),
                           extract(t18, 'CHR18'),
                           extract(t21, 'CHR21'),
                           extract(tX, 'CHRX'),
                           ])

    p = ggplot(aneuploid, aes(x='SNP_FETAL_PCT', y='TVALUE', color='factor(MODEL)')) \
        + geom_point(size=0.3, alpha=0.5) \
        + geom_smooth(method=ols_smoother, se=False, fullrange=True) \
        + scale_x_continuous(name='Fetal Fraction % (SNP)', breaks=np.linspace(0, 10, 11), limits=(0, 10)) \
        + scale_y_continuous(name='T-Value', breaks=np.linspace(0, 15, 16), limits=(0, 15)) \
        + facet_wrap('CONTIG', ncol=2) \
        + theme_bw() \
        + theme(axis_text_y=element_text(size=6)) \
        + theme(legend_position='bottom') \
        + geom_hline(yintercept=4, linetype='dashed', color='red') \
        + labs(color='')

    if xlim:
        p+=scales.xlim(xlim[0], xlim[1])

    p.save(outpath / 'lod.png', format='png', dpi=500)


    t13 = joint_calls2.loc[joint_calls2['KNOWN_PLOIDY'].str.contains('TRISOMY 13', na=False)]
    t18 = joint_calls2.loc[joint_calls2['KNOWN_PLOIDY'].str.contains('TRISOMY 18', na=False)]
    t21 = joint_calls2.loc[joint_calls2['KNOWN_PLOIDY'].str.contains('TRISOMY 21', na=False)]
    tX = joint_calls2.loc[(joint_calls2['KNOWN_PLOIDY'].str.contains('XXX', na=False)) | (joint_calls2['KNOWN_PLOIDY'].str.contains('X0', na=False))]



    chroms = ['13', '18', '21', 'X']
    pos_data = [t13, t18, t21, tX]

    for chr, pos in zip(chroms, pos_data):

        if chr != 'X':

            p1 = ggplot(pos, aes(f'CHR{chr}_TVALUE', f'CHR{chr}_TVALUE_3.1')) \
                 + geom_point(alpha=0.4) + geom_abline(aes(slope=1, intercept=0), color='red') \
                 + labs(x='Model 4f T-Value', y='Model 3.1 T-Value', title=f'CHR{chr} TVALUE for Positives: Original vs Model 4f') \
                 + theme_bw() + plotnine.scales.xlim(0, 8) + plotnine.scales.ylim(0, 8) + geom_hline(yintercept=4, linetype="dashed") + geom_vline(xintercept=4, linetype="dashed")

        else:
            p1 = ggplot(pos, aes(f'CHR{chr}_TVALUE', f'CHR{chr}_TVALUE_3.1', color='KNOWN_PLOIDY')) \
                 + geom_point(alpha=0.4) + geom_abline(aes(slope=1, intercept=0), color='red') \
                 + labs(x='Model 4f T-Value', y='Model 3.1 T-Value', title=f'CHR{chr} TVALUE for Positives: Original vs Model 4f') \
                 + theme_bw() + plotnine.scales.xlim(0, 8) + plotnine.scales.ylim(0, 8)

        pos.to_csv(outpath / f't{chr}.tsv', sep='\t', index=None)

        p1.save(outpath / f'chr{chr}_zoom.png', format='png', dpi=500)


def arg_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('model',    type=str,  help='candidate model name')
    parser.add_argument('manifest', type=Path, help='path to manifest file')
    parser.add_argument('calls',    type=Path, help='path to calls file')
    parser.add_argument('out_path', type=Path, help='output file name')

    return parser

def cli():
    parser = arg_parser()
    args = parser.parse_args()

    manifest = pd.read_csv(args.manifest, sep='\t', header=0)
    model_candidate = pd.read_csv(args.calls,sep='\t', header=0)
    manifest['MODEL'] = 'Model 3.1'
    model_candidate['MODEL'] = f'Model {args.model}'

    manifest = manifest.loc[manifest['SAMPLE_ID'].isin(model_candidate['SAMPLE_ID'])] #subset based on candidate samples
    all_vert = pd.concat([manifest, model_candidate], join='inner') # for LOD plot

    manifest.rename(columns={'CHR13_TVALUE': 'CHR13_TVALUE_3.1',
                             'CHR18_TVALUE': 'CHR18_TVALUE_3.1',
                             'CHR21_TVALUE': 'CHR21_TVALUE_3.1',
                             'CHRX_TVALUE': 'CHRX_TVALUE_3.1',}, inplace=True)

    all_hort = model_candidate.join(manifest[['SAMPLE_ID', 'CHR13_TVALUE_3.1', 'CHR18_TVALUE_3.1', 'CHR21_TVALUE_3.1', 'CHRX_TVALUE_3.1']].set_index('SAMPLE_ID'), on='SAMPLE_ID')
    all_hort.to_csv('temp.tsv', sep='\t', index=None)

    plot_stuff(all_vert, all_hort, args.out_path)






if __name__ == '__main__':
    cli()
