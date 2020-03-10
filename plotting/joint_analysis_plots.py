import pandas as pd
import plotnine        as p9
import statsmodels.api as sm
import numpy as np
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
            'TVALUE':  data[f'{contig}_TVALUE'],
            'SNP_FETAL_PCT'    :  data['SNP_FETAL_PCT'],
            'RUN_NUM': data['RUN_NUM'],
        },
        index=data.index,
    )


def plot_stuff(joint_calls, outfile, xlim=None):
    t13 = joint_calls.loc[joint_calls['KNOWN_PLOIDY'].str.contains('TRISOMY 13', na=False)]
    t18 = joint_calls.loc[joint_calls['KNOWN_PLOIDY'].str.contains('TRISOMY 18', na=False)]
    t21 = joint_calls.loc[joint_calls['KNOWN_PLOIDY'].str.contains('TRISOMY 21', na=False)]
    tX = joint_calls.loc[joint_calls['KNOWN_PLOIDY'].str.contains('XXX', na=False)]

    aneuploid = pd.concat([extract(t13, 'CHR13'),
                           extract(t18, 'CHR18'),
                           extract(t21, 'CHR21'),
                           extract(tX, 'CHRX'),
                           ])

    p = p9.ggplot(aneuploid, p9.aes(x='SNP_FETAL_PCT', y='TVALUE', color='factor(RUN_NUM)')) \
        + p9.geom_point(size=0.3, alpha=0.5) \
        + p9.geom_smooth(method=ols_smoother, se=False, fullrange=True) \
        + p9.scale_x_continuous(name='Fetal Fraction % (SNP)', breaks=np.linspace(0, 10, 11), limits=(0, 10)) \
        + p9.scale_y_continuous(name='T-Value', breaks=np.linspace(0, 15, 16), limits=(0, 15)) \
        + p9.facet_wrap('CONTIG', ncol=2) \
        + p9.theme_bw() \
        + p9.theme(axis_text_y=p9.element_text(size=6)) \
        + p9.geom_hline(yintercept=4, linetype='dashed', color='red') \
        + p9.labs(color='')

    if xlim:
        p+=p9.scales.xlim(xlim[0], xlim[1])

    p.save(outfile, format='png', dpi=500)


manifest     = pd.read_csv('/mnt/bfx_projects/nipt_lifecycle/data/wip/manifest.tsv', sep='\t', header=0, index_col=0)
joint_calls  = pd.read_csv('/mnt/bfx_projects/nipt_lifecycle/analysis/wip/joint_testing/joint_run_fits.tsv', sep='\t', header=0, index_col=0)

joint_calls = joint_calls.join(manifest[['RUN_NUM', 'ANALYSIS_DATETIME']], sort=False, how='left')

joint_calls['INDIVIDUAL_ID'] = joint_calls.index.str.split('_').str[-1]
joint_calls.drop(labels=['KNOWN_PLOIDY'], axis=1, inplace=True)
joint_calls.loc[joint_calls['RUN_NUM'].isnull() == True, 'RUN_NUM'] = '1+2'

joint_calls.sort_values(by=['INDIVIDUAL_ID', 'RUN_NUM'], inplace=True)
#joint_calls['SNP_FETAL_PCT'].fillna(method='ffill', inplace=True) # make RUN2 the FF for joint

# average the fetal fraction for join runs
joint_calls['SNP_FETAL_PCT'] =joint_calls.apply(lambda row: row['SNP_FETAL_PCT'] if not pd.isna(row['SNP_FETAL_PCT']) else joint_calls.loc[joint_calls['INDIVIDUAL_ID'] == row['INDIVIDUAL_ID'], 'SNP_FETAL_PCT'].mean(), axis=1)

joint_calls.to_csv('temp.tsv', sep='\t')

# merge in known ploidy fromm manifest
manifest.reset_index(inplace=True)
joint_calls.reset_index(inplace=True)

individuals = manifest.drop_duplicates(subset='INDIVIDUAL_ID', keep='first')
joint_calls = joint_calls.merge(individuals[['INDIVIDUAL_ID', 'KNOWN_PLOIDY']], how='left', left_on='INDIVIDUAL_ID', right_on='INDIVIDUAL_ID')
joint_calls.set_index('SAMPLE_ID', inplace=True)

# list of individuals that were rerun
reruns = manifest.loc[(manifest['RERUN'] == True) & (manifest['CONTROL_SAMPLE'] == 'Test')]
rerun_individual_id = reruns['INDIVIDUAL_ID'].drop_duplicates(keep='first')

#plot
plot_stuff(joint_calls, 'all_calls.png')
joint_calls.to_csv('joint_output_data.tsv', sep='\t')

only_reruns = joint_calls.loc[joint_calls['INDIVIDUAL_ID'].isin(rerun_individual_id)]
plot_stuff(only_reruns, 'joint_calls.png')

plot_stuff(joint_calls.loc[joint_calls['SNP_FETAL_PCT'] > 6], 'all_calls_ff_6.png')
plot_stuff(only_reruns.loc[only_reruns['SNP_FETAL_PCT'] > 6], 'joint_calls_ff_6.png')
