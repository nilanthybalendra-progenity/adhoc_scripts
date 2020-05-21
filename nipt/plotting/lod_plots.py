import pandas as pd
import plotnine        as p9
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
            'TVALUE':  data[f'{contig}_TVALUE'],
            'SNP_FETAL_PCT'    :  data['SNP_FETAL_PCT'],
            'MODEL'            :  data['MODEL'],
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

    p = p9.ggplot(aneuploid, p9.aes(x='SNP_FETAL_PCT', y='TVALUE', color='factor(MODEL)')) \
        + p9.geom_point(size=0.3, alpha=0.5) \
        + p9.geom_smooth(method=ols_smoother, se=False, fullrange=True) \
        + p9.scale_x_continuous(name='Fetal Fraction % (SNP)', breaks=np.linspace(0, 10, 11), limits=(0, 10)) \
        + p9.scale_y_continuous(name='T-Value', breaks=np.linspace(0, 15, 16), limits=(0, 15)) \
        + p9.facet_wrap('CONTIG', ncol=2) \
        + p9.theme_bw() \
        + p9.theme(axis_text_y=p9.element_text(size=6)) \
        + p9.theme(legend_position='bottom') \
        + p9.geom_hline(yintercept=4, linetype='dashed', color='red') \
        + p9.labs(color='')

    if xlim:
        p+=p9.scales.xlim(xlim[0], xlim[1])

    p.save(outfile, format='png', dpi=500)

main_dir = Path('/mnt/ruo_rw/rnd/SCRUM_Outputs/NIPT_9002/BFX-1130_NB_model4o/')
manifest = pd.read_csv(main_dir /'manifest_sept2019.tsv', sep='\t', header=0)
model_4o = pd.read_csv(main_dir / 'calls_model4o.tsv', sep='\t', header=0)


manifest['MODEL'] = 'Model 3.1'
model_4o['MODEL'] = 'Model 4.0'
# merge in known ploidy from manifest
#model_4o = model_4o.join(manifest['KNOWN_PLOIDY'], sort=False, how='left', on='SAMPLE_ID')

all = pd.concat([manifest, model_4o], join='inner')


#plot
#all.to_csv('/mnt/ruo_rw/rnd/staff/nilanthy.balendra/tools/adhoc_scripts/sept_plots/sept.tsv', sep='\t')
plot_stuff(all, main_dir / 'LOD_update.png')


# only_reruns = joint_calls.loc[joint_calls['INDIVIDUAL_ID'].isin(rerun_individual_id)]
# plot_stuff(only_reruns, 'joint_calls.png')
#
# plot_stuff(joint_calls.loc[joint_calls['SNP_FETAL_PCT'] > 6], 'all_calls_ff_6.png')
# plot_stuff(only_reruns.loc[only_reruns['SNP_FETAL_PCT'] > 6], 'joint_calls_ff_6.png')
