"""This script is used to propose new LQ thresholds for SOP20, these threshold are:

1. Py/FFsnp ratio for ['FETAL EUPLOIDY, MALE', 'MATERNAL XXX, FETAL EUPLOIDY, MALE', 'FETAL XXY'] calls
--> currently set at 0.5

2. Py/FFsnp ratio for FETAL XYY calls
---> currently set at 1

3. FFx/FFsnp ratio for ['FETAL EUPLOIDY, MALE', 'MATERNAL XXX, FETAL EUPLOIDY, MALE', 'FETAL XYY'] calls
---> currently set at 0.5

4. FFx/FFsnp ratio for ['FETAL XXX', 'MATERNAL XXX, FETAL EUPLOIDY, FEMALE', 'FETAL X0', 'MATERNAL X0, FETAL EUPLOIDY, FEMALE'] calls
---> currently set at 0.8

"""
import numpy as np
import pandas as pd
import plotnine          as p9

from pandas.api.types import CategoricalDtype
from pathlib import Path
from scipy         import stats
from sklearn.linear_model import LinearRegression


call_subset_1 = ['FETAL EUPLOIDY, MALE', 'MATERNAL XXX, FETAL EUPLOIDY, MALE', 'FETAL XXY']
call_subset_3 = ['FETAL EUPLOIDY, MALE', 'MATERNAL XXX, FETAL EUPLOIDY, MALE', 'FETAL XYY']
call_subset_4 = ['FETAL XXX', 'MATERNAL XXX, FETAL EUPLOIDY, FEMALE', 'FETAL X0', 'MATERNAL X0, FETAL EUPLOIDY, FEMALE']

def determine_reasons_y(data, r_slope, r_int=0, l_slope=None, l_int=None, snp_cap=None, y_cap=None):
    """This script is used to apply LQ rules specific to the subset1 and XXY calls"""
    data = data.copy(deep=True)
    data['YLOW_REASON'] = np.where((data['CHRY_TVALUE'] < 8) & (data['CHRY_TVALUE'] >= 4), True, False)
    data['RIGHT_RATIO'] = np.where((data['CHRY_PLOIDY'] - r_slope*data['SNP_100'] <= r_int) & ((data['SNP_100'] < snp_cap) | (data['CHRY_PLOIDY'] < y_cap)), True, False)
    data['LEFT_RATIO'] = np.where((data['CHRY_PLOIDY'] - l_slope*data['SNP_100'] >= l_int) & ((data['CHRY_PLOIDY'] < y_cap)), True, False)

    data['LQ_REASON'] = 'None (HQ)'
    data.loc[(data['YLOW_REASON'] == True) & (data['RIGHT_RATIO'] == True), 'LQ_REASON'] = 'Ty + RATIO RIGHT'
    data.loc[(data['YLOW_REASON'] == True) & (data['LEFT_RATIO'] == True), 'LQ_REASON'] = 'Ty + RATIO LEFT)'
    data.loc[(data['YLOW_REASON'] == True) & (data['RIGHT_RATIO'] == False)  & (data['LEFT_RATIO'] == False), 'LQ_REASON'] = 'Ty'
    data.loc[(data['YLOW_REASON'] == False) & (data['RIGHT_RATIO'] == True)  & (data['LEFT_RATIO'] == False), 'LQ_REASON'] = 'RATIO RIGHT'
    data.loc[(data['YLOW_REASON'] == False) & (data['RIGHT_RATIO'] == False)  & (data['LEFT_RATIO'] == True), 'LQ_REASON'] = 'RATIO LEFT'
    data.loc[(data['YLOW_REASON'] == False) & (data['RIGHT_RATIO'] == True)  & (data['LEFT_RATIO'] == True), 'LQ_REASON'] = 'Both Ratios??'


    return data


def determine_reasons_x(data, r_slope, r_int=0, l_slope=None, l_int=None, snp_cap=None, x_cap=None):
    """This script is used to apply LQ rules specific to the subset3 calls"""
    data = data.copy(deep=True)
    data['LQ_REASON'] = 'None (HQ)'

    data['RIGHT_RATIO'] = np.where((data['CHRX_FF_100'] - r_slope*data['SNP_100'] <= r_int) & ((data['SNP_100'] < snp_cap) | (data['CHRX_FF_100'] < x_cap)), True, False)
    data['LEFT_RATIO']  = np.where((data['CHRX_FF_100'] - l_slope*data['SNP_100'] >= l_int) & ((data['CHRX_FF_100'] < x_cap)), True, False)
    
    data.loc[data['RIGHT_RATIO'] == True, 'LQ_REASON'] = 'RATIO RIGHT'
    data.loc[data['LEFT_RATIO'] == True, 'LQ_REASON'] = 'RATIO LEFT'

    return data


def determine_reasons_x2(data, r_slope, r_int=0):
    """This script is used to apply LQ rules specific to the subset4 calls"""
    data = data.copy(deep=True)
    data['XT_REASON'] = False
    data['LQ_REASON'] = 'None (HQ)'

    data.loc[(data['CHRXY_CALL'].isin(['FETAL XXX', 'MATERNAL XXX, FETAL EUPLOIDY, FEMALE'])) & (data['CHRX_TVALUE'] <  5), 'XT_REASON'] = True
    data.loc[(data['CHRXY_CALL'].isin(['FETAL X0', 'MATERNAL X0, FETAL EUPLOIDY, FEMALE']))   & (data['CHRX_TVALUE'] > -6), 'XT_REASON'] = True

    data['RIGHT_RATIO'] = np.where(data['CHRX_FF_100'] - r_slope*data['SNP_100'] <= r_int, True, False)

    data['LQ_REASON'] = 'NEITHER (HQ)'
    data.loc[(data['XT_REASON'] == True) & (data['RIGHT_RATIO'] == True), 'LQ_REASON']   = 'BOTH'
    data.loc[(data['XT_REASON'] == True) & (data['RIGHT_RATIO'] == False), 'LQ_REASON']  = 'X TVALUE'
    data.loc[(data['XT_REASON'] == False) & (data['RIGHT_RATIO'] == True), 'LQ_REASON'] = 'FF_RATIO'

    return data



def main():

    main_dir = Path('/mnt/bfx_projects/nipt_lifecycle/analysis/modelA_B/fits_A4/')
    out_dir = Path('/mnt/bfx_projects/nipt_lifecycle/analysis/modelA_B/fits_A4/ff6_compare_calls/plots/')

    calls_ff5_sop19 = pd.read_csv(main_dir / 'call_test_A4.tsv', sep='\t')
    calls_ff6_sop20 = pd.read_csv(main_dir / 'ff6_compare_calls' / 'call_test_A4_FF6_SOP20.tsv', sep='\t')
    manifest = pd.read_csv('/mnt/bfx_projects/nipt_lifecycle/data/manifest.tsv', sep='\t')

    calls_ff5_sop19 = calls_ff5_sop19.merge(manifest[['SAMPLE_ID', 'SAMPLE_TYPE', 'ANALYSIS_DATETIME', 'RERUN', 'IS_FAIL']], on='SAMPLE_ID', how='left')
    calls_ff6_sop20 = calls_ff6_sop20.merge(manifest[['SAMPLE_ID', 'SAMPLE_TYPE', 'ANALYSIS_DATETIME', 'RERUN', 'IS_FAIL']], on='SAMPLE_ID', how='left')

    calls_ff5_sop19['SNP_100'] = calls_ff5_sop19['SNP_FETAL_PCT'] / 100
    calls_ff6_sop20['SNP_100'] = calls_ff6_sop20['SNP_FETAL_PCT'] / 100

    calls_ff5_sop19['CHRX_FF_100'] = calls_ff5_sop19['CHRX_FETAL_PCT'] / 100
    calls_ff6_sop20['CHRX_FF_100'] = calls_ff6_sop20['CHRX_FETAL_PCT'] / 100

    calls_ff5_sop19_clinical = calls_ff5_sop19.loc[calls_ff5_sop19['SAMPLE_TYPE'] == 'Test']
    calls_ff6_sop20_clinical = calls_ff6_sop20.loc[calls_ff6_sop20['SAMPLE_TYPE'] == 'Test']

    #1. Py/FFsnp ratio for ['FETAL EUPLOIDY, MALE', 'MATERNAL XXX, FETAL EUPLOIDY, MALE', 'FETAL XXY'] calls
    calls_ff6_sop20_clinical_subset_1 = calls_ff6_sop20_clinical.loc[calls_ff6_sop20_clinical['CHRXY_CALL'].isin(call_subset_1)]
    calls_ff6_sop20_LQ_1 = determine_reasons_y(calls_ff6_sop20_clinical_subset_1, 0.75, r_int=-0.025, l_slope=2, l_int=0, snp_cap=0.05, y_cap=0.05)

    calls_ff6_sop20_LQ_1.to_csv(out_dir / 'test.tsv', sep='\t', index=None)

    # plot
    plot = (
        p9.ggplot(calls_ff6_sop20_LQ_1, p9.aes(x='SNP_100', y='CHRY_PLOIDY', color='LQ_REASON')) 
        + p9.geom_point(size=0.2, alpha=0.5)
        + p9.xlab('SNP Fetal Fraction')
        + p9.ylab(f'chrY Ploidy')
        + p9.geom_hline(yintercept = 0.05, alpha=0.5, linetype='dotted')
        + p9.geom_vline(xintercept = 0.05, alpha=0.5, linetype='dotted')
        + p9.geom_abline(intercept=-0.025, slope=0.75, alpha=0.5, linetype = 'dashed')
        + p9.geom_abline(intercept=0, slope=2, alpha=0.5, linetype = 'dashed')
        + p9.ggtitle('Proposed thresholds: chrY Ploidy vs. SNP fetal fraction')
        + p9.scales.ylim(0, 0.2)
        + p9.theme_bw()
    )

    plot.save(out_dir / 'lq1.png', format='png', dpi=300, verbose=False)

    #2. Py/FFsnp ratio for FETAL XYY calls
    calls_ff6_sop20_clinical_subset_1 = calls_ff6_sop20_clinical.loc[calls_ff6_sop20_clinical['CHRXY_CALL'] == 'FETAL XYY']
    calls_ff6_sop20_LQ_2 = determine_reasons_y(calls_ff6_sop20_clinical_subset_1, 0.75, -0.025, 2, 0, 0.05, 0.05)

    # plot
    plot = (
        p9.ggplot(calls_ff6_sop20_LQ_2, p9.aes(x='SNP_100', y='CHRY_PLOIDY', color='LQ_REASON')) 
        + p9.geom_point(size=0.2, alpha=0.5)
        + p9.xlab('SNP Fetal Fraction')
        + p9.ylab(f'chrY Ploidy')
        + p9.geom_hline(yintercept = 0.05, alpha=0.5, linetype='dotted')
        + p9.geom_vline(xintercept = 0.05, alpha=0.5, linetype='dotted')
        + p9.geom_abline(intercept=-0.025, slope=0.75, alpha=0.5, linetype = 'dashed')
        + p9.geom_abline(intercept=0, slope=1, alpha=0.5, color='red')
        + p9.geom_abline(intercept=0, slope=2, alpha=0.5, linetype = 'dashed')
        + p9.ggtitle('Proposed thresholds: chrY Ploidy vs. SNP fetal fraction, XYY calls')
        + p9.scales.ylim(0, 0.2)
        + p9.theme_bw()
    )

    plot.save(out_dir / 'lq2.png', format='png', dpi=300, verbose=False)

    #3. Py/FFsnp ratio for ['FETAL EUPLOIDY, MALE', 'MATERNAL XXX, FETAL EUPLOIDY, MALE', 'FETAL XXY'] calls
    calls_ff6_sop20_clinical_subset_3 = calls_ff6_sop20_clinical.loc[calls_ff6_sop20_clinical['CHRXY_CALL'].isin(call_subset_3)]
    calls_ff6_sop20_LQ_3 = determine_reasons_x(calls_ff6_sop20_clinical_subset_3, 1, r_int=-0.05, l_slope=1.7, l_int=0.025, snp_cap=0.05, x_cap=0.05)

    # plot
    plot = (
        p9.ggplot(calls_ff6_sop20_LQ_3, p9.aes(x='SNP_100', y='CHRX_FF_100', color='LQ_REASON')) 
        + p9.geom_point(size=0.2, alpha=0.5)
        + p9.xlab('SNP Fetal Fraction')
        + p9.ylab(f'chrX Fetal Fraction')
        + p9.geom_hline(yintercept = 0.05, alpha=0.5, linetype='dotted')
        + p9.geom_vline(xintercept = 0.05, alpha=0.5, linetype='dotted')
        + p9.geom_abline(intercept=-0.05, slope=1, alpha=0.5, linetype = 'dashed')
        + p9.geom_abline(intercept=0.025, slope=1.7, alpha=0.5, linetype = 'dashed')
        + p9.ggtitle('Proposed thresholds: chrX Fetal Fraction vs. SNP fetal fraction - Male Calls')
        + p9.scales.ylim(0, 0.2)
        + p9.theme_bw()
    )

    plot.save(out_dir / 'lq3.png', format='png', dpi=300, verbose=False)


    # 4. FFx/FFsnp ratio for ['FETAL XXX', 'MATERNAL XXX, FETAL EUPLOIDY, FEMALE', 'FETAL X0', 'MATERNAL X0, FETAL EUPLOIDY, FEMALE'] calls
    calls_ff6_sop20_clinical_subset_4 = calls_ff6_sop20_clinical.loc[calls_ff6_sop20_clinical['CHRXY_CALL'].isin(call_subset_4)]
    calls_ff6_sop20_LQ_4 = determine_reasons_x2(calls_ff6_sop20_clinical_subset_4, 0.8)

    calls_ff5_sop19_clinical_subset_4 = calls_ff5_sop19_clinical.loc[calls_ff5_sop19_clinical['CHRXY_CALL'].isin(call_subset_4)]
    calls_ff5_sop19_LQ_4 = determine_reasons_x2(calls_ff5_sop19_clinical_subset_4, 0.8)

    # plot
    plot = (
        p9.ggplot(calls_ff6_sop20_LQ_4, p9.aes(x='SNP_100', y='CHRX_FF_100', color='LQ_REASON')) 
        + p9.geom_point(size=0.2, alpha=0.5)
        + p9.xlab('SNP Fetal Fraction')
        + p9.ylab(f'chrX Fetal Fraction')
        + p9.geom_hline(yintercept = 0.05, alpha=0.5, linetype='dotted')
        + p9.geom_vline(xintercept = 0.05, alpha=0.5, linetype='dotted')
        + p9.geom_abline(intercept=0, slope=0.8, alpha=0.5, linetype = 'dashed')
        + p9.ggtitle('FF6: chrX Fetal Fraction vs. SNP fetal fraction - Female Calls')
        + p9.scales.ylim(0, 0.2)
        + p9.theme_bw()
    )

    plot.save(out_dir / 'lq4_ff6.png', format='png', dpi=300, verbose=False)

    # plot facet by call
    plot = (
        p9.ggplot(calls_ff6_sop20_LQ_4, p9.aes(x='SNP_100', y='CHRX_FF_100', color='LQ_REASON')) 
        + p9.geom_point(size=0.2, alpha=0.5)
        + p9.xlab('SNP Fetal Fraction')
        + p9.ylab(f'chrX Fetal Fraction')
        + p9.geom_hline(yintercept = 0.05, alpha=0.5, linetype='dotted')
        + p9.geom_vline(xintercept = 0.05, alpha=0.5, linetype='dotted')
        + p9.geom_abline(intercept=0, slope=0.8, alpha=0.5, linetype = 'dashed')
        + p9.facet_wrap('CHRXY_CALL')
        + p9.ggtitle('FF6: chrX Fetal Fraction vs. SNP fetal fraction - Female Calls')
        + p9.scales.ylim(0, 0.2)
        + p9.theme_bw()
    )

    plot.save(out_dir / 'lq4_facet_ff6.png', format='png', dpi=300, verbose=False)


    # plot facet by call
    plot = (
        p9.ggplot(calls_ff6_sop20_LQ_4, p9.aes(x='SNP_100', y='CHRX_FF_100', color='LQ_REASON')) 
        + p9.geom_point(size=0.2, alpha=0.5)
        + p9.xlab('SNP Fetal Fraction')
        + p9.ylab(f'chrX Fetal Fraction')
        + p9.geom_hline(yintercept = 0.05, alpha=0.5, linetype='dotted')
        + p9.geom_vline(xintercept = 0.05, alpha=0.5, linetype='dotted')
        + p9.geom_abline(intercept=0, slope=0.8, alpha=0.5, linetype = 'dashed')
        + p9.facet_wrap('CHRXY_CALL')
        + p9.scales.ylim(0, 0.2)
        + p9.scales.xlim(0, 0.35)
        + p9.ggtitle('FF6: chrX Fetal Fraction vs. SNP fetal fraction - Female Calls')
        + p9.theme_bw()
    )

    plot.save(out_dir / 'lq4_facet_zoom_ff6.png', format='png', dpi=300, verbose=False)

    # plot the same as above for FF5
    plot = (
        p9.ggplot(calls_ff5_sop19_LQ_4, p9.aes(x='SNP_100', y='CHRX_FF_100', color='LQ_REASON')) 
        + p9.geom_point(size=0.2, alpha=0.5)
        + p9.xlab('SNP Fetal Fraction')
        + p9.ylab(f'chrX Fetal Fraction')
        + p9.geom_hline(yintercept = 0.05, alpha=0.5, linetype='dotted')
        + p9.geom_vline(xintercept = 0.05, alpha=0.5, linetype='dotted')
        + p9.geom_abline(intercept=0, slope=0.8, alpha=0.5, linetype = 'dashed')
        + p9.ggtitle('FF5: chrX Fetal Fraction vs. SNP fetal fraction - Female Calls')
        + p9.scales.ylim(0, 0.2)
        + p9.theme_bw()
    )

    plot.save(out_dir / 'lq4_ff5.png', format='png', dpi=300, verbose=False)

    # plot facet by call - FF5
    plot = (
        p9.ggplot(calls_ff5_sop19_LQ_4, p9.aes(x='SNP_100', y='CHRX_FF_100', color='LQ_REASON')) 
        + p9.geom_point(size=0.2, alpha=0.5)
        + p9.xlab('SNP Fetal Fraction')
        + p9.ylab(f'chrX Fetal Fraction')
        + p9.geom_hline(yintercept = 0.05, alpha=0.5, linetype='dotted')
        + p9.geom_vline(xintercept = 0.05, alpha=0.5, linetype='dotted')
        + p9.geom_abline(intercept=0, slope=0.8, alpha=0.5, linetype = 'dashed')
        + p9.facet_wrap('CHRXY_CALL')
        + p9.ggtitle('FF5: chrX Fetal Fraction vs. SNP fetal fraction - Female Calls')
        + p9.scales.ylim(0, 0.2)
        + p9.theme_bw()
    )

    plot.save(out_dir / 'lq4_facet_ff5.png', format='png', dpi=300, verbose=False)


    # plot facet by call - FF5
    plot = (
        p9.ggplot(calls_ff5_sop19_LQ_4, p9.aes(x='SNP_100', y='CHRX_FF_100', color='LQ_REASON')) 
        + p9.geom_point(size=0.2, alpha=0.5)
        + p9.xlab('SNP Fetal Fraction')
        + p9.ylab(f'chrX Fetal Fraction')
        + p9.geom_hline(yintercept = 0.05, alpha=0.5, linetype='dotted')
        + p9.geom_vline(xintercept = 0.05, alpha=0.5, linetype='dotted')
        + p9.geom_abline(intercept=0, slope=0.8, alpha=0.5, linetype = 'dashed')
        + p9.facet_wrap('CHRXY_CALL')
        + p9.scales.ylim(0, 0.2)
        + p9.scales.xlim(0, 0.35)
        + p9.ggtitle('FF5: chrX Fetal Fraction vs. SNP fetal fraction - Female Calls')
        + p9.theme_bw()
    )

    plot.save(out_dir / 'lq4_facet_zoom_ff5.png', format='png', dpi=300, verbose=False)












if __name__ == '__main__':
    main()
