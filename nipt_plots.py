from pathlib import Path

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from scipy import stats


def nipt_null_histogram(data, title, output_file, plot_x=False):

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


def main():

    suffix = 'o'

    fit_file_path = Path(f'calls_model4{suffix}.tsv')
    data = pd.read_csv(fit_file_path, sep='\t', header=0)

    nipt_null_histogram(data, f'Model4{suffix}: Null Distribution', f'model4{suffix}.png')
    nipt_null_histogram(data, f'Model4{suffix}: Null Distribution', f'model4{suffix}_x.png', plot_x=True)


if __name__ == '__main__':
    main()
