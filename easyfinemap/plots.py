"""Visualisation of results."""

import matplotlib
import matplotlib.colors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import gridspec
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

from .constant import ColName


def locus_plot(indf: pd.DataFrame, figsize: tuple = (9, 4), size: int = 8, save: bool = False, outprefix: str = ''):
    """Plot the locus plot."""
    sumstat = indf[[ColName.CHR, ColName.BP, ColName.P, 'R2']].copy()
    sumstat['BP'] = sumstat['BP'] / 1e6
    sumstat['P'] = -np.log10(sumstat['P'])
    tmap = matplotlib.colors.LinearSegmentedColormap.from_list(
        'terrain_map_white',
        np.array([[255, 5, 29], [255, 165, 45], [0, 254, 55], [123, 206, 247], [0, 3, 125]]) / 255,
        N=5,
    )

    fig = plt.figure(figsize=figsize)
    gs = gridspec.GridSpec(1, 1)
    ax = plt.subplot(gs[0])
    with_ld = sumstat[sumstat['R2'] > 0]
    ax.scatter(
        sumstat[sumstat['R2'] == -1]['BP'],
        sumstat[sumstat['R2'] == -1]['P'],
        c='grey',
        s=size,
    )
    plt.scatter(
        with_ld['BP'],
        with_ld['P'],
        c=with_ld['R2'],
        cmap=tmap.reversed(),
        s=size,
    )

    ax.set_ylabel('-log10P')
    ax.set_xlabel(f'Position on chromosome {sumstat[ColName.CHR].unique()[0]} (Mb)')

    cbaxes = inset_axes(
        ax,
        width="2%",
        height="28%",
        loc=2,
    )
    plt.colorbar(cax=cbaxes, ticks=[0, 0.2, 0.4, 0.6, 0.8, 1.0])

    # ax.set_xlim([sumstat['BP'].min() - 0.5, sumstat['BP'].max() + 0.5])
    # ax.set_ylim([0, sumstat['P'].max() + 1])
    if save:
        fig.savefig(f'{outprefix}.locus.pdf', bbox_inches='tight')
