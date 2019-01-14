import matplotlib
import seaborn
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

import refgenome
import utils


def hex_to_rgb(h):
    if h is None:
        return np.array((0, 0, 0))
    h = h.lstrip('#')
    return np.array(tuple(np.uint8(int(h[i:i+2], 16)) for i in (0, 2 ,4)))


color_reference = {0:'#3182BD', 1:'#9ECAE1', 2:'#CCCCCC', 3:'#FDCC8A', 4:'#FC8D59', 5:'#E34A33', 6:'#B30000', 7:'#980043', 8:'#DD1C77', 9:'#DF65B0', 10:'#C994C7', 11:'#D4B9DA'}


def get_cn_cmap(cn_data):
    min_cn = int(cn_data.min())
    max_cn = int(cn_data.max())
    assert min_cn - cn_data.min() == 0
    assert max_cn - cn_data.max() == 0
    color_list = []
    for cn in range(min_cn, max_cn+1):
        if cn > max(color_reference.keys()):
            cn = max(color_reference.keys())
        color_list.append(color_reference[cn])
    return ListedColormap(color_list)


def plot_cbar(ax):
    ax.imshow(np.array([np.arange(len(color_reference))]).T[::-1], cmap=cmap, aspect=1)
    ax.set_xticks([])
    ax.set_yticks(np.arange(len(color_reference)))
    ax.set_yticklabels(np.arange(len(color_reference))[::-1])


def plot_clustered_cell_cn_matrix(ax, cn_data, cn_field_name, raw=False, max_cn=13):
    plot_data = cn_data.merge(utils.chrom_idxs)
    plot_data = plot_data.set_index(['chr_index', 'start', 'cell_id', 'cluster_id'])[cn_field_name].unstack(level=[2, 3]).fillna(0)
    plot_data = plot_data.sort_index(axis=1, level=1)
    if max_cn is not None:
        plot_data[plot_data > max_cn] = max_cn

    mat_chrom_idxs = plot_data.index.get_level_values(0).values
    chrom_boundaries = np.array([0] + list(np.where(mat_chrom_idxs[1:] != mat_chrom_idxs[:-1])[0]) + [plot_data.shape[0] - 1])
    chrom_sizes = chrom_boundaries[1:] - chrom_boundaries[:-1]
    chrom_mids = chrom_boundaries[:-1] + chrom_sizes / 2

    mat_cluster_ids = plot_data.columns.get_level_values(1).values
    cluster_boundaries = np.array([0] + list(np.where(mat_cluster_ids[1:] != mat_cluster_ids[:-1])[0]) + [plot_data.shape[1] - 1])
    cluster_sizes = cluster_boundaries[1:] - cluster_boundaries[:-1]
    cluster_mids = cluster_boundaries[:-1] + cluster_sizes / 2

    cmap = None
    if not raw:
        cmap = get_cn_cmap(plot_data.values)

    im = ax.imshow(plot_data.T, aspect='auto', cmap=cmap)

    ax.set(xticks=chrom_mids)
    ax.set(xticklabels=utils.chrom_names)


def plot_cell_cn_profile(ax, cn_data, value_field_name, cn_field_name, max_cn=13):
    plot_data = cn_data.copy()
    plot_data = plot_data[plot_data['chr'].isin(refgenome.info.chromosomes)]

    plot_data.set_index('chr', inplace=True)
    plot_data['chromosome_start'] = refgenome.info.chromosome_start
    plot_data.reset_index(inplace=True)

    plot_data['start'] = plot_data['start'] + plot_data['chromosome_start']
    plot_data['end'] = plot_data['end'] + plot_data['chromosome_start']

    ax.scatter(
        plot_data['start'], plot_data[value_field_name],
        c=plot_data[cn_field_name], s=1,
        cmap=get_cn_cmap(plot_data[cn_field_name].values),
    )

    ax.set_xlim((-0.5, refgenome.info.chromosome_end.max()))
    ax.set_xlabel('chromosome')
    ax.set_xticks([0] + list(refgenome.info.chromosome_end.values))
    ax.set_xticklabels([])
    ax.set_ylim((-0.5, max_cn))
    ax.set_yticks(np.arange(0, max_cn, 2))
    ax.xaxis.tick_bottom()
    ax.yaxis.tick_left()
    ax.xaxis.set_minor_locator(matplotlib.ticker.FixedLocator(refgenome.info.chromosome_mid))
    ax.xaxis.set_minor_formatter(matplotlib.ticker.FixedFormatter(refgenome.info.chromosomes))

    seaborn.despine(ax=ax, offset=10, trim=True)


def plot_cluster_cn_matrix(ax, cn_data, cn_field_name):
    plot_data = cn_data.merge(utils.chrom_idxs)
    plot_data = plot_data.groupby(['chr_index', 'start', 'cluster_id'])[cn_field_name].median().astype(int)
    plot_data = plot_data.unstack(level=2).fillna(0)
    plot_data = plot_data.sort_index(axis=1, level=1)

    mat_chrom_idxs = plot_data.index.get_level_values(0).values
    chrom_boundaries = np.array([0] + list(np.where(mat_chrom_idxs[1:] != mat_chrom_idxs[:-1])[0]) + [plot_data.shape[0] - 1])
    chrom_sizes = chrom_boundaries[1:] - chrom_boundaries[:-1]
    chrom_mids = chrom_boundaries[:-1] + chrom_sizes / 2

    im = ax.imshow(plot_data.T, aspect='auto', cmap=get_cn_cmap(plot_data.values))

    ax.set(xticks=chrom_mids)
    ax.set(xticklabels=utils.chrom_names)
    ax.set(yticks=range(len(plot_data.columns.values)))
    ax.set(yticklabels=plot_data.columns.values)

    for val in chrom_boundaries[:-1]:
        ax.axvline(x=val, linewidth=1, color='black', zorder=100)


