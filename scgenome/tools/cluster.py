import logging
import sklearn.cluster
import umap
import pandas as pd
import numpy as np
import anndata as ad
from natsort import natsorted

from anndata import AnnData
from typing import Dict

import scgenome.cncluster
import scgenome.preprocessing.transform


def cluster_cells(adata: AnnData, layer_name='copy', method='kmeans') -> AnnData:
    """ Cluster cells by copy number.

    Parameters
    ----------
    adata : AnnData
        copy number data
    layer_name : str, optional
        layer with copy number data to plot, None for X, by default 'state'
    method : str, optional
        clustering method, by default 'kmeans'

    Returns
    -------
    AnnData
        copy number data with additional `cluster_id` column
    """

    if method == 'kmeans':
        cluster_cells_kmeans(adata, layer_name=layer_name)


def cluster_cells_kmeans(adata: AnnData, layer_name='copy', min_k=2, max_k=100) -> AnnData:
    """ Cluster cells by copy number using kmeans.

    Parameters
    ----------
    adata : AnnData
        copy number data
    layer_name : str, optional
        layer with copy number data to plot, None for X, by default 'state'
    min_k : int, optional
        minimum number of clusters, by default 2
    max_k : int, optional
        maximum number of clusters, by default 100

    Returns
    -------
    AnnData
        copy number data with additional `cluster_id` column
    """

    if layer_name is not None:
        X = adata.layers[layer_name]
    else:
        X = adata.X

    X = scgenome.preprocessing.transform.fill_missing(X)

    ks = range(min_k, max_k + 1)

    logging.info(f'trying with max k={max_k}')

    kmeans = []
    bics = []
    for k in ks:
        logging.info(f'trying with k={k}')
        model = sklearn.cluster.KMeans(n_clusters=k, init="k-means++").fit(X)
        bic = scgenome.cncluster.compute_bic(model, X)
        kmeans.append(model)
        bics.append(bic)

    opt_k = np.array(bics).argmax()
    logging.info(f'selected k={opt_k}')

    model = kmeans[opt_k]

    adata.obs['cluster_id'] = model.labels_

    # store information on the clustering parameters
    adata.uns['kmeans'] = {}
    adata.uns['kmeans']['params'] = dict(
        opt_k=opt_k,
        min_k=min_k,
        max_k=max_k,
        layer_name=layer_name,
    )

    return adata


def aggregate_clusters(
        adata: AnnData,
        agg_X: Dict,
        agg_layers: Dict=None,
        agg_obs: Dict=None,
        cluster_col: str='cluster_id') -> AnnData:
    """ Aggregate copy number by cluster to create cluster CN matrix

    Parameters
    ----------
    adata : AnnData
        copy number data
    agg_X : Dict
        [description]
    agg_layers : Dict, optional
        [description], by default None
    agg_obs : Dict, optional
        [description], by default None
    cluster_col : str, optional
        column with cluster ids, by default 'cluster_id'

    Returns
    -------
    AnnData
        aggregsated cluster copy number
    """

    X = (
        adata
            .to_df()
            .set_index(adata.obs[cluster_col].astype(str))
            .groupby(level=0)
            .agg(agg_X)
            .sort_index())

    layer_data = None
    if agg_layers is not None:
        layer_data = {}
        for layer_name in agg_layers:
            layer_data[layer_name] = (
                adata
                    .to_df(layer=layer_name)
                    .set_index(adata.obs[cluster_col].astype(str))
                    .groupby(level=0)
                    .agg(agg_layers[layer_name])
                    .sort_index())

    if agg_obs is not None:
        obs_data = {}
        for obs_name in agg_obs:
            obs_data[obs_name] = (
                adata.obs
                    .set_index(adata.obs[cluster_col].astype(str))[obs_name]
                    .groupby(level=0)
                    .agg(agg_obs[obs_name])
                    .sort_index())
        obs_data = pd.DataFrame(obs_data)

    adata = ad.AnnData(
        X,
        obs=obs_data,
        var=adata.var,
        layers=layer_data,
    )

    return adata


def aggregate_clusters_hmmcopy(adata: AnnData) -> AnnData:
    """ Aggregate hmmcopy copy number by cluster to create cluster CN matrix

    Parameters
    ----------
    adata : AnnData
        hmmcopy copy number data

    Returns
    -------
    AnnData
        aggregsated cluster copy number
    """

    agg_X = np.sum

    agg_layers = {
        'copy': np.nanmean,
        'state': np.nanmedian,
    }

    agg_obs = {
        'total_reads': np.nansum,
    }

    return aggregate_clusters(adata, agg_X, agg_layers, agg_obs, cluster_col='cluster_id')
