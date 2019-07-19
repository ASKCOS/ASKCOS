from __future__ import print_function

import rdkit.Chem as Chem
from rdkit.Chem import AllChem
import numpy as np

import hdbscan
import sklearn.cluster as cluster

def group_results(original, outcomes, **kwargs):
    '''Cluster the similar transformed outcomes together
    Input:
    original:       string,             smiles string
    outcomes:       list of string,     smiles strings of outcomes

    Optional Input:
    feature:        string,             'original', 'outcomes', 'all'. Only use features
                                        presenting in original, outcomes or both.
    fingerprint:    function object,    f(smiles_string:str) -> list of integer bits
    cluster_method: string,             cluster method: hdbscan, kmeans

    Return:
                    list of integer,    cluster indices for outcomes
    '''
    fp_generator = kwargs.get(
        'fingerprint',
        lambda x: AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(x), 1, nBits=512)
        )
    cluster_method = kwargs.get('cluster_method', 'kmeans')
    feature = kwargs.get('feature', 'original')

    # calc fingerprint
    original_fp = np.array(fp_generator(original))
    outcomes_fp = np.array([fp_generator(i) for i in outcomes])

    diff_fp = original_fp - outcomes_fp

    if 'original' == feature:
        diff_fp = diff_fp.clip(min=0)
    elif 'outcomes' == feature:
        diff_fp = diff_fp.clip(max=0)
    elif 'all' == feature:
        pass
    else:
        raise Exception('Fatal error: feature={} is not recognized.'.format(feature))

    # calc culster indices
    res = []
    if 'hdbscan' == cluster_method:
        clusterer = hdbscan.HDBSCAN(min_cluster_size=5, gen_min_span_tree=False)
        clusterer.fit(diff_fp)
        res = clusterer.labels_
        # non-clustered inputs have id -1, make them appear as individual clusters
        max_cluster = np.amax(res)
        for i in range(len(res)):
            if res[i] == -1:
                max_cluster += 1
                res[i] = max_cluster
    elif 'kmeans' == cluster_method:
        for cluster_size in range(len(diff_fp)):
            kmeans = cluster.KMeans(n_clusters=cluster_size+1).fit(diff_fp)
            if kmeans.inertia_ < 1:
                break
        res = kmeans.labels_
    else:
        raise Exception('Fatal error: cluster_method={} is not recognized.'.format(cluster_method))

    res = [int(i) for i in res]
    return res, feature, cluster_method
