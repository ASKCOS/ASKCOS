from __future__ import print_function

import rdkit.Chem as Chem
from rdkit.Chem import AllChem
import numpy as np

import hdbscan
import sklearn.cluster as cluster

def group_results(original, outcomes, **kwargs):
    '''Cluster the similar transformed outcomes together
    
    Args:
        original (str): SMILES string of original target molecule
        outcomes (list of str): List containing SMILES strings of outcomes to be clustered
        feature (str, optional): Only use features disappearing from 'original', appearing in 'outcomes' or a combination of 'all'. (default: {'original'})
        cluster_method (str, optional): Method to use for clustering ['kmeans', 'hbdscan'] (default: {'kmeans'})
        fp_type (str, optional): Type of fingerprinting method to use. (default: {'morgan'})
        fp_length (int, optional): Fixed-length folding of fingerprint. (default: {512})
        fp_radius (int, optional): Radius to use for fingerprint. (default: {1})
        scores (list or float, optional): Listof precursor outcome scores to number clusters i.e. - cluster 1 contains precursor outcome with best score. (default: {None})

    Returns:
        list of int: Cluster indices for outcomes, 0-based
    '''
    fp_type = kwargs.get('fp_type', 'morgan')
    fp_length = kwargs.get('fp_length', 512)
    fp_radius = kwargs.get('fp_radius', 1)
    if fp_type == 'morgan':
        fp_generator = lambda x: AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(x), fp_radius, nBits=fp_length)
    else:
        raise Exception('Fatal error: fingerprint type {} is not supported.'.format(fingerprint))
    
    cluster_method = kwargs.get('cluster_method', 'kmeans')
    feature = kwargs.get('feature', 'original')
    scores = kwargs.get('scores')

    if not outcomes:
        return []

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

    if scores is not None:
        if len(scores) != len(res):
            raise Exception('Fatal error: length of score ({}) and smiles ({}) are different.'.format(len(scores), len(res)))
        best_cluster_score = {}
        for cluster_id, score in zip(res, scores):
            best_cluster_score[cluster_id] = max(
                best_cluster_score.get(cluster_id, -float('inf')),
                score
            )
        print('best_cluster_score: ', best_cluster_score)
        new_order = list(sorted(best_cluster_score.items(), key=lambda x: -x[1]))
        order_mapping = {new_order[n][0]: n for n in range(len(new_order))}
        print('order_mapping: ', order_mapping)
        res = [order_mapping[n] for n in res]
    
    return res
