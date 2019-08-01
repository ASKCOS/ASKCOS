from django.http import JsonResponse
from makeit.utilities.cluster import group_results
import rdkit.Chem as Chem
from rdkit.Chem import AllChem

def cluster(request):
    '''Cluster the similar transformed outcomes together
    Method:     GET
    API parameters:
    original:       string,             smiles string
    outcomes:       string,             smiles strings of outcomes, separated by semicolon

    Optional Input:
    feature:        string,             'original', 'outcomes', 'all'. Only use features
                                        presenting in original, outcomes or both.
    fingerprint:    string,             'morgan'
    fpradius:       int,                default is 1
    fpnbits:        int,                default is 512
    clustermethod:  string,             cluster method: 'hdbscan', 'kmeans'

    Return:
                    list of integer,    cluster indices for outcomes

    Testcase:
        %3B=';'
        curl -k 'https://localhost/api/cluster/' -d 'original=CCOC&outcomes=CCO%3BCC'
    '''
    iserr = False
    err_msg = ''
    resp = {}

    if len(request.POST) == 0:
        resp['error'] = 'No input data. Need POST request.'
        return JsonResponse(resp)

    resp['request_POST'] = dict(**request.POST)

    # get parameters
    original = request.POST.get('original', None)
    outcomes = request.POST.get('outcomes', None)
    feature  = request.POST.get('feature', 'original')
    fp_name  = request.POST.get('fingerprint', 'morgan')
    fpradius = int(request.POST.get('fpradius', 1))
    fpnbits  = int(request.POST.get('fpnbits', 512))
    cluster_method = request.POST.get('clustermethod', 'kmeans')

    # error checking
    if original is None:
        iserr = True
        err_msg += 'Error: original is empty. '

    if outcomes is None:
        iserr = True
        err_msg += 'Error: outcomes is empty. '
    else:
        outcomes = outcomes.split(';')

    if feature not in ['original', 'outcomes', 'all']:
        iserr = True
        err_msg += 'Error: unrecognized feature name. '

    if 'morgan' == fp_name:
        fp_generator = lambda x: AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(x), fpradius, nBits=fpnbits)
    else:
        iserr = True
        err_msg += 'Error: unrecognized fingerprint name. '

    if cluster_method not in ['hdbscan', 'kmeans']:
        iserr = True
        err_msg += 'Error: unrecognized clustermethod name. '

    if iserr:
        resp['error'] = err_msg
        return JsonResponse(resp)

    idx, feature, cluster_method = group_results(
        original,
        outcomes,
        feature=feature,
        fingerprint=fp_generator,
        cluster_method=cluster_method
        )
    resp['group_id'] = idx
    resp['feature'] = feature
    resp['clustermethod'] = cluster_method
    resp['fingerprint'] = fp_name
    resp['fpradius'] = fpradius
    resp['fpnbits'] = fpnbits
    return JsonResponse(resp)
