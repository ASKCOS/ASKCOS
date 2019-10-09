from django.http import JsonResponse
import makeit.utilities.buyable.pricer as pricer

Pricer = pricer.Pricer()
Pricer.load()

def price(request):
    resp = {}
    resp['request'] = dict(**request.GET)
    smiles = request.GET.get('smiles', None)
    try:
        p = Pricer.lookup_smiles(smiles, alreadyCanonical=True)
    except Exception as e:
        resp['error'] = str(e)
        res.revoke()
        return JsonResponse(resp, status=400)
    resp['price'] = p
    return JsonResponse(resp)