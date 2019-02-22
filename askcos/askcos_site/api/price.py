from django.http import JsonResponse
import makeit.utilities.buyable.pricer as pricer

Pricer = pricer.Pricer()
Pricer.load()

def price(request):
    resp = {}
    resp['request'] = dict(**request.GET)
    smiles = request.GET.get('smiles', None)
    p = Pricer.lookup_smiles(smiles, alreadyCanonical=True)
    resp['price'] = p
    return JsonResponse(resp)