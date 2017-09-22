from django.shortcuts import render, HttpResponse, redirect
from django.core.urlresolvers import reverse
from django.contrib.auth.decorators import login_required
from django.http import JsonResponse
from django.conf import settings
import django.contrib.auth.views

from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas

from ..utils import ajax_error_wrapper, resolve_smiles
from ..forms import DrawingInputForm

@ajax_error_wrapper
def ajax_smiles_to_image(request):
    '''Takes an Ajax call with a smiles string
    and returns the HTML for embedding an image'''

    smiles = request.GET.get('smiles', None)
    print('SMILES from Ajax: {}'.format(smiles))
    smiles = resolve_smiles(smiles)
    if smiles is None:
        return JsonResponse({'err': True})
    print('Resolved smiles -> {}'.format(smiles))

    url = reverse('draw_smiles', kwargs={'smiles':smiles})
    data = {
        'html': '<img src="' + url + '">',
        'smiles': smiles,
    }

    return JsonResponse(data)

@ajax_error_wrapper
def ajax_rxn_to_image(request):
    '''Takes an Ajax call with a rxn smiles string
    and returns the HTML for embedding an image'''

    reactants = request.GET.get('reactants', '')
    product = request.GET.get('product', '')

    reactants = resolve_smiles(reactants)
    product = resolve_smiles(product)
    smiles = reactants + '>>' + product
    print('RXN SMILES from Ajax: {}'.format(smiles))
    url = reverse('draw_reaction', kwargs={'smiles':smiles})
    data = {
        'html': '<img src="' + url + '">',
        'smiles': smiles,
        'reactants': reactants,
        'product': product,
    }
    return JsonResponse(data)

@login_required
def draw_smiles(request, smiles):
    '''
    Returns a png response for a target smiles
    '''
    from makeit.retro.draw import MolsSmilesToImage
    response = HttpResponse(content_type = 'image/jpeg')
    MolsSmilesToImage(str(smiles)).save(response, 'jpeg', quality=70)
    return response

@login_required
def draw_template(request, template):
    '''
    Returns a png response for a reaction SMARTS template
    '''
    from makeit.retro.draw import TransformStringToImage
    response = HttpResponse(content_type = 'image/jpeg')
    TransformStringToImage(str(template)).save(response, 'jpeg', quality=70)
    return response

@login_required
def draw_reaction(request, smiles):
    '''
    Returns a png response for a SMILES reaction string
    '''
    from makeit.retro.draw import ReactionStringToImage
    response = HttpResponse(content_type = 'image/jpeg')
    ReactionStringToImage(str(smiles)).save(response, 'jpeg', quality=70)
    return response

@login_required
def draw(request):
    '''
    Landing page for al draw_*_page functions
    '''
    context = {}

    if request.method == 'POST':
        context['form'] = DrawingInputForm(request.POST)
        if not context['form'].is_valid():
            context['err'] = 'Could not parse!'
        else:
            # Identify target
            text = context['form'].cleaned_data['text']
            try:
                if 'mol' in request.POST:
                    #text = resolve_smiles(text)
                    context['image_url'] = reverse('draw_smiles', kwargs={'smiles':text})
                    context['label_title'] = 'Molecule SMILES'
                    context['label'] = text
                elif 'rxn' in request.POST:
                    #text = '>>'.join([resolve_smiles(frag) for frag in text.split('>>')])
                    context['image_url'] = reverse('draw_reaction', kwargs={'smiles':text})
                    context['label_title'] = 'Reaction SMILES'
                    context['label'] = text
                elif 'tform' in request.POST:
                    context['image_url'] = reverse('draw_template', kwargs={'template':text})
                    context['label_title'] = 'Template SMARTS'
                    context['label'] = text
                else:
                    context['err'] = 'Did not understand request'

            except Exception as e:
                context['err'] = e

    else:
        context['form'] = DrawingInputForm()

    return render(request, 'image.html', context)

def draw_fig(request, fig):
    '''
    Returns a png response for figure object
    '''
    response = HttpResponse(content_type = 'img/png')
    canvas = FigureCanvas(fig)
    canvas.print_png(response)
    return response