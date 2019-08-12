from django.shortcuts import render, HttpResponse, redirect
from django.urls import reverse
from django.contrib.auth.decorators import login_required
from django.http import JsonResponse
from django.conf import settings
import django.contrib.auth.views

from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas

from ..utils import ajax_error_wrapper, resolve_smiles
from ..forms import DrawingInputForm

import re

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

    url = reverse('draw_smiles', kwargs={'smiles': smiles})
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
    url = reverse('draw_reaction', kwargs={'smiles': smiles})
    data = {
        'html': '<img src="' + url + '">',
        'smiles': smiles,
        'reactants': reactants,
        'product': product,
    }
    return JsonResponse(data)


# @login_required
def draw_smiles(request, smiles):
    '''
    Returns a png response for a target smiles
    '''
    from makeit.utilities.io.draw import MolsSmilesToImage, MakeBackgroundTransparent
    isTransparent = request.GET.get('transparent', 'False')
    response = HttpResponse(content_type='image/png')
    if isTransparent.lower() in ['true', 't', 'yes', 'y', '1']:
        MakeBackgroundTransparent(MolsSmilesToImage(str(smiles))).save(response, 'png', quality=70)
    else:
        MolsSmilesToImage(str(smiles)).save(response, 'png', quality=70)
    return response


#@login_required
def draw_template(request, template):
    '''
    Returns a png response for a reaction SMARTS template
    '''
    from makeit.utilities.io.draw import TransformStringToImage
    response = HttpResponse(content_type='image/jpeg')
    TransformStringToImage(str(template)).save(response, 'png', quality=70)
    return response


# @login_required
def draw_reaction(request, smiles):
    '''
    Returns a png response for a SMILES reaction string
    '''
    from makeit.utilities.io.draw import ReactionStringToImage
    response = HttpResponse(content_type='image/jpeg')
    ReactionStringToImage(str(smiles)).save(response, 'png', quality=70)
    return response

def draw_smiles_highlight(request, smiles, reacting_atoms, bonds='False'):
    '''
    Returns a svg xml with atoms highlighted
    '''
    from makeit.utilities.io.draw import MolsSmilesToImageHighlight
    from ast import literal_eval
    reacting_atoms = literal_eval(reacting_atoms)
    #TODO has to be a better way to evaluate string to true or false
    bonds = bonds.lower() in ['true', '1', 't', 'y', 'yes']
    res = MolsSmilesToImageHighlight(smiles, reacting_atoms=reacting_atoms, bonds=bonds, clear_map=True)
    IsTransparent = request.GET.get('transparent', 'False')
    if IsTransparent.lower() in ['true', '1', 't', 'y', 'yes']:
        # <svg:rect style='opacity:1.0;fill:#FFFFFF;stroke:none' width='300' height='300' x='0' y='0'> </svg:rect>
        # replace #FFFFFF with none
        res = re.sub(r"(.*)<svg:rect(.*)style='([^']*)'(.*)>(.*)</svg:rect>(.*)",
               r"\1<svg:rect\2style='opacity:0.0;fill:none;stroke:none'\4>\5</svg:rect>\6", res)
    response = HttpResponse(res, content_type='image/svg+xml')
  
    return response


#@login_required
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
                    context['image_url'] = reverse(
                        'draw_smiles', kwargs={'smiles': text})
                    context['label_title'] = 'Molecule SMILES'
                    context['label'] = text
                elif 'rxn' in request.POST:
                    #text = '>>'.join([resolve_smiles(frag) for frag in text.split('>>')])
                    context['image_url'] = reverse(
                        'draw_reaction', kwargs={'smiles': text})
                    context['label_title'] = 'Reaction SMILES'
                    context['label'] = text
                elif 'tform' in request.POST:
                    context['image_url'] = reverse(
                        'draw_template', kwargs={'template': text})
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
    response = HttpResponse(content_type='img/png')
    canvas = FigureCanvas(fig)
    canvas.print_png(response)
    return response
