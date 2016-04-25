from django.shortcuts import render, HttpResponse, redirect
from django.core.urlresolvers import reverse
from django.contrib.auth.decorators import login_required
import django.contrib.auth.views
from forms import SmilesInputForm

# Define transformer
from db import db_client
from django.conf import settings
database = db_client[settings.RETRO_TRANSFORMS['database']]
collection = database[settings.RETRO_TRANSFORMS['collection']]
import makeit.retro.transformer as transformer 
Transformer = transformer.Transformer()
Transformer.load(collection)
print('Loaded {} templates'.format(Transformer.num_templates))

def index(request):
	'''
	Homepage
	'''
	return render(request, 'index.html')

def login(request):
	'''
	User login
	'''
	return django.contrib.auth.views.login(request, template_name = 'login.html')

def logout(request):
	'''
	User logout
	'''
	return django.contrib.auth.views.logout(request, template_name = 'logout.html')

@login_required
def retro(request):
	'''
	Retrosynthesis homepage
	'''
	context = {}

	if request.method == 'POST':
		context['form'] = SmilesInputForm(request.POST) 
		if not context['form'].is_valid():
			context['err'] = 'Could not parse!'
		else:
			# Identify target
			smiles = context['form'].cleaned_data['smiles']
			return redirect('retro_target', smiles = smiles)
	else:
		context['form'] = SmilesInputForm()

	# Define suggestions
	context['suggestions'] = [
		{'name': 'Diphenhydramine', 'smiles': 'CN(C)CCOC(c1ccccc1)c2ccccc2'},
		{'name': 'Fluconazole', 'smiles': 'OC(Cn1cncn1)(Cn2cncn2)c3ccc(F)cc3F'},
		{'name': 'Nevirapine', 'smiles': 'Cc1ccnc2N(C3CC3)c4ncccc4C(=O)Nc12'},
		{'name': 'Atropine', 'smiles': 'CN1C2CCC1CC(C2)OC(=O)C(CO)c3ccccc3'},
		{'name': 'Diazepam', 'smiles': 'CN1C(=O)CN=C(c2ccccc2)c3cc(Cl)ccc13'},
		{'name': 'Hydroxychloroquine', 'smiles': 'CCN(CCO)CCCC(C)Nc1ccnc2cc(Cl)ccc12'},
		{'name': 'Ibuprofen', 'smiles': 'CC(C)Cc1ccc(cc1)C(C)C(O)=O'},
		{'name': 'Tramadol', 'smiles': 'CN(C)C[C@H]1CCCC[C@@]1(C2=CC(=CC=C2)OC)O'},
		{'name': 'Lamivudine', 'smiles': 'NC1=NC(=O)N(C=C1)[C@@H]2CS[C@H](CO)O2'},
		{'name': 'Pregabalin', 'smiles': 'CC(C)C[C@H](CN)CC(O)=O'},
		{'name': 'Naproxen', 'smiles': 'COc1ccc2cc(ccc2c1)C(C)C(O)=O'},
		{'name': 'Imatinib', 'smiles': 'CN1CCN(CC1)Cc2ccc(cc2)C(=O)Nc3ccc(C)c(Nc4nccc(n4)c5cccnc5)c3'},
		{'name': 'Quinapril', 'smiles': 'CCOC(=O)[C@H](CCc1ccccc1)N[C@@H](C)C(=O)N2Cc3ccccc3C[C@H]2C(O)=O'},
		{'name': 'Atorvastatin', 'smiles': 'CC(C)c1n(CC[C@@H](O)C[C@@H](O)CC(O)=O)c(c2ccc(F)cc2)c(c3ccccc3)c1C(=O)Nc4ccccc4'},
		{'name': 'Bortezomib', 'smiles': 'CC(C)C[C@@H](NC(=O)[C@@H](Cc1ccccc1)NC(=O)c2cnccn2)B(O)O'},
		{'name': 'Itraconazole', 'smiles': 'CCC(C)N1N=CN(C1=O)c2ccc(cc2)N3CCN(CC3)c4ccc(OC[C@H]5CO[C@@](Cn6cncn6)(O5)c7ccc(Cl)cc7Cl)cc4'}
	]
	return render(request, 'retro.html', context)

@login_required
def retro_target(request, smiles, max_n = 50):
	'''
	Given a target molecule, render page
	'''

	# Render form with target
	context = {}
	context['form'] = SmilesInputForm({'smiles': smiles})

	# Look up target
	smiles_img = reverse('draw_smiles', kwargs={'smiles':smiles})
	context['target'] = {
		'smiles': smiles,
		'img': smiles_img
	}

	# Perform retrosynthesis
	result = Transformer.perform_retro(smiles)
	context['precursors'] = result.return_top(n = 50)

	# Change 'tform' field to be reaction SMARTS, not ObjectID from Mongo
	for (i, precursor) in enumerate(context['precursors']):
		context['precursors'][i]['tforms'] = \
			[Transformer.lookup_id(_id) for _id in precursor['tforms']]

	return render(request, 'retro.html', context)

@login_required
def draw_smiles(request, smiles):
	'''
	Returns a png response for a target smiles
	'''
	from makeit.retro.draw import MolsSmilesToImage
	response = HttpResponse(content_type = 'img/png')
	MolsSmilesToImage(str(smiles)).save(response, 'png')
	return response

@login_required
def draw_smiles_page(request, smiles):
	'''
	Same as draw_smiles but loads as an indepdent page
	'''
	context = {
		'image_url': reverse('draw_smiles', kwargs={'smiles':smiles}),
		'label_title': 'SMILES',
		'label': smiles,
	}
	return render(request, 'image.html', context)

@login_required
def draw_template(request, template):
	'''
	Returns a png response for a reaction SMARTS template
	'''
	from makeit.retro.draw import TransformStringToImage
	response = HttpResponse(content_type = 'img/png')
	TransformStringToImage(str(template)).save(response, 'png')
	return response

@login_required
def draw_template_page(request, template):
	'''
	Same as draw_template but loads as an indepdent page
	'''
	context = {
		'image_url': reverse('draw_template', kwargs={'template':template}),
		'label_title': 'SMARTS template',
		'label': template,
	}
	return render(request, 'image.html', context)

@login_required
def draw_reaction(request, smiles):
	'''
	Returns a png response for a SMILES reaction string
	'''
	from makeit.retro.draw import ReactionStringToImage
	response = HttpResponse(content_type = 'img/png')
	ReactionStringToImage(str(smiles)).save(response, 'png')
	return response

@login_required
def draw_reaction_page(request, smiles):
	'''
	Same as draw_reaction but loads as an indepdent page
	'''
	context = {
		'image_url': reverse('draw_reaction', kwargs={'smiles':smiles}),
		'label_title': 'Reaction SMILES',
		'label': smiles,
	}
	return render(request, 'image.html', context)