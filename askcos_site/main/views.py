from django.shortcuts import render, HttpResponse, redirect
from django.core.urlresolvers import reverse
from django.contrib.auth.decorators import login_required
import django.contrib.auth.views
from forms import SmilesInputForm, DrawingInputForm, is_valid_smiles
from bson.objectid import ObjectId
import time

### Retro transformer
from db import db_client
from django.conf import settings
database = db_client[settings.RETRO_TRANSFORMS['database']]
collection = database[settings.RETRO_TRANSFORMS['collection']]
import makeit.retro.transformer as transformer 
RetroTransformer = transformer.Transformer(
	parallel = settings.RETRO_TRANSFORMS['parallel'], 
	nb_workers = settings.RETRO_TRANSFORMS['nb_workers'],
)
mincount_retro = settings.RETRO_TRANSFORMS['mincount']
RetroTransformer.load(collection, mincount = mincount_retro, get_retro = True, get_synth = False)
print('Loaded {} retro templates'.format(RetroTransformer.num_templates))
RETRO_FOOTNOTE = 'Using {} retrosynthesis templates (mincount {}) from {}/{}'.format(RetroTransformer.num_templates,
	mincount_retro, settings.RETRO_TRANSFORMS['database'], settings.RETRO_TRANSFORMS['collection'])

### Forward transformer 
database = db_client[settings.SYNTH_TRANSFORMS['database']]
collection = database[settings.SYNTH_TRANSFORMS['collection']]
SynthTransformer = transformer.Transformer()
mincount_synth = settings.SYNTH_TRANSFORMS['mincount']
SynthTransformer.load(collection, mincount = mincount_synth, get_retro = False, get_synth = True)
print('Loaded {} forward templates'.format(SynthTransformer.num_templates))
SYNTH_FOOTNOTE = 'Using {} forward templates (mincount {}) from {}/{}'.format(SynthTransformer.num_templates,
	mincount_synth, settings.SYNTH_TRANSFORMS['database'], settings.SYNTH_TRANSFORMS['collection'])

### Databases
db = db_client[settings.REACTIONS['database']]
REACTION_DB = db[settings.REACTIONS['collection']]
RETRO_LIT_FOOTNOTE = 'Searched {} known reactions from literature'.format(REACTION_DB.count())

db = db_client[settings.INSTANCES['database']]
INSTANCE_DB = db[settings.INSTANCES['collection']]

db = db_client[settings.CHEMICALS['database']]
CHEMICAL_DB = db[settings.CHEMICALS['collection']]

db = db_client[settings.BUYABLES['database']]
BUYABLE_DB = db[settings.BUYABLES['collection']]

### Prices
print('Loading prices...')
import makeit.retro.pricer as pricer
Pricer = pricer.Pricer()
Pricer.load(CHEMICAL_DB, BUYABLE_DB)
print('Loaded known prices')

### Literaturue transformer
import makeit.retro.transformer_onlyKnown as transformer_onlyKnown
TransformerOnlyKnown = transformer_onlyKnown.TransformerOnlyKnown()
TransformerOnlyKnown.load(CHEMICAL_DB, REACTION_DB)

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
			if 'retro_lit' in request.POST: return redirect('retro_lit_target', smiles = smiles)
			if 'retro' in request.POST: 
				if is_valid_smiles(smiles):
					return redirect('retro_target', smiles = smiles)
				else:
					context['err'] = 'Invalid SMILES string: {}'.format(smiles)
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

	context['footnote'] = RETRO_FOOTNOTE
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
	startTime = time.time()
	result = RetroTransformer.perform_retro(smiles)
	context['precursors'] = result.return_top(n = 50)
	elapsedTime = time.time() - startTime

	# Change 'tform' field to be reaction SMARTS, not ObjectID from Mongo
	# Also add up total number of examples
	for (i, precursor) in enumerate(context['precursors']):
		context['precursors'][i]['tforms'] = \
			[dict(RetroTransformer.lookup_id(_id), **{'id':str(_id)}) for _id in precursor['tforms']]
		context['precursors'][i]['mols'] = []
		for smiles in precursor['smiles_split']:
			ppg = Pricer.lookup_smiles(smiles, alreadyCanonical = True)
			context['precursors'][i]['mols'].append({
				'smiles': smiles,
				'ppg': '${}/g'.format(ppg) if ppg else 'cannot buy'
			})

	context['footnote'] = RETRO_FOOTNOTE + ', generated results in %0.3f seconds' % elapsedTime
	return render(request, 'retro.html', context)

@login_required
def retro_lit_target(request, smiles, max_n = 50):
	'''
	Given a target molecule, render page
	'''

	# Render form with target
	context = {}

	# Perform retrosynthesis
	result = TransformerOnlyKnown.perform_retro(smiles)
	if result:
		context['precursors'] = result.return_top(n = 50)
		# Erase 'tform' field - we have specific RX IDs, not templates
		# Also add up total number of examples
		for (i, precursor) in enumerate(context['precursors']):
			context['precursors'][i]['rxid'] = context['precursors'][i]['tforms'][0]
			del context['precursors'][i]['tforms']
			context['precursors'][i]['mols'] = []
			for smiles in precursor['smiles_split']:
				ppg = Pricer.lookup_smiles(smiles, alreadyCanonical = True)
				context['precursors'][i]['mols'].append({
					'smiles': smiles,
					'ppg': '${}/g'.format(ppg) if ppg else 'cannot buy'
				})
		smiles = result.target_smiles
	else:
		context['err'] = 'Could not find a database match for {}'.format(smiles)

	context['form'] = SmilesInputForm({'smiles': smiles})
	# Look up target
	smiles_img = reverse('draw_smiles', kwargs={'smiles':smiles})
	context['target'] = {
		'smiles': smiles,
		'img': smiles_img
	}
	context['lit_only'] = True
	context['footnote'] = RETRO_LIT_FOOTNOTE
	return render(request, 'retro.html', context)

@login_required
def synth(request):
	'''
	Forward synthesis homepage
	'''
	context = {}

	if request.method == 'POST':
		context['form'] = SmilesInputForm(request.POST) 
		if not context['form'].is_valid():
			context['err'] = 'Could not parse!'
		else:
			# Identify target
			smiles = context['form'].cleaned_data['smiles']
			if is_valid_smiles(smiles):
				return redirect('synth_target', smiles = smiles)
			else:
				context['err'] = 'Invalid SMILES string: {}'.format(smiles)
	else:
		context['form'] = SmilesInputForm()

	context['footnote'] = SYNTH_FOOTNOTE
	return render(request, 'synth.html', context)

@login_required
def synth_target(request, smiles, max_n = 50):
	'''
	Given a set of reactants as a single SMILES string, find products
	'''

	# Render form with target
	context = {}
	context['form'] = SmilesInputForm({'smiles': smiles})

	# Look up target
	smiles_img = reverse('draw_smiles', kwargs={'smiles':smiles})
	context['target'] = {
		'smiles': smiles,
		'img': smiles_img,
	}

	# Perform forward synthesis
	result = SynthTransformer.perform_forward(smiles)
	context['products'] = result.return_top(n = 50)
	#print(context['products'])
	# Change 'tform' field to be reaction SMARTS, not ObjectID from Mongo
	# (add "id" field to be string version of ObjectId "_id")
	for (i, product) in enumerate(context['products']):
		context['products'][i]['tforms'] = \
			[dict(SynthTransformer.lookup_id(_id), **{'id':str(_id)}) for _id in product['tforms']]

	context['footnote'] = SYNTH_FOOTNOTE
	return render(request, 'synth.html', context)

@login_required
def template_target(request, id):
	'''
	Examines a template from its id in the database
	where id == str(_id)

	Reaxys templates refer to instances, but we need the 
	reactions for visualization
	'''
	context = {}

	transform = RetroTransformer.lookup_id(ObjectId(id))
	if not transform: transform = SynthTransformer.lookup_id(ObjectId(id))

	if not transform:
		context['err'] = 'Transform not found'
		return render(request, 'template.html', context)
	context['template'] = transform
	reference_ids = transform['references']

	references = []; docs = []
	context['total_references'] = len(reference_ids)
	for i, ref_id in enumerate(reference_ids):	
		doc = None
		try:
			doc = REACTION_DB.find_one({'_id': int(ref_id.split('-')[0])})
		except Exception as e:
			print(e) 
		if doc: 
			docs.append(doc)
			references.append({
				'label': ref_id,
				'reaction_smiles': doc['RXN_SMILES'],
			})
		else:
			print('Could not find doc with example')

	from makeit.retro.conditions import average_template_list
	context['suggested_conditions'] = average_template_list(INSTANCE_DB, CHEMICAL_DB, reference_ids)

	context['references'] = references[:15]
	return render(request, 'template.html', context)

@login_required
def rxid_target(request, rxid):
	'''
	Examines a reaction record from its id in the database
	where rxid == str(rxid)
	'''
	context = {
		'rxid': rxid
	}

	doc = REACTION_DB.find_one({'_id': int(rxid)})
	if not doc:
		context['err'] = 'RXID not found!'
	else:
		context['rxn_smiles'] = doc['RXN_SMILES']
		context['num_instances'] = doc['RX_NVAR']

	reference_ids = ['{}-{}'.format(rxid, i + 1) for i in range(doc['RX_NVAR'])]
	from makeit.retro.conditions import average_template_list
	context['suggested_conditions'] = average_template_list(INSTANCE_DB, CHEMICAL_DB, reference_ids)

	return render(request, 'rxid.html', context)

@login_required
def price_smiles(request, smiles):
	response = HttpResponse(content_type = 'text/plain')
	ppg = Pricer.lookup_smiles(smiles, alreadyCanonical = True)
	if ppg:
		response.write('${}/g'.format(ppg))
	else:
		response.write('cannot buy')
	return response

@login_required
def price_xrn(request, xrn):
	response = HttpResponse(content_type = 'text/plain')
	ppg = Pricer.lookup_xrn(xrn)
	if ppg:
		response.write('${}/g'.format(ppg))
	else:
		response.write('cannot buy')
	return response

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
	
	OUTDATED: USE PRIMARY /draw/ PAGE (8-31-2016)
	'''
	context = {
		'image_url': reverse('draw_smiles', kwargs={'smiles':smiles}),
		'label_title': 'SMILES',
		'label': smiles,
	}
	if request.method == 'POST':
		context['form'] = DrawingInputForm(request.POST) 
	else:
		context['form'] = DrawingInputForm()
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

	OUTDATED: USE PRIMARY /draw/ PAGE (8-31-2016)
	'''
	context = {
		'image_url': reverse('draw_template', kwargs={'template':template}),
		'label_title': 'SMARTS template',
		'label': template,
	}
	if request.method == 'POST':
		context['form'] = DrawingInputForm(request.POST) 
	else:
		context['form'] = DrawingInputForm()
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

	OUTDATED: USE PRIMARY /draw/ PAGE (8-31-2016)
	'''
	context = {
		'image_url': reverse('draw_reaction', kwargs={'smiles':smiles}),
		'label_title': 'Reaction SMILES',
		'label': smiles,
	}
	if request.method == 'POST':
		context['form'] = DrawingInputForm(request.POST) 
	else:
		context['form'] = DrawingInputForm()
	return render(request, 'image.html', context)

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
					context['image_url'] = reverse('draw_smiles', kwargs={'smiles':text})
					context['label_title'] = 'Molecule SMILES'
					context['label'] = text
				elif 'rxn' in request.POST:
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

@login_required
def draw_synthesis_tree(request):
	'''
	Draw a synthesis tree
	'''

	context = {}

	def chem_dict(smiles, children = []):
		return {
			'is_chemical': True,
			'smiles' : smiles,
			'img' : reverse('draw_smiles', kwargs={'smiles':smiles}),
			'ppg' : Pricer.lookup_smiles(smiles),
			'children': children
		}

	def rxn_dict(info, children = []):
		return {
			'is_reaction': True,
			'info': info,
			'children': children
		}

	target = 'CCCOCCC'
	tree = chem_dict(target, children = [
		rxn_dict('rxn1', children = [
			chem_dict('CCCO'),
			chem_dict('CCC[Br]')
		]),
		rxn_dict('rxn2', children = [
			chem_dict('CCCO'),
			chem_dict('CCC[Cl]', children = [
				rxn_dict('Needs source of [Cl]', children = [
					chem_dict('CCCO'),
				])
			])
		]),
	])

	context['target'] = target
	context['tree'] = tree

	return render(request, 'tree.html', context)