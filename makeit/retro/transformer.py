from __future__ import print_function
from global_config import USE_STEREOCHEMISTRY
import rdkit.Chem as Chem          
from rdkit.Chem import AllChem
import numpy as np
from multiprocessing import Pool, cpu_count
from functools import partial # used for passing args to multiprocessing

class Transformer:
    '''
    The Transformer class defines an object which can be used to perform
    one-step retrosyntheses for a given molecule.
    '''

    def __init__(self, parallel=False, nb_workers=None):
        self.source = None
        self.templates = []
        self.has_synth = False
        self.has_retro = False
        self.parallel = parallel 
        self.nb_workers = nb_workers

        # Make sure that nb_workers is okay if we're using a parallel implementation
        if self.parallel and self.nb_workers == None:
            self.nb_workers = cpu_count() - 1
        if self.nb_workers <= 1:
            self.parallel = False 
            self.nb_workers = 1

        # # MEMORY LEAK DEBUG
        # from pympler.tracker import SummaryTracker
        # self.tracker = SummaryTracker()

    def load(self, collection, mincount=4, get_retro=True, get_synth=True, 
        lowe=False, refs=False, efgs=False):
        '''
        Loads the object from a MongoDB collection containing transform
        template records.
        '''
        # Save collection source
        self.source = collection

        # Save get_retro/get_synth:
        if get_retro: self.has_retro = True
        if get_synth: self.has_synth = True

        if mincount and 'count' in collection.find_one(): 
            filter_dict = {'count': { '$gte': mincount}}
        else: 
            filter_dict = {}

        # Look for all templates in collection
        to_retrieve = ['_id', 'reaction_smarts', 'necessary_reagent', 'count', 'intra_only']
        if refs:
            to_retrieve.append('references')
        if efgs:
            to_retrieve.append('efgs')
        for document in collection.find(filter_dict, to_retrieve):
            # Skip if no reaction SMARTS
            if 'reaction_smarts' not in document: continue
            reaction_smarts = str(document['reaction_smarts'])
            if not reaction_smarts: continue

            # Define dictionary
            template = {
                'name':                 document['name'] if 'name' in document else '',
                'reaction_smarts':      reaction_smarts,
                'incompatible_groups':  document['incompatible_groups'] if 'incompatible_groups' in document else [],
                'reference':            document['reference'] if 'reference' in document else '',
                'references':           document['references'] if 'references' in document else [],
                'rxn_example':          document['rxn_example'] if 'rxn_example' in document else '',
                'explicit_H':           document['explicit_H'] if 'explicit_H' in document else False,
                '_id':                  document['_id'] if '_id' in document else -1,
                'product_smiles':       document['product_smiles'] if 'product_smiles' in document else [], 
                'necessary_reagent':    document['necessary_reagent'] if 'necessary_reagent' in document else '',       
                'efgs':                 document['efgs'] if 'efgs' in document else None,
                'intra_only':           document['intra_only'] if 'intra_only' in document else False,
            }

            # Frequency/popularity score
            if 'count' in document: 
                template['count'] = document['count']
            elif 'popularity' in document:
                template['count'] = document['popularity']
            else:
                template['count'] = 1

            # Define reaction in RDKit and validate
            if get_retro:
                try:
                    # Force reactants and products to be one molecule (not really, but for bookkeeping)
                    reaction_smarts_retro = '(' + reaction_smarts.replace('>>', ')>>(') + ')'
                    rxn = AllChem.ReactionFromSmarts(str(reaction_smarts_retro))
                    #if rxn.Validate() == (0, 0):
                    if rxn.Validate()[1] == 0: 
                        template['rxn'] = rxn
                    else:
                        template['rxn'] = None
                except Exception as e:
                    print('Couldnt load retro: {}: {}'.format(reaction_smarts_retro, e))
                    template['rxn'] = None

            # Define forward version, too
            if get_synth:
                try:
                    if lowe:
                        reaction_smarts_synth = '(' + reaction_smarts.split('>')[2] + ')>>(' + reaction_smarts.split('>')[0] + ')'
                    else:
                        reaction_smarts_synth = '(' + reaction_smarts.replace('>>', ')>>(') + ')'
                    rxn_f = AllChem.ReactionFromSmarts(reaction_smarts_synth)
                    #if rxn_f.Validate() == (0, 0):
                    if rxn_f.Validate()[1] == 0:
                        template['rxn_f'] = rxn_f
                    else:
                        template['rxn_f'] = None
                except Exception as e:
                    print('Couldnt load forward: {}: {}'.format(reaction_smarts_synth, e))
                    template['rxn_f'] = None

            # Need to have either a retro or forward reaction be valid
            if get_retro and get_synth:
                if not template['rxn'] and not template['rxn_f']: continue
            elif get_retro:
                if not template['rxn']: continue
            elif get_synth: 
                if not template['rxn_f']: continue

            # Add to list
            self.templates.append(template)
        self.num_templates = len(self.templates)
        self.reorder()

    def reorder(self):
        '''
        Re-orders the list of templates (self.templates) according to 
        field 'count' in descending order. This means we will apply the
        most popular templates first
        '''
        self.templates[:] = [x for x in sorted(self.templates, key = lambda z: z['count'], reverse = True)]

    def top_templates(self, mincount=0):
        '''Generator to return only top templates. 
        Assumes templates are already sorted'''
        for template in self.templates: 
            if template['count'] < mincount: 
                break
            yield template

    def perform_retro(self, smiles, mincount=0):
        '''
        Performs a one-step retrosynthesis given a SMILES string of a
        target molecule by applying each transformation template
        sequentially.
        '''

        # Define mol to operate on
        mol = Chem.MolFromSmiles(smiles)
        smiles = Chem.MolToSmiles(mol, isomericSmiles = USE_STEREOCHEMISTRY) # to canonicalize

        # Initialize results object
        result = RetroResult(smiles)

        if not self.parallel:
            # Try each in turn
            for template in self.top_templates(mincount=mincount):
                for precursor in apply_one_retrotemplate(mol, smiles, template):
                    result.add_precursor(precursor)
        else:
            # Parallel implementation - define what the starting mol is
            apply_one_retrotemplate_to_mol = partial(apply_one_retrotemplate, mol, smiles)
            pool = Pool(self.nb_workers)
            for precursors in pool.imap_unordered(apply_one_retrotemplate_to_mol, self.top_templates(mincount=mincount), chunksize = 50):
                for precursor in precursors:
                    result.add_precursor(precursor)
            pool.close()
            pool.join()

        return result

    def perform_forward(self, smiles, stop_if=None, progbar=False, singleonly=False, mincount=0):
        '''
        Performs a forward synthesis (i.e., reaction enumeration) given
        a SMILES string by applying each transformation template in 
        reverse sequentially

        stop_if - can be used for testing product matching based on 
        if the isomericSmiles matches with one of the products. It terminates
        early instead of going through all of the templates and returns True.
        '''

        # Define pseudo-molecule (single molecule) to operate on
        mol = Chem.MolFromSmiles(smiles)
        smiles = '.'.join(sorted(Chem.MolToSmiles(mol, isomericSmiles = USE_STEREOCHEMISTRY).split('.')))

        # Initialize results object
        result = ForwardResult(smiles)

        # Draw?
        if progbar:
            from tqdm import tqdm
            generator = tqdm(self.top_templates(mincount=mincount))
        else:
            generator = self.top_templates(mincount=mincount)

        # Try each in turn
        for template in generator:
            # Perform
            try:
                outcomes = template['rxn_f'].RunReactants([mol])
            except Exception as e:
                #print('Forward warning: {}'.format(e))
                continue
                #print('Retro version of reaction: {}'.format(template['reaction_smarts']))
            if not outcomes: continue
            for j, outcome in enumerate(outcomes):
                try:
                    for x in outcome:
                        x.UpdatePropertyCache()
                        Chem.SanitizeMol(x)
                except Exception as e:
                    #print(e)
                    continue
                smiles_list = []
                for x in outcome: 
                    smiles_list.extend(Chem.MolToSmiles(x, isomericSmiles = USE_STEREOCHEMISTRY).split('.'))
                # Reduce to largest (longest) product only?
                if singleonly: smiles_list = [max(smiles_list, key = len)]
                product = ForwardProduct(
                    smiles_list = sorted(smiles_list),
                    template_id = template['_id'],
                    num_examples = template['count'],
                )
                if '.'.join(product.smiles_list) == smiles: continue # no transformation
                
                # Early termination?
                if stop_if:
                    if stop_if in product.smiles_list: 
                        print('Found true product - skipping remaining templates to apply')
                        return True
                # If not early termination, we want to keep all products
                else:
                    result.add_product(product)
        # Were we trying to stop early?
        if stop_if: 
            return False
        # Otherwise, return the full list of products
        return result

    def lookup_id(self, template_id):
        '''
        Find the reaction smarts for this template_id
        '''
        for template in self.templates:
            if template['_id'] == template_id:
                return template

class ForwardResult:
    '''
    A class to store the results of a one-step forward synthesis.
    '''

    def __init__(self, smiles):
        self.smiles = smiles 
        self.products = []

    def add_product(self, product):
        '''
        Adds a product to the product set if it is a new product
        '''
        # Check if it is new or old
        for old_product in self.products:
            if product.smiles_list == old_product.smiles_list:
                # Just add this template_id and score
                old_product.template_ids |= set(product.template_ids)
                old_product.num_examples += product.num_examples
                return
        # New!
        self.products.append(product)

    def return_top(self, n = 50):
        '''
        Returns the top n products as a list of dictionaries, 
        sorted by descending score
        '''
        top = []
        for (i, product) in enumerate(sorted(self.products, key = lambda x: x.num_examples, reverse = True)):
            top.append({
                'rank': i + 1,
                'smiles': '.'.join(product.smiles_list),
                'smiles_split': product.smiles_list,
                'num_examples': product.num_examples,
                'tforms': sorted(list(product.template_ids)),
                })
            if i + 1 == n: 
                break
        return top

class ForwardProduct:
    '''
    A class to store a single forward product for reaction enumeration
    '''
    def __init__(self, smiles_list = [], template_id = -1, num_examples = 0):
        self.smiles_list = smiles_list
        self.template_ids = set([template_id])
        self.num_examples = num_examples

class RetroResult:
    '''
    A class to store the results of a one-step retrosynthesis.
    '''
    def __init__(self, target_smiles):
        self.target_smiles = target_smiles
        self.precursors = []

    def add_precursor(self, precursor):
        '''
        Adds a precursor to the retrosynthesis result if it is a new
        and unique product
        '''
        # Check if the precursor set is new or old
        for old_precursor in self.precursors:
            if precursor.smiles_list == old_precursor.smiles_list:
                # Just need to add the fact that this template_id can make it
                old_precursor.template_ids |= set(precursor.template_ids)
                old_precursor.num_examples += precursor.num_examples
                return
        # New! Need to score and add to list
        precursor.score()
        self.precursors.append(precursor)

    def return_top(self, n = 50):
        '''
        Returns the top n precursors as a list of dictionaries, 
        sorted by descending score
        '''
        top = []
        for (i, precursor) in enumerate(sorted(self.precursors, \
                key = lambda x: (x.retroscore, x.num_examples), reverse = True)):
            top.append({
                'rank': i + 1,
                'smiles': '.'.join(precursor.smiles_list),
                'smiles_split': precursor.smiles_list,
                'score': precursor.retroscore,
                'num_examples': precursor.num_examples,
                'tforms': sorted(list(precursor.template_ids)),
                'necessary_reagent': precursor.necessary_reagent,
                })
            if i + 1 == n: 
                break
        return top

class RetroPrecursor:
    '''
    A class to store a single set of precursor(s) for a retrosynthesis
    does NOT contain the target molecule information
    '''
    def __init__(self, smiles_list = [], template_id = -1, num_examples = 0, necessary_reagent = ''):
        self.retroscore = 0
        self.num_examples = num_examples
        self.smiles_list = smiles_list
        self.template_ids = set([template_id])
        self.necessary_reagent = necessary_reagent

    def score(self):
        '''
        Calculate the score of this step as the worst of all precursors,
        plus some penalty for a large necessary_reagent
        '''
        necessary_reagent_atoms = self.necessary_reagent.count('[') / 2.
        scores = []
        for smiles in self.smiles_list:
            x = Chem.MolFromSmiles(smiles)
            total_atoms = x.GetNumHeavyAtoms()
            ring_bonds = sum([b.IsInRing() - b.GetIsAromatic() for b in x.GetBonds()])
            chiral_centers = len(Chem.FindMolChiralCenters(x))

            scores.append(
                - 2.00 * np.power(total_atoms, 1.5) \
                - 1.00 * np.power(ring_bonds, 1.5) \
                - 0.00 * np.power(chiral_centers, 2.0)
            )

        self.retroscore = np.min(scores) - 4.00 * np.power(necessary_reagent_atoms, 2.0)

def apply_one_retrotemplate(mol, smiles, template):
    '''Takes a mol object and applies a single template, returning
    a list of precursors'''
    precursors = []
    try:
        if template['product_smiles']:
            react_mol = Chem.MolFromSmiles(smiles + '.' + '.'.join(template['product_smiles']))
        else:
            react_mol = mol 
        outcomes = template['rxn'].RunReactants([react_mol])
    except Exception as e:
        print('warning: {}'.format(e))
        print(template['reaction_smarts'])
        return []
    for j, outcome in enumerate(outcomes):
        if template['intra_only'] and len(outcome) > 1:
            continue
        try:
            for x in outcome:
                x.UpdatePropertyCache()
                Chem.SanitizeMol(x)
        except Exception as e:
            #print(e) # fail quietly
            continue
        smiles_list = []
        for x in outcome: 
            smiles_list.extend(Chem.MolToSmiles(x, isomericSmiles = USE_STEREOCHEMISTRY).split('.'))
        precursor = RetroPrecursor(
            smiles_list = sorted(smiles_list),
            template_id = template['_id'],
            num_examples = template['count'],
            necessary_reagent = template['necessary_reagent']
        )
        if '.'.join(precursor.smiles_list) == smiles: continue # no transformation
        precursors.append(precursor)
    return precursors