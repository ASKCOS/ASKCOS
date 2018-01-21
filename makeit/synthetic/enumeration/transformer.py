from __future__ import print_function
import makeit.global_config as gc
USE_STEREOCHEMISTRY = False
import rdkit.Chem as Chem
from rdkit.Chem import AllChem
import numpy as np
import os
import cPickle as pickle
from functools import partial  # used for passing args to multiprocessing
from makeit.utilities.io.logging import MyLogger
from makeit.synthetic.enumeration.results import ForwardResult, ForwardProduct
from pymongo import MongoClient
from makeit.interfaces.template_transformer import TemplateTransformer
from makeit.interfaces.forward_enumerator import ForwardEnumerator
from makeit.prioritization.templates.popularity import PopularityTemplatePrioritizer
from makeit.prioritization.templates.relevance import RelevanceTemplatePrioritizer
from makeit.prioritization.default import DefaultPrioritizer
from makeit.utilities.reactants import clean_reactant_mapping
from makeit.utilities.outcomes import summarize_reaction_outcome

forward_transformer_loc = 'forward_transformer'


class ForwardTransformer(TemplateTransformer, ForwardEnumerator):
    '''
    The Transformer class defines an object which can be used to perform
    one-step retrosyntheses for a given molecule.
    '''

    def __init__(self, mincount=0, TEMPLATE_DB=None, loc=False, done=None, celery=False):
        '''
        Initialize a transformer.
        TEMPLATE_DB: indicate the database you want to use (def. none)
        loc: indicate that local file data should be read instead of online data (def. false)
        '''

        self.done = done
        self.mincount = mincount
        self.templates = []
        self.id_to_index = {}
        self.celery = celery
        self.template_prioritizers = {}
        self.TEMPLATE_DB = TEMPLATE_DB

        super(ForwardTransformer, self).__init__()

    def template_count(self):
        return len(self.templates)

    def get_outcomes(self, smiles, mincount, template_prioritization, start_at=-1, end_at=-1,
                     singleonly=True, stop_if=False, template_count=10000):
        '''
        Each candidate in self.result.products is of type ForwardProduct
        '''
        self.get_template_prioritizers(template_prioritization)
        # Get sorted by popularity during loading.
        if template_prioritization == gc.popularity:
            prioritized_templates = self.templates
        else:
            prioritized_templates = self.template_prioritizer.get_priority(
                (self.templates, smiles), template_count)
        self.mincount = mincount
        self.start_at = start_at
        self.singleonly = singleonly
        self.stop_if = stop_if

        if end_at == -1 or end_at >= len(self.templates):
            self.end_at = len(self.templates)
        else:
            self.end_at = end_at
         # Define mol to operate on
        mol = Chem.MolFromSmiles(smiles)
        clean_reactant_mapping(mol)
        reactants_smiles = Chem.MolToSmiles(mol)
        smiles = Chem.MolToSmiles(
            mol, isomericSmiles=USE_STEREOCHEMISTRY)  # to canonicalize
        # Initialize results object
        if self.celery:
            result = []
        else:
            result = ForwardResult(smiles)
        for i in range(self.start_at, self.end_at):

            # only use templates between the specified boundaries.
            template = prioritized_templates[i]
            if template['count'] > mincount:
                products = self.apply_one_template(
                    mol, smiles, template, singleonly=singleonly, stop_if=stop_if)
                if self.celery:
                    for product in products:
                        result.append({'smiles_list': product.smiles_list,
                                       'smiles': product.smiles,
                                       'edits': product.edits,
                                       'template_ids': product.template_ids,
                                       'num_examples': product.num_examples
                                       })
                else:
                    result.add_products(products)

        return (smiles, result)

    def load(self, chiral=False, lowe=False, refs=False, efgs=False, rxn_ex=False, worker_no=0):
        '''
        Loads and parses the template database to a useable one
        Chrial and rxn_ex are not used, but for compatibility with retro_transformer
        '''
        if worker_no == 0:
            MyLogger.print_and_log('Loading synthetic transformer, including all templates with more than {} hits'.format(
                self.mincount), forward_transformer_loc)

        self.load_templates(False, lowe=lowe, refs=refs, efgs=efgs)

        if worker_no == 0:
            MyLogger.print_and_log('Synthetic transformer has been loaded - using {} templates'.format(
                self.num_templates), forward_transformer_loc)

    def apply_one_template(self, mol, smiles, template, singleonly=True, stop_if=False):
        '''
        Takes a mol object and applies a single template. 
        '''

        try:
            if template['product_smiles']:
                react_mol = Chem.MolFromSmiles(
                    smiles + '.' + '.'.join(template['product_smiles']))
            else:
                react_mol = mol

            outcomes = template['rxn_f'].RunReactants([react_mol])

        except Exception as e:
            if gc.DEBUG:
                MyLogger.print_and_log('Failed transformation for {} because of {}'.format(
                    template['reaction_smarts'], e), forward_transformer_loc, level=1)
            return []

        results = []
        if not outcomes:
            pass
        else:
            for outcome in outcomes:
                smiles_list = []
                # all products represented as single mol by transforms
                outcome = outcome[0]

                try:
                    outcome.UpdatePropertyCache()
                    Chem.SanitizeMol(outcome)
                except Exception as e:
                    if gc.DEBUG:
                        MyLogger.print_and_log('Non-sensible molecule constructed by template {}'.format(
                            template['reaction_smarts']), forward_transformer_loc, level=1)
                    continue
                [a.SetProp(str('molAtomMapNumber'), a.GetProp(str('old_molAtomMapNumber')))
                    for a in outcome.GetAtoms()
                    if str('old_molAtomMapNumber') in a.GetPropsAsDict()]

                # Reduce to largest (longest) product only
                candidate_smiles = Chem.MolToSmiles(
                    outcome, isomericSmiles=True)
                smiles_list = candidate_smiles.split('.')
                if singleonly:
                    candidate_smiles = max(
                        candidate_smiles.split('.'), key=len)
                outcome = Chem.MolFromSmiles(candidate_smiles)

                # Find what edits were made
                edits = summarize_reaction_outcome(react_mol, outcome)

                # Remove mapping before matching
                [x.ClearProp(str('molAtomMapNumber')) for x in outcome.GetAtoms()
                    if x.HasProp(str('molAtomMapNumber'))]  # remove atom mapping from outcome

                # Overwrite candidate_smiles without atom mapping numbers
                candidate_smiles = Chem.MolToSmiles(
                    outcome, isomericSmiles=True)

                product = ForwardProduct(
                    smiles_list=sorted(smiles_list),
                    smiles=candidate_smiles,
                    template_id=str(template['_id']),
                    num_examples=template['count'],
                    edits=edits
                )

                if candidate_smiles == smiles:
                    continue  # no transformation
                if stop_if:
                    if stop_if in product.smiles_list:
                        print(
                            'Found true product - skipping remaining templates to apply')
                        return True
                else:
                    results.append(product)
            # Were we trying to stop early?
            if stop_if:
                return False

        return results

if __name__ == '__main__':
    MyLogger.initialize_logFile()
    ft = ForwardTransformer(mincount=10)
    ft.load()

    template_count = ft.template_count()
    smiles = 'NC(=O)[C@H](CCC=O)N1C(=O)c2ccccc2C1=O'
    for batch_size in range(100, 1000, 100):
        print()
        print(batch_size)
        outcomes = []
        i = 0
        for start_at in range(0, template_count, batch_size):
            i += 1
            outcomes.append(ft.get_outcomes(smiles, 100, start_at=start_at,
                                            end_at=start_at+batch_size, template_prioritization=gc.popularity))
        print('Ran {} batches of {} templates'.format(i, batch_size))
        unique_res = ForwardResult(smiles)

        for smiles, result in outcomes:
            unique_res.add_products(result.products)
        print(len(unique_res.products))
