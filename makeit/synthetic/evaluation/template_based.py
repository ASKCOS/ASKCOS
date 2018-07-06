import makeit.global_config as gc
import rdkit.Chem as Chem
import json
import time
from pymongo import MongoClient
import numpy as np
import os
import sys
from celery.result import allow_join_result
from multiprocessing import Queue, Process, Manager
import Queue as VanillaQueue
from makeit.interfaces.scorer import Scorer
from makeit.synthetic.enumeration.transformer import ForwardTransformer
from makeit.synthetic.enumeration.results import ForwardResult, ForwardProduct
from makeit.synthetic.evaluation.template_based_aux import build
import makeit.utilities.contexts as context_cleaner
from makeit.utilities.outcomes import summarize_reaction_outcome
from makeit.utilities.descriptors import edits_to_vectors, edits_to_tensor, edit_vector_lengths
from makeit.utilities.parsing import parse_list_to_smiles
from makeit.utilities.io.logging import MyLogger
from makeit.utilities.reactants import clean_reactant_mapping
from askcos_site.askcos_celery.treeevaluator.forward_trans_worker import get_outcomes, template_count
from operator import itemgetter
template_nn_scorer_loc = 'template_based'


class TemplateNeuralNetScorer(Scorer):

    def __init__(self, celery=False, mincount=0, chunk_size=500, 
                template_prioritization=gc.popularity,
                forward_transformer=None, nproc=1):
        self.model = None
        self.running = False
        self.F_atom = edit_vector_lengths()['atoms']
        self.F_bond = edit_vector_lengths()['bonds']
        self.celery = celery
        self.nproc = nproc
        self.mincount = mincount
        self.chunk_size = chunk_size
        self.template_prioritization = template_prioritization
        self.solvent_name_to_smiles = {}
        self.solvent_smiles_to_params = {}
        self.pending_results = []

        if self.celery:
            tc = template_count.apply_async()
            while not tc.ready():
                time.sleep(0.25)
            self.template_count = tc.get()

        else:
            self.expansion_queue = Queue()
            self.results_queue = Queue()
            if forward_transformer:
                self.forward_transformer = forward_transformer
            else:
                MyLogger.print_and_log(
                    'Cannot evaluate a reaction without a forward transformer. Exiting...', template_nn_scorer_loc, level=3)
            self.template_count = self.forward_transformer.template_count()
            self.workers = []
            self.manager = Manager()
            self.done = self.manager.Value('i', 0)
            self.paused = self.manager.Value('i', 0)
            self.idle = self.manager.list()
            
        if self.celery:
            def expand(reactants_smiles, start_at, end_at):
                self.pending_results.append(get_outcomes.apply_async(args=(reactants_smiles, self.mincount, start_at, end_at,
                                                                           self.template_prioritization),
                                                                     kwargs={'template_count':self.template_count}))
        else:
            def expand(reactants_smiles, start_at, end_at):
                self.expansion_queue.put((reactants_smiles, start_at, end_at))
        self.expand = expand

        if self.celery:
            def waiting_for_results():
                time.sleep(0.05)
                return self.pending_results != []
        else:
            def waiting_for_results():
                waiting = [self.expansion_queue.empty()]
                waiting.append(self.results_queue.empty())
                waiting += self.idle
                return (not all(waiting))
        self.waiting_for_results = waiting_for_results

        if self.celery:
            def get_ready_result(is_ready):

                for i in is_ready:
                    (smiles, outcomes) = self.pending_results[i].get(timeout=0.2)
                    self.pending_results[i].forget()
                    result = ForwardResult(smiles)
                    for outcome in outcomes:
                        result.add_product(ForwardProduct(smiles_list=outcome['smiles_list'], smiles=outcome['smiles'],
                                                          template_ids=outcome[
                                                              'template_ids'],
                                                          num_examples=outcome[
                                                              'num_examples'],
                                                          edits=outcome['edits']))
                    yield result, is_ready
        else:
            def get_ready_result(is_ready):
                while not self.results_queue.empty():
                    [result, start_at, end_at] = self.results_queue.get(0.2)
                    yield result, []

        self.get_ready_result = get_ready_result

        if self.celery:
            def prepare():
                self.running = True
        else:
            def prepare():
                if gc.DEBUG:
                    MyLogger.print_and_log('Template based scorer spinning off {} child processes'.format(
                    self.nproc), template_nn_scorer_loc)
                for i in range(self.nproc):
                    p = Process(target=self.work, args=(i,))
                    self.workers.append(p)
                    p.start()
                self.running = True
                
        self.prepare = prepare
        
        if self.celery:
            def stop_expansion():
                if self.pending_results != []:
                    # OPTION 1 - REVOKE TASKS, WHICH GETS SENT TO ALL WORKERS REGARDLESS OF TYPE
                    [res.revoke() for res in pending_results]
                self.running = False
                
        else:
            def stop_expansion():
                if not self.running:
                    return
                self.done.value = 1
                if gc.DEBUG:
                    MyLogger.print_and_log(
                    'Terminating forward template expansion process.', template_nn_scorer_loc)

                for p in self.workers:
                    if p and p.is_alive():
                        p.terminate()
                if gc.DEBUG:
                    MyLogger.print_and_log(
                    'All forward template expansion processes done.', template_nn_scorer_loc)
                self.running = False
        self.stop_expansion = stop_expansion
        
    def load(self, SOLVENT_DB, folder="", worker_no = 0):
        '''Load a neural network scoring model'''
        if worker_no==0:
            MyLogger.print_and_log('Starting to load scorer...', template_nn_scorer_loc)

        # First load neural network
        if not folder:
            MyLogger.print_and_log(
                'Cannot load neural network without the directory in which the parameters are saved. Exiting...', template_nn_scorer_loc, level=3)
        # Get model args
        ARGS_FPATH = os.path.join(folder, 'args.json')
        with open(ARGS_FPATH, 'r') as fid:
            args = json.load(fid)

        N_h2 = int(args['Nh2'])
        N_h1 = int(args['Nh1'])
        N_h3 = int(args['Nh3'])
        N_hf = int(args['Nhf'])
        l2v = float(args['l2'])
        lr = float(args['lr'])
        context_weight = float(args['context_weight'])
        enhancement_weight = float(args['enhancement_weight'])
        optimizer = args['optimizer']
        inner_act = args['inner_act']
        TARGET_YIELD = False

        self.model = build(F_atom=self.F_atom, F_bond=self.F_bond, N_h1=N_h1,
                           N_h2=N_h2, N_h3=N_h3, N_hf=N_hf, l2v=l2v, inner_act=inner_act,
                           context_weight=context_weight, enhancement_weight=enhancement_weight, TARGET_YIELD=TARGET_YIELD,
                           absolute_score=True)

        WEIGHTS_FPATH = os.path.join(folder, 'weights.h5')
        self.model.load_weights(WEIGHTS_FPATH, by_name=True)

        # Now load solvent information

        for doc in SOLVENT_DB.find():
            try:
                if doc['_id'] == u'default':
                    self.solvent_name_to_smiles['default'] = doc['_id']
                else:
                    self.solvent_name_to_smiles[doc['name']] = doc['_id']

                self.solvent_smiles_to_params[doc['_id']] = doc
            except KeyError:
                MyLogger.print_and_log('Solvent doc {} missing a name'.format(
                    doc), template_nn_scorer_loc, level=1)
        if worker_no == 0:
            MyLogger.print_and_log('Scorer has been loaded.', template_nn_scorer_loc)

    def reset(self):
        if self.celery:
            self.pending_results = []
        else:
            self.done = self.manager.Value('i', 0)
            self.paused = self.manager.Value('i', 0)
            self.idle = self.manager.list()
            self.results_queue = Queue()
            self.workers = []

    def work(self, i):
        while True:
            # If done, stop
            if self.done.value:
                MyLogger.print_and_log(
                    'Worker {} saw done signal, terminating'.format(i), template_nn_scorer_loc)
                break
            # If paused, wait and check again
            if self.paused.value:
                #print('Worker {} saw pause signal, sleeping for 1 second'.format(i))
                time.sleep(1)
                continue
            # Grab something off the queue
            try:
                (reactants_smiles, start_at, end_at) = self.expansion_queue.get(
                    timeout=0.5)  # short timeout
                self.idle[i] = False
                (smiles, result) = self.forward_transformer.get_outcomes(reactants_smiles, self.mincount, self.template_prioritization,
                                                                         start_at=start_at, end_at=end_at, template_count=self.template_count)
                self.results_queue.put([result, start_at, end_at])
                #print('Worker {} added children of {} (ID {}) to results queue'.format(i, smiles, _id))
            except VanillaQueue.Empty:
                #print('Queue {} empty for worker {}'.format(j, i))
                pass
            except Exception as e:
                print e
            # Wait briefly to allow the results_queue to properly update
            time.sleep(0.5)
            self.idle[i] = True

    def initialize(self, reactants, batch_size):
        if batch_size >= self.template_count or batch_size == 0:
            self.expand(reactants, -1, -1)
            # wait until result to continue
            while(self.results_queue.empty()):
                time.sleep(0.5)
        else:
            for start_at in range(0, self.template_count, batch_size):
                self.expand(reactants, start_at, start_at + batch_size)

    def set_template_count(self, template_count):
        self.template_count = template_count

    def evaluate(self, reactants_smiles, contexts, **kwargs):
        self.reset()
        self.nproc = kwargs.pop('nproc', 1)
        batch_size = kwargs.pop('batch_size', 250)
        if not self.celery:
            for i in range(self.nproc):
                self.idle.append(True)
                self.expansion_queue = Queue()

        mol = Chem.MolFromSmiles(reactants_smiles)
        if mol is None: 
            MyLogger.print_and_log('Reactants smiles not parsible: {}'.format(
                    reactants_smiles), template_nn_scorer_loc, level=1)
            return [[{
                        'rank': 1,
                        'outcome': '',
                        'score': 0,
                        'prob': 0,
                        }]]

        clean_reactant_mapping(mol)

        reactants_smiles = Chem.MolToSmiles(mol)
        with allow_join_result():
            self.template_prioritization = kwargs.pop('template_prioritization', gc.popularity)
            self.prepare()
            self.initialize(reactants_smiles, batch_size)
            (all_results, candidate_edits) = self.get_candidate_edits(reactants_smiles)
            reactants = Chem.MolFromSmiles(reactants_smiles)
            atom_desc_dict = edits_to_vectors(
                [], reactants, return_atom_desc_dict=True)
            candidate_tensor = edits_to_tensor(
                candidate_edits, reactants, atom_desc_dict)

            if not candidate_tensor:
                return [[{
                        'rank': 1,
                        'outcome': '',
                        'score': 0,
                        'prob': 0,
                        }]]

            all_outcomes = []
            for context in contexts:
                if context == []:
                    all_outcomes.append({'rank': 0.0,
                                         'outcome': None,
                                         'score': 0.0,
                                         'prob': 0.0,
                                         })
                    continue
                # prediction_context = context_cleaner.clean_context(context) ## move this step to tree evaluator
                prediction_context = context
                context_tensor = context_cleaner.context_to_edit(
                    prediction_context, self.solvent_name_to_smiles, self.solvent_smiles_to_params)
                if not context_tensor: 
                    all_outcomes.append({'rank': 0.0,
                                         'outcome': None,
                                         'score': 0.0,
                                         'prob': 0.0,
                                         })
                    continue
                scores = self.model.predict(candidate_tensor + context_tensor)
                probs = scores

                if kwargs.pop('soft_max', True):
                    probs = softmax(scores)

                this_outcome = sorted(zip(all_results, scores[0], probs[
                                      0]), key=lambda x: x[2], reverse=True)

                # Convert to outcome dict, canonicalizing by SMILES
                outcome_dict = {}
                for i, outcome in enumerate(this_outcome):
                    try:
                        outcome_smiles = outcome[0].smiles
                    except AttributeError:
                        outcome_smiles = outcome[0]['smiles']
                    if outcome_smiles not in outcome_dict:
                        outcome_dict[outcome_smiles] = {
                            'rank': i + 1,
                            'outcome': outcome[0],
                            'score': float(outcome[1]),
                            'prob': float(outcome[2]),
                        }
                    else: # just add probability
                        outcome_dict[outcome_smiles]['prob'] += float(outcome[2])

                all_outcomes.append(sorted(outcome_dict.values(), key=lambda x: x['prob'], reverse=True))

            return all_outcomes

    def get_candidate_edits(self, smiles):

        candidate_edits = []
        stiched_result = ForwardResult(smiles)
        rct_temp = Chem.MolFromSmiles(smiles)
        [a.ClearProp('molAtomMapNumber') for a in rct_temp.GetAtoms()]
        split_smiles = Chem.MolToSmiles(rct_temp).split('.')
        print('SPLIT SMILES FOR GET_CANDIDATE_EDITS: {}'.format(split_smiles))
        all_results = []
        is_ready = [i for (i, res) in enumerate(
            self.pending_results) if res.ready()]
        while self.waiting_for_results():
            try:
                for result, is_ready in self.get_ready_result(is_ready):
                    for product in result.products:
                        if product.smiles not in split_smiles:
                            stiched_result.add_product(product)
                    '''
                    products = result.get_products()
                    for product in products:
                        self.add_product(all_results, candidate_edits, product)
                    '''
                time.sleep(0.5)
                self.pending_results = [res for (i, res) in enumerate(
                    self.pending_results) if i not in is_ready]
                is_ready = [i for (i, res) in enumerate(
                    self.pending_results) if res.ready()]
            except Exception as e:
                print e
                pass

        for product in stiched_result.products:
            all_results.append(product.as_dict())
            candidate_edits.append((product.get_smiles(), product.get_edits()))

        self.stop_expansion()
        return (all_results, candidate_edits)


def softmax(x):
    """Compute softmax values for each sets of scores in x."""
    e_x = np.exp(x - np.max(x))
    return e_x / e_x.sum()
if __name__ == '__main__':
    MyLogger.initialize_logFile()
    celery = False
    db_client = MongoClient(gc.MONGO['path'], gc.MONGO[
                            'id'], connect=gc.MONGO['connect'])
    db = db_client[gc.SOLVENTS['database']]
    SOLVENT_DB = db[gc.SOLVENTS['collection']]
    if not celery:
        forw_trans = ForwardTransformer(mincount=25, celery=False)
        forw_trans.load()
        scorer = TemplateNeuralNetScorer(forward_transformer=forw_trans)
    else:
        scorer = TemplateNeuralNetScorer(celery=celery)

    scorer.load(SOLVENT_DB, gc.PREDICTOR['trained_model_path'])
    res = scorer.evaluate('CN1C2CCC1CC(O)C2.O=C(O)C(CO)c1ccccc1', [
                          [80.0, u'', u'', u'', -1, 50.0]], batch_size=100, nproc=8)
    for re in res[0]:
        print re['outcome'].smiles + " {}".format(re['prob'])
    print 'done'
