from makeit.retrosynthetic.transformer import RetroTransformer
from makeit.utilities.buyable.pricer import Pricer
from multiprocessing import Process, Manager, Queue, Pool
from celery.result import allow_join_result
from pymongo import MongoClient

# from makeit.mcts.cost import Reset, score_max_depth, MinCost, BuyablePathwayCount
# from makeit.mcts.misc import get_feature_vec, save_sparse_tree
# from makeit.mcts.misc import value_network_training_states
from makeit.retrosynthetic.mcts.nodes import Chemical, Reaction, ChemicalTemplateApplication
from makeit.utilities.io.logger import MyLogger
from makeit.utilities.io import model_loader
from makeit.utilities.formats import chem_dict, rxn_dict


from collections import defaultdict 

import makeit.global_config as gc
import sys
is_py2 = sys.version[0] == '2'
if is_py2:
    import Queue as VanillaQueue
    import cPickle as pickle
else:
    import queue as VanillaQueue
    import pickle as pickle
import multiprocessing as mp
import numpy as np
import traceback
import itertools
import random
import time 
import gzip 
import sys
import os

os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'

treebuilder_loc = 'mcts_tree_builder'

VIRTUAL_LOSS = 1000000

WAITING = 0
DONE = 1


class MCTS:

    def __init__(self, retroTransformer=None, pricer=None, max_branching=20, max_depth=3, expansion_time=60,
                 celery=False, mincount=25, chiral=True, mincount_chiral=10,
                 template_prioritization=gc.relevance, precursor_prioritization=gc.relevanceheuristic,
                 chemhistorian=None, nproc=4, num_active_pathways=None):
        """Class for retrosynthetic tree expansion using a depth-first search

        Initialization of an object of the TreeBuilder class sets default values
        for various settings and loads transformers as needed (i.e., based on 
        whether Celery is being used or not). Most settings are overridden
        by the get_buyable_paths method anyway.

        Keyword Arguments:
            retroTransformer {None or RetroTransformer} -- RetroTransformer object
                to be used for expansion when *not* using Celery. If none, 
                will be initialized using the model_loader.load_Retro_Transformer
                function (default: {None})
            pricer {Pricer} -- Pricer object to be used for checking stop criteria
                (buyability). If none, will be initialized using default settings
                from the global configuration (default: {None})
            max_branching {number} -- Maximum number of precursor suggestions to
                add to the tree at each expansion (default: {20})
            max_depth {number} -- Maximum number of reactions to allow before
                stopping the recursive expansion down one branch (default: {3})
            expansion_time {number} -- Time (in seconds) to allow for expansion
                before searching the generated tree for buyable pathways (default: {240})
            celery {bool} -- Whether or not Celery is being used. If True, then 
                the TreeBuilder relies on reservable retrotransformer workers
                initialized separately. If False, then retrotransformer workers
                will be spun up using multiprocessing (default: {False})
            nproc {number} -- Number of retrotransformer processes to fork for
                faster expansion (default: {1})
            mincount {number} -- Minimum number of precedents for an achiral template
                for inclusion in the template library. Only used when retrotransformers
                need to be initialized (default: {25})
            mincount_chiral {number} -- Minimum number of precedents for a chiral template
                for inclusion in the template library. Only used when retrotransformers
                need to be initialized. Chiral templates are necessarily more specific,
                so we generally use a lower threshold than achiral templates (default: {10})
            chiral {bool} -- Whether or not to pay close attention to chirality. When 
                False, even achiral templates can lead to accidental inversion of
                chirality in non-reacting parts of the molecule. It is highly 
                recommended to keep this as True (default: {True})
            template_prioritization {string} -- Strategy used for template
                prioritization, as a string. There are a limited number of available
                options - consult the global configuration file for info (default: {gc.popularity})
            precursor_prioritization {string} -- Strategy used for precursor
                prioritization, as a string. There are a limited number of available
                options - consult the global configuration file for info (default: {gc.heuristic})
        """

        self.celery = celery
        self.mincount = mincount
        self.mincount_chiral = mincount_chiral
        self.max_depth = max_depth  
        self.max_branching = max_branching
        self.expansion_time = expansion_time
        self.template_prioritization = template_prioritization
        if self.template_prioritization != gc.relevance:
            raise ValueError('Cannot do MCTS without relevance template prioritization!')
        self.precursor_prioritization = precursor_prioritization
        self.nproc = nproc
        self.chiral = chiral
        self.max_cum_template_prob = 1

        if num_active_pathways is None:
            num_active_pathways = self.nproc
        self.num_active_pathways = num_active_pathways

        ## Pricer
        if pricer:
            self.pricer = pricer
        else:
            self.pricer = Pricer()
            self.pricer.load()

        # Initialize vars, reset dicts, etc.
        self.reset()


        # Get template relevance model - need for target to get things started
        # NOTE: VERY IMPORTANT TO NOT USE TENSORFLOW!! OTHERWISE FORKED PROCESSES HANG
        # THIS SHOULD BE ABLE TO BE FIXED
        from makeit.prioritization.templates.relevance import RelevanceTemplatePrioritizer
        template_prioritizer = RelevanceTemplatePrioritizer(use_tf=False)
        template_prioritizer.load_model()
        self.template_prioritizer = template_prioritizer
        

        # When not using Celery, need to ensure retroTransformer initialized
        if not self.celery:
            if retroTransformer:
                self.retroTransformer = retroTransformer
            else:
                self.retroTransformer = model_loader.load_Retro_Transformer(mincount=self.mincount,
                                                                            mincount_chiral=self.mincount_chiral,
                                                                            chiral=self.chiral)
                self.retroTransformer.load_fast_filter()
                # don't load template prioritizer until later, TF doesn't like forking


        if self.celery:
            def expand(_id, smiles, template_idx): # TODO: make Celery workers
                # Chiral transformation or heuristic prioritization requires
                # same database. _id is _id of active pathway
                self.pending_results.append(tb_c_worker.apply_one_template_by_idx.apply_async(
                    args=(_id, smiles, template_idx),
                    kwargs={'template_count': self.template_count,
                            'max_cum_prob': self.max_cum_template_prob,
                            'apply_fast_filter': self.apply_fast_filter,
                            'filter_threshold': self.filter_threshold},
                    priority=int(depth), # TODO: figure out priority
                    queue=self.private_worker_queue,
                ))
                self.status[(chem_smi, template_idx)] = WAITING
                self.active_pathways_pending[_id] += 1
        else:
            def expand(_id, smiles, template_idx):
                self.expansion_queue.put((_id, smiles, template_idx))
                self.status[(smiles, template_idx)] = WAITING
                self.active_pathways_pending[_id] += 1
        self.expand = expand

        self.status = {}

        # Define method to start up parallelization.
        if self.celery:
            def prepare():
                try:
                    if self.chiral:
                        request = tb_c_worker.reserve_worker_pool.delay()
                        self.private_worker_queue = request.get(timeout=10)
                    else:
                        request = tb_worker.reserve_worker_pool.delay()
                        self.private_worker_queue = request.get(timeout=10)
                except Exception as e:
                    request.revoke()
                    raise IOError(
                        'Did not find an available pool of workers! Try again later ({})'.format(e))
        else:
            def prepare():
                MyLogger.print_and_log('Tree builder spinning off {} child processes'.format(self.nproc), treebuilder_loc)
                for i in range(self.nproc):
                    p = Process(target=self.work, args=(i,))
                    # p.daemon = True
                    self.workers.append(p)
                    p.start()
        self.prepare = prepare

        # Define method to get a processed result.
        if self.celery:
            def get_ready_result():
                # Update which processes are ready
                self.is_ready = [i for (i, res) in enumerate(self.pending_results) if res.ready()]
                for i in self.is_ready:
                    yield self.pending_results[i].get(timeout=0.1)
                    self.pending_results[i].forget()
                self.pending_results = [res for (i, res) in enumerate(self.pending_results) if i not in self.is_ready]
        else:
            def get_ready_result():
                while not self.results_queue.empty():
                    yield self.results_queue.get(timeout=0.5)
        self.get_ready_result = get_ready_result

        # Define how first target is set.
        def set_initial_target(_id, leaves): # i = index of active pathway
            for leaf in leaves:
                if leaf in self.status: # already being worked on
                    continue
                chem_smi, template_idx = leaf
                self.expand(_id, chem_smi, template_idx)     
        self.set_initial_target = set_initial_target

        # Define method to stop working.
        if self.celery:
            def stop():
                pass
        else:
            def stop():
                if not self.running:
                    return
                self.done.value = 1
                #MyLogger.print_and_log('Terminating tree building process.', treebuilder_loc)
                for p in self.workers:
                    if p and p.is_alive():
                        p.terminate()
                #MyLogger.print_and_log('All tree building processes done.', treebuilder_loc)
                self.running = False
        self.stop = stop

    def get_price(self, chem_smi):
        ppg = self.pricer.lookup_smiles(chem_smi, alreadyCanonical=True)
        return ppg
        # if ppg:
        #   return 0.0
        # else:
        #   return None

    def ResetVisitCount(self):
        for chem_key in self.Chemicals: 
            self.Chemicals[chem_key].visit_count = 0
        for rxn_key in self.Reactions: 
            self.Reactions[rxn_key].visit_count = 0
            self.Reactions[rxn_key].successes = []
            self.Reactions[rxn_key].rewards = []


    def coordinate(self):

        while not all(self.initialized):
            MyLogger.print_and_log('Waiting for workers to initialize...', treebuilder_loc)
            time.sleep(2)
        start_time = time.time()
        elapsed_time = time.time() - start_time
        next = 1
        MyLogger.print_and_log('Starting cooridnation loop', treebuilder_loc)
        while (elapsed_time < self.expansion_time): # and self.waiting_for_results():

            if (int(elapsed_time)//5 == next):
                next += 1
                print ("Worked for {}/{} s".format(int(elapsed_time*10)/10.0, self.expansion_time))
                print ("... current min-price {}".format(self.Chemicals[self.smiles].price))
                print ("... |C| = {} |R| = {}".format(len(self.Chemicals), len(self.status)))
                for _id in range(self.num_active_pathways):
                    print('Active pathway {}: {}'.format(_id, self.active_pathways[_id]))
                print('Expansion empty? {}'.format(self.expansion_queue.empty()))
                print('results_queue empty? {}'.format(self.results_queue.empty()))
                print('Active pathway pending? {}'.format(self.active_pathways_pending))
                print('All idle? {}'.format(self.idle))
                # print(self.expansion_queue.qsize()) # TODO: make this Celery compatible
                # print(self.results_queue.qsize())

                # for _id in range(self.nproc):
                #   print(_id, self.expansion_queues[_id].qsize(), self.results_queues[_id].qsize())
                # time.sleep(2)

            for all_outcomes in self.get_ready_result():
                # Record that we've gotten a result for the _id of the active pathway
                _id = all_outcomes[0][0]
                self.active_pathways_pending[_id] -= 1

                # Result of applying one template_idx to one chem_smi can be multiple eoutcomes
                for (_id, chem_smi, template_idx, reactants, filter_score) in all_outcomes:                        
                    # print('coord pulled {} result from result queue'.format(chem_smi))
                    self.status[(chem_smi, template_idx)] = DONE
                    # R = self.Chemicals[chem_smi].reactions[template_idx]
                    C = self.Chemicals[chem_smi]
                    CTA = C.template_idx_results[template_idx] # TODO: make sure CTA created
                    CTA.waiting = False

                    # Any proposed reactants?
                    if len(reactants) == 0: 
                        CTA.valid = False # no precursors, reaction failed
                        continue

                    # Get reactants SMILES
                    reactant_smiles = '.'.join([smi for (smi, _, _, _) in reactants])
                    # TODO: check if banned reaction
                    matched_prev = False
                    for prev_tid, prev_cta in C.template_idx_results.items():
                        if reactant_smiles in prev_cta.reactions: 
                            prev_R = prev_cta.reactions[reactant_smiles]
                            matched_prev = True
                            # Now merge the two...
                            prev_R.tforms.append(template_idx)
                            prev_R.template_score = max(C.prob[template_idx], prev_R.template_score)
                            CTA.reactions[reactant_smiles] = prev_R
                            break
                    if matched_prev:
                        continue # don't make a new reaction

                    # Define reaction using product SMILES, template_idx, and reactants SMILES
                    R = Reaction(chem_smi, template_idx)
                    R.plausibility = filter_score # fast filter score
                    R.template_score = C.prob[template_idx] # template relevance
                    #for smi, prob, value in reactants:
                    for (smi, top_probs, top_indeces, value) in reactants: # all precursors
                        R.reactant_smiles.append(smi)
                        if smi not in self.Chemicals:
                            self.Chemicals[smi] = Chemical(smi)
                            self.Chemicals[smi].set_template_relevance_probs(top_probs, top_indeces, value)
                            ppg = self.get_price(smi)
                            if ppg is not None and ppg > 0:
                                self.Chemicals[smi].set_price(ppg)
                    R.estimate_price = sum([self.Chemicals[smi].estimate_price for smi in R.reactant_smiles])

                    # Add this reaction result to CTA (key = reactant smiles)
                    CTA.reactions[reactant_smiles] = R
                       
            # See if this rollout is done (TODO: make this Celery compatible)
            for _id in range(self.num_active_pathways):
                if self.active_pathways_pending[_id] == 0: # this expansion step is done
                    # This expansion step is done = record!
                    self.update(self.smiles, self.active_pathways[_id])

                    # Set new target
                    leaves, pathway = self.select_leaf()
                    self.active_pathways[_id] = pathway 
                    self.set_initial_target(_id, leaves)


            # if self.expansion_queue.empty() and self.results_queue.empty() and all(self.idle):
            #     # print('All idle and queues empty')
            #     self.update(self.smiles, self.active_pathway)
            #     self.active_pathway = {}

            # Set new target (THIS IS OLD)
            # if len(self.active_pathway) == 0:
            #     # print('Finding a new active pathway')
            #     leaves, pathway = self.select_leaf()
            #     self.active_pathway = pathway
            #     self.set_initial_target(leaves)

            # for _id in range(self.nproc):
            #     if self.expansion_queues[_id].empty() and self.results_queues[_id].empty() and self.idle[_id]:
            #         self.update(self.smiles, self.pathways[_id])
            #         self.pathways[_id] = {}

            # for _id in range(self.nproc):
            #     if len(self.pathways[_id]) == 0:
            #         leaves, pathway = self.select_leaf()
            #         # if len(self.Chemicals) > 30:
            #         # print('###############', _id, leaves, pathway)
            #         self.pathways[_id] = pathway
            #         self.set_initial_target(_id, leaves)

            elapsed_time = time.time() - start_time

            if self.Chemicals[self.smiles].price != -1 and self.time_for_first_path == -1:
                self.time_for_first_path = elapsed_time
                MyLogger.print_and_log('Found the first pathway after {:.2f} seconds'.format(elapsed_time), treebuilder_loc)

        self.stop()

        for _id in range(self.num_active_pathways):
            self.update(self.smiles, self.active_pathways[_id])
        self.active_pathways = [{} for _id in range(self.num_active_pathways)]
        # print(self.active_pathway)

    def work(self, i):
        # with tf.device('/gpu:%d' % (i % self.ngpus)):
        #     self.model = RLModel()
        #     self.model.load(MODEL_PATH)

        # Load models that are required
        self.retroTransformer.get_template_prioritizers(gc.relevance)
        self.initialized[i] = True

        while True:
            # If done, stop
            if self.done.value:
                # print 'Worker {} saw done signal, terminating'.format(i)
                break
            
            # Grab something off the queue
            if not self.expansion_queue.empty():
                try:
                    self.idle[i] = False
                    (_id, smiles, template_idx) = self.expansion_queue.get(timeout=0.1)  # short timeout
                
                    # print('{} grabbed {} and {} from queue'.format(_id, smiles, template_idx))
                    try:
                        all_outcomes = self.retroTransformer.apply_one_template_by_idx(_id, smiles, template_idx) # TODO: add settings
                    except Exception as e:
                        print(e)
                        all_outcomes = [(_id, smiles, template_idx, [], 0.0)]
                    # print('{} applied one template and got {}'.format(i, all_outcomes))
                    # all_outcomes = list of (_id, smiles, template_idx, reactants, filter_score)
               
                    
                    self.results_queue.put(all_outcomes)
                    # print('{} put {} outcomes on queue'.format(i, len(all_outcomes)))

                except VanillaQueue.Empty:
                    self.idle[i] = True
                    pass # looks like someone got there first...

            
            # time.sleep(0.01)
            self.idle[i] = True

    def UCB(self, chem_smi, c_exploration=0.2, path=[]):
        '''
        Can either select an unapplied template to apply, or select a specific reactant to expand further (?)
        TODO: check these changes...
        '''
        rxn_scores = []

        C = self.Chemicals[chem_smi]
        product_visits = C.visit_count
        max_estimate_price = 0

        for template_idx in C.template_idx_results:
            CTA = C.template_idx_results[template_idx] 
            if CTA.waiting or not CTA.valid:
                continue

            for reactants_smi in CTA.reactions:
                R = CTA.reactions[reactants_smi]

                if len(set(R.reactant_smiles) & set(path)) > 0: # avoid cycles
                    continue
                if R.done:
                    continue
                max_estimate_price = max(max_estimate_price, R.estimate_price)
                Q_sa = - R.estimate_price
                try:
                    U_sa = c_exploration * C.prob[template_idx] * np.sqrt(product_visits) / (1 + R.visit_count)
                except:
                    print(chem_smi, product_visits)
                score = Q_sa + U_sa
                rxn_scores.append((score, template_idx, reactants_smi))

        # unexpanded template - find most relevant template that hasn't been tried
        num_branches = len(rxn_scores)
        if num_branches < self.max_branching or chem_smi == self.smiles:
            for template_idx in C.top_indeces: 
                if template_idx not in C.template_idx_results:
                    Q_sa = - (max_estimate_price + 0.1)
                    U_sa = c_exploration * C.prob[template_idx] * np.sqrt(product_visits) / 1
                    score = Q_sa + U_sa
                    rxn_scores.append((score, template_idx, None)) # record estimated score if we were to actually apply that template
                    # TODO: figure out if this "None" makes sense for the reactants smiles 
                    break

        if len(rxn_scores) > 0:
            sorted_rxn_scores = sorted(rxn_scores, key=lambda x: x[0], reverse=True)
            best_rxn_score, selected_template_idx, selected_reactants_smi = sorted_rxn_scores[0] # get next best template to apply
        else:
            selected_template_idx, selected_reactants_smi = None, None

        return selected_template_idx, selected_reactants_smi


    def select_leaf(self, c_exploration=1.):
        #start_time = time.time()
        pathway = {}
        leaves = []
        queue = VanillaQueue.Queue()
        queue.put((self.smiles, 0, [self.smiles]))

        while not queue.empty():
            chem_smi, depth, path = queue.get()
            if depth >= self.max_depth or chem_smi in pathway: # don't go too deep or recursively
                continue
            template_idx, reactants_smi = self.UCB(chem_smi, c_exploration=c_exploration, path=path)
            if template_idx is None:
                continue
            
            # Only grow pathway when we have picked a specific reactants_smi (?)
            if reactants_smi is not None:
                pathway[chem_smi] = (template_idx, reactants_smi) # TODO: figure out if reactants_smi==None case is an issue
            else:
                pathway[chem_smi] = template_idx # still record template selection
            
            C = self.Chemicals[chem_smi]
            C.visit_count += VIRTUAL_LOSS

            # print('Looking at chemical C: {}'.format(C))
            if template_idx not in C.template_idx_results:
                # print('Creating CTA for {} and {}'.format(chem_smi, template_idx))
                C.template_idx_results[template_idx] = ChemicalTemplateApplication(chem_smi, template_idx)
                CTA = C.template_idx_results[template_idx]

                # TODO: figure out VIRTUAL_LOSS for R.visit_count change?
                # C.reactions[template_idx] = Reaction(chem_smi, template_idx)
                # R = C.reactions[template_idx]
                # R.visit_count += VIRTUAL_LOSS
                leaves.append((chem_smi, template_idx))

            else:
                # Can we assume that the reactants_smi exists in this CTA? I guess so...
                CTA = C.template_idx_results[template_idx]

                if reactants_smi: # if we choose a specific reaction, not just a template...
                    if reactants_smi in CTA.reactions:

                        R = CTA.reactions[reactants_smi]
                        R.visit_count += VIRTUAL_LOSS

                        for smi in R.reactant_smiles:
                            assert smi in self.Chemicals
                            # if self.Chemicals[smi].purchase_price == -1:
                            if not self.Chemicals[smi].done:
                                queue.put((smi, depth+1, path+[smi]))
                        if R.done:
                            C.visit_count += R.visit_count
                            R.visit_count += R.visit_count

        return leaves, pathway


    def update(self, chem_smi, pathway, depth=0):
        
        if depth == 0:
            for smi in pathway:
                if type(pathway[smi]) == tuple:
                    (template_idx, reactants_smi) = pathway[smi]
                else:
                    (template_idx, reactants_smi) = (pathway[smi], None)
                C = self.Chemicals[smi]
                CTA = C.template_idx_results[template_idx]
                C.visit_count -= (VIRTUAL_LOSS - 1)
                if reactants_smi:
                    R = CTA.reactions[reactants_smi]
                    R.visit_count -= (VIRTUAL_LOSS - 1)

        if (chem_smi not in pathway) or (depth >= self.max_depth):
            return

        if type(pathway[chem_smi]) == tuple:
            (template_idx, reactants_smi) = pathway[chem_smi]
        else:
            (template_idx, reactants_smi) = (pathway[chem_smi], None)

        C = self.Chemicals[chem_smi]
        CTA = C.template_idx_results[template_idx]
        if CTA.waiting: # haven't actually expanded
            return 

        if reactants_smi:
            R = CTA.reactions[reactants_smi]
            if R.valid and (not R.done):
                R.done = all([self.Chemicals[smi].done for smi in R.reactant_smiles])

                for smi in R.reactant_smiles:
                    self.update(smi, pathway, depth+1)
                
                estimate_price = sum([self.Chemicals[smi].estimate_price for smi in R.reactant_smiles])
                R.update_estimate_price(estimate_price)
                C.update_estimate_price(estimate_price)

                price_list = [self.Chemicals[smi].price for smi in R.reactant_smiles]
                if all([price != -1 for price in price_list]):
                    price = sum(price_list)
                    R.price = price
                    if R.price < C.price or C.price == -1:
                        C.price = R.price

        if sum(len(CTA.reactions) for tid,CTA in C.template_idx_results.items()) >= self.max_branching:
            print('{} hit max branching, checking if "done"'.format(chem_smi))
            C.done = all([(R.done or (not R.valid)) for rsmi,R in CTA.reactions.items() for tid,CTA in C.template_idx_results.items()])

        # if C.price != -1 and C.price < C.estimate_price:
        #   C.estimate_price = C.price


    def full_update(self, chem_smi, depth=0, path=[]):

        C = self.Chemicals[chem_smi]
        C.pathway_count = 0

        if C.purchase_price != -1:
            C.pathway_count = 1
            return

        if depth > self.max_depth:
            return

        prefix = '    '* depth

        for template_idx in C.template_idx_results:
            CTA = C.template_idx_results[template_idx]
            for reactants_smi in CTA.reactions:
                R = CTA.reactions[reactants_smi]
                R.pathway_count = 0
                if (not R.valid) or len(set(R.reactant_smiles) & set(path)) > 0:
                    continue
                for smi in R.reactant_smiles:
                    self.full_update(smi, depth+1, path+[chem_smi])
                price_list = [self.Chemicals[smi].price for smi in R.reactant_smiles]
                if all([price != -1 for price in price_list]):
                    price = sum(price_list)
                    R.price = price
                    if R.price < C.price or C.price == -1:
                        C.price = R.price
                        C.best_template = template_idx
                    R.pathway_count = np.prod([self.Chemicals[smi].pathway_count for smi in R.reactant_smiles])
                    # if R.pathway_count != 0:
                    #   print(prefix + '  Reac %d: '%template_idx + str(R.reactant_smiles) + ' %d paths'%R.pathway_count)
                else:
                    R.pathway_count = 0

                # print(prefix + str(R.reactant_smiles) + ' - %d' % R.pathway_count)

        C.pathway_count = 0
        for tid,CTA in C.template_idx_results.items():
            for rct_smi,R in CTA.reactions.items():
                C.pathway_count += R.pathway_count

        # if C.pathway_count != 0:
        #   print(prefix + chem_smi + ' %d paths, price: %.1f' % (C.pathway_count, C.price))


    def build_tree(self):

        self.running = True

        # Reset dicts
        Chemicals, Reactions = {}, {}
        self.Chemicals = Chemicals
        self.Reactions = Reactions
        
        # Define first chemical node (target)
        max_cum_prob = 0.995 # TODO: make setting
        template_count = 100 # TODO: make setting
        probs, indeces = self.template_prioritizer.get_topk_from_smi(self.smiles, k=template_count)
        truncate_to = np.argwhere(np.cumsum(probs) >= max_cum_prob)
        if len(truncate_to):
            truncate_to = truncate_to[0][0] + 1 # Truncate based on max_cum_prob?
        else:
            truncate_to = template_count
        value = 1 # current value assigned to precursor (note: may replace with real value function)
        self.Chemicals[self.smiles] = Chemical(self.smiles)
        self.Chemicals[self.smiles].set_template_relevance_probs(probs[:truncate_to], indeces[:truncate_to], value)
        MyLogger.print_and_log('Calculating initial probs for target', treebuilder_loc)

        # First selection is all the same
        leaves, pathway = self.select_leaf()
        for _id in range(self.num_active_pathways):
            self.active_pathways[_id] = pathway
            self.set_initial_target(_id, leaves)
        MyLogger.print_and_log('Set initial leaves for active pathways', treebuilder_loc)
        
        # Coordinate workers.
        self.prepare()
        self.coordinate()

        # Do a final pass to get counts
        MyLogger.print_and_log('Doing final update of pathway counts / prices', treebuilder_loc)
        self.full_update(self.smiles)
        C = self.Chemicals[self.smiles]

        print("Finished working.")
        print("=== found %d pathways (overcounting duplicate templates)" % C.pathway_count)
        print("=== time for fist pathway: %.2fs" % self.time_for_first_path)
        print("=== min price: %.1f" % C.price)
        print("---------------------------")
        return # self.Chemicals, C.pathway_count, self.time_for_first_path


    def tree_status(self):
        """Summarize size of tree after expansion

        Returns:
            num_chemicals {int} -- number of chemical nodes in the tree
            num_reactions {int} -- number of reaction nodes in the tree
        """

        num_chemicals = len(self.Chemicals)
        num_reactions = len(self.status)
        return (num_chemicals, num_reactions, [])


    def reset(self):
        '''Prepare for a new expansion'''
        if self.celery:
            # general parameters in celery format
            # TODO: anything goes here?
            pass
        else:
            self.workers = []
            self.running = False
            self.manager = Manager()
            self.done = self.manager.Value('i', 0)
            self.idle = self.manager.list()
            self.initialized = self.manager.list()
            for i in range(self.nproc):
                self.idle.append(True)
                self.initialized.append(False)
            self.expansion_queue = Queue()
            self.results_queue = Queue()
            
        self.status = {}
        self.active_pathways = [{} for _id in range(self.num_active_pathways)]
        self.active_pathways_pending = [0 for _id in range(self.num_active_pathways)]
        self.pathway_count = 0 
        self.mincost = 10000.0        
        self.Chemicals = {} # new
        self.Reactions = {} # new
        self.time_for_first_path = -1


    def return_trees(self):

        def IDDFS():
            """Perform an iterative deepening depth-first search to find buyable
            pathways.
                        
            Yields:
                nested dictionaries defining synthesis trees
            """
            for path in DLS_chem(self.smiles, depth=0, headNode=True):
                yield chem_dict(self.smiles, children=path, **{})


        def DLS_chem(chem_smi, depth, headNode=False):
            """Expand at a fixed depth for the current node chem_id."""
            C = self.Chemicals[chem_smi]
            if C.purchase_price != -1:
                yield []

            if depth > self.max_depth:
                return

            done_children_of_this_chemical = []
            for tid, CTA in C.template_idx_results.items():
                if CTA.waiting:
                    continue
                for rct_smi, R in CTA.reactions.items():
                    if (not R.valid) or R.price == -1:
                        continue
                    rxn_smiles = '.'.join(sorted(R.reactant_smiles)) + '>>' + chem_smi
                    if rxn_smiles not in done_children_of_this_chemical: # necessary to avoid duplicates
                        for path in DLS_rxn(chem_smi, tid, rct_smi, depth):
                            yield [rxn_dict(tid, rxn_smiles, children=path, **{})]
                        done_children_of_this_chemical.append(rxn_smiles)


        def DLS_rxn(chem_smi, template_idx, rct_smi, depth):
            """Return children paths starting from a specific rxn_id"""
            # TODO: add in auxiliary information about templates, etc.
            R = self.Chemicals[chem_smi].template_idx_results[template_idx].reactions[rct_smi]

            # rxn_list = []
            # for smi in R.reactant_smiles:
            #     rxn_list.append([chem_dict(smi, children=path, **{}) for path in DLS_chem(smi, depth+1)])
                
            # return [rxns[0] for rxns in itertools.product(rxn_list)]

            ###################
            # To get recursion working properly with generators, need to hard-code these cases? Unclear
            # whether itertools.product can actually work with generators. Seems like it can't work
            # well...

            # Only one reactant? easy!
            if len(R.reactant_smiles) == 1:
                chem_smi0 = R.reactant_smiles[0]
                for path in DLS_chem(chem_smi0, depth+1):
                    yield [
                        chem_dict(chem_smi0, children=path, **{})
                    ]

            # Two reactants? want to capture all combinations of each node's
            # options
            elif len(R.reactant_smiles) == 2:
                chem_smi0 = R.reactant_smiles[0]
                chem_smi1 = R.reactant_smiles[1]
                for path0 in DLS_chem(chem_smi0, depth+1):
                    for path1 in DLS_chem(chem_smi1, depth+1):
                        yield [
                            chem_dict(chem_smi0, children=path0, **{}),
                            chem_dict(chem_smi1, children=path1, **{}),
                        ]

            # Three reactants? This is not elegant...
            elif len(R.reactant_smiles) == 3:
                chem_smi0 = R.reactant_smiles[0]
                chem_smi1 = R.reactant_smiles[1]
                chem_smi2 = R.reactant_smiles[2]
                for path0 in DLS_chem(chem_smi0, depth+1):
                    for path1 in DLS_chem(chem_smi1, depth+1):
                        for path2 in DLS_chem(chem_smi2, depth+1):
                            yield [
                                chem_dict(chem_smi0, children=path0, **{}),
                                chem_dict(chem_smi1, children=path1, **{}),
                                chem_dict(chem_smi2, children=path2, **{}),
                            ]

            # I am ashamed
            elif len(R.reactant_smiles) == 4:
                chem_smi0 = R.reactant_smiles[0]
                chem_smi1 = R.reactant_smiles[1]
                chem_smi2 = R.reactant_smiles[2]
                chem_smi3 = R.reactant_smiles[3]
                for path0 in DLS_chem(chem_smi0, depth+1):
                    for path1 in DLS_chem(chem_smi1, depth+1):
                        for path2 in DLS_chem(chem_smi2, depth+1):
                            for path3 in DLS_chem(chem_smi3, depth+1):
                                yield [
                                    chem_dict(chem_smi0, children=path0, **{}),
                                    chem_dict(chem_smi1, children=path1, **{}),
                                    chem_dict(chem_smi2, children=path2, **{}),
                                    chem_dict(chem_smi3, children=path3, **{}),
                                ]

            else:
                print('Too many reactants! Only have cases 1-4 programmed')
                print('There probably are not any real 5 component reactions')
                print(R.reactant_smiles)


        MyLogger.print_and_log('Retrieving trees...', treebuilder_loc)
        trees = [tree for tree in IDDFS()]

        # Sort by some metric
        def number_of_starting_materials(tree):
            if tree != []:
                if tree['children']:
                    return sum(number_of_starting_materials(tree_child) for tree_child in tree['children'][0]['children'])
            return 1.0
        def number_of_reactions(tree):
            if tree != []:
                if tree['children']:
                    return 1.0 + max(number_of_reactions(tree_child) for tree_child in tree['children'][0]['children'])
            return 0.0

        MyLogger.print_and_log('Sorting {} trees...'.format(len(trees)), treebuilder_loc)
        trees = sorted(trees, key=lambda x: number_of_reactions(x))

        return self.tree_status(), trees 


    # TODO: use these settings...
    def get_buyable_paths(self, 
                            smiles, 
                            max_depth=10,
                            max_branching=25,
                            expansion_time=120,
                            nproc=4,
                            num_active_pathways=None,
                            chiral=True,
                            max_trees=5000,
                            max_ppg=1e10,
                            known_bad_reactions=[],
                            forbidden_molecules=[],
                            template_count=100,
                            max_cum_template_prob=0.995, 
                            max_natom_dict=defaultdict(lambda: 1e9, {'logic': None}),
                            min_chemical_history_dict={'as_reactant':1e9, 'as_product':1e9,'logic':None},
                            apply_fast_filter= True, 
                            filter_threshold=0.75,
                            **kwargs):

        

        self.smiles = smiles 
        self.max_depth = max_depth
        self.expansion_time = expansion_time
        self.nproc = nproc
        if num_active_pathways is None:
            num_active_pathways = nproc
        self.num_active_pathways = num_active_pathways

        self.reset()
        
        MyLogger.print_and_log('Starting search for {}'.format(smiles), treebuilder_loc)
        self.build_tree()

        return self.return_trees()


if __name__ == '__main__':

    random.seed(1)
    np.random.seed(1)
    MyLogger.initialize_logFile()
    simulation_time = 30.

    # Load tree builder 
    NCPUS = 4
    print("There are {} processes available ... ".format(NCPUS))
    Tree = MCTS(nproc=NCPUS, mincount=gc.RETRO_TRANSFORMS_CHIRAL['mincount'], 
        mincount_chiral=gc.RETRO_TRANSFORMS_CHIRAL['mincount_chiral'])


    ####################################################################################
    ############################# DEBUGGING ############################################
    ####################################################################################

    smiles = "C1=CC(=C(C=C1F)F)C(CN2C=NC=N2)(CN3C=NC=N3)O"
    import rdkit.Chem as Chem 
    smiles = Chem.MolToSmiles(Chem.MolFromSmiles(smiles), True)
    status, paths = Tree.get_buyable_paths(smiles,
                                        nproc=NCPUS,
                                        expansion_time=simulation_time)
    print(status)
    for path in paths[:5]:
        print(path)
    print('Total num paths: {}'.format(len(paths)))

    ####################################################################################
    ############################# TESTING ##############################################
    ####################################################################################

    f = open('test_smiles.txt')
    N = 5
    smiles_list = [line.strip() for line in f]

    ########### STAGE 1 - PROCESS ALL CHEMICALS
    with open('chemicals.pkl', 'wb') as fid:
        for _id, smiles in enumerate(smiles_list[:N]): 
            smiles = Chem.MolToSmiles(Chem.MolFromSmiles(smiles), True)
            status, paths = Tree.get_buyable_paths(smiles,
                                                nproc=NCPUS,
                                                expansion_time=simulation_time)
    
            pickle.dump((Tree.Chemicals, Tree.time_for_first_path, paths), fid)

    ########### STAGE 2 - ANALYZE RESULTS
    success = 0
    total = 0
    first_time = []
    pathway_count = []
    min_price = []
    with open('chemicals.pkl', 'rb') as fid:
        for _id, smiles in enumerate(smiles_list[:N]): 
            smiles = Chem.MolToSmiles(Chem.MolFromSmiles(smiles), True)
            (Chemicals, ftime, paths) = pickle.load(fid)

            total += 1
            if Chemicals[smiles].price != -1:
                success += 1
                first_time.append(ftime)
                pathway_count.append(len(paths))
                min_price.append(Chemicals[smiles].price)

        print('After looking at chemical index {}'.format(_id))
        print('Success ratio: %f (%d/%d)' % (float(success)/total), success, total)        
        print('average time for first pathway: %f' % np.mean(first_time))
        print('average number of pathways:     %f' % np.mean(pathway_count))
        print('average minimum price:          %f' % np.mean(min_price))
    
