import makeit.global_config as gc
from multiprocessing import Process, Manager, Queue
from evaluator import Evaluator
import Queue as VanillaQueue
import time
from makeit.utilities.io.logging import MyLogger
from makeit.utilities.io import model_loader
from celery.result import allow_join_result
from askcos_site.askcos_celery.contextrecommender.cr_coordinator import get_context_recommendations
from askcos_site.askcos_celery.treeevaluator.scoring_coordinator import evaluate
from makeit.prioritization.contexts.probability import ProbabilityContextPrioritizer
from makeit.prioritization.contexts.rank import RankContextPrioritizer
from makeit.prioritization.default import DefaultPrioritizer
treeEvaluator_loc = 'tree_evaluator'


class TreeEvaluator():
    '''
    Class for the evaluation of the found retrosynthetic tree
    '''

    def __init__(self, celery=False, rank_inclusion=10, prob_inclusion=0.2, max_contexts=10, single_solv=True,
                 with_smiles=True, context_recommender=None, batch_size=500):
        self.celery = celery
        self.single_solv = single_solv
        self.with_smiles = with_smiles
        self.batch_size = batch_size
        self.evaluator = Evaluator(celery=self.celery)
        self.rank_threshold = rank_inclusion
        self.prob_threshold = prob_inclusion
        self.recommender = context_recommender

        if not self.celery:
            self.context_recommender = model_loader.load_Context_Recommender(
                context_recommender, max_contexts=max_contexts)

        if self.celery:
            def get_contexts(rxn, n):
                res = get_context_recommendations.apply_async(args=(rxn,),
                                                              kwargs={'n': n, 'singleSlvt': self.single_solv,
                                                                      'with_smiles': self.with_smiles,
                                                                      'context_recommender': self.recommender})
                return res.get(120)
        else:
            def get_contexts(rxn, n):
                return self.context_recommender.get_n_conditions(rxn, n=n, singleSlvt=self.single_solv, with_smiles=self.with_smiles)

        self.get_contexts = get_contexts

        if self.celery:
            def evaluate_reaction(reactant_smiles, target, contexts, worker_no = 0):
                res = evaluate.apply_async(args=(reactant_smiles, target, contexts),
                                           kwargs={'mincount': self.mincount, 'forward_scorer': self.forward_scorer})
                return res.get(120)
        else:
            def evaluate_reaction(reactant_smiles, target, contexts, worker_no = 0):
                return self.evaluator.evaluate(reactant_smiles, target, contexts, mincount=self.mincount,
                                               forward_scorer=self.forward_scorer, nproc=self.nproc, 
                                               batch_size=self.batch_size, worker_no = worker_no)

        self.evaluate_reaction = evaluate_reaction

    def get_context_prioritizer(self, context_method):
        if context_method == gc.probability:
            self.context_prioritizer = ProbabilityContextPrioritizer()
        elif context_method == gc.rank:
            self.context_prioritizer = RankContextPrioritizer()
        else:
            MyLogger.print_and_log(
                'Specified prioritization method does not exist. Using default method.', treeEvaluator_loc, level=1)
            self.context_prioritizer = DefaultPrioritizer()
        self.context_prioritizer.load_model()

    def get_top_context(self, evaluation):
        return self.context_prioritizer.get_priority(evaluation)[0]

    def reset(self):
        if self.celery:
            self.evaluation_dict = {}
        else:
            self.evaluation_queue = Queue()
            self.results_queue = Queue()
            self.workers = []
            self.manager = Manager()
            self.done = self.manager.Value('i', 0)
            self.paused = self.manager.Value('i', 0)
            self.idle = self.manager.list()
            self.evaluation_dict = self.manager.dict()
        self.scored_trees = []

    def is_plausible(self, result):
        prob = result['target']['prob']
        rank = result['target']['rank']
        return prob > self.prob_threshold and rank < self.rank_threshold

    def score_step(self, template_score, forward_score):
        if self.tree_scorer == gc.templateonly:
            return template_score
        elif self.tree_scorer == gc.forwardonly:
            return forward_score
        elif self.tree_scorer == gc.product:
            return template_score * forward_score
        else:
            MyLogger.print_and_log(
                'Specified tree scoring method is not implemented. Returning product', treeEvaluator_loc, level=2)
            return template_score * forward_score

    def work(self, i):
        while True:
            # If done, stop
            if self.done.value:
                MyLogger.print_and_log(
                    'Worker {} saw done signal, terminating'.format(i), treeEvaluator_loc)
                break
            # If paused, wait and check again
            if self.paused.value:
                #print('Worker {} saw pause signal, sleeping for 1 second'.format(i))
                time.sleep(1)
                continue
            # Grab something off the queue
            try:
                tree = self.evaluation_queue.get(timeout=0.5)  # short timeout
                self.idle[i] = False
                plausible, score = self.evaluate_tree(tree, context_recommender='', context_scoring_method='',
                                                      forward_scoring_method='', tree_scoring_method='',
                                                      rank_threshold=5, prob_threshold=0.2, is_target=False,
                                                      mincount=25, nproc=1, batch_size=500, n=10, worker_no = i)
                self.results_queue.put([tree, plausible, score])
                #print('Worker {} added children of {} (ID {}) to results queue'.format(i, smiles, _id))
            except VanillaQueue.Empty:
                #print('Queue {} empty for worker {}'.format(j, i))
                pass
            except Exception as e:
                print e
            self.idle[i] = True

    def evaluate_trees(self, tree_list, context_recommender='', context_scoring_method='', forward_scoring_method='',
                       tree_scoring_method='', rank_threshold=5, prob_threshold=0.2, mincount=25, nproc=1,
                       batch_size=500, n=10, nproc_t=1, parallel=False):
        self.reset()
        self.recommender = context_recommender
        self.get_context_prioritizer(context_scoring_method)
        self.rank_threshold = rank_threshold
        self.prob_threshold = prob_threshold
        self.mincount = mincount
        self.nproc = nproc
        self.batch_size = batch_size
        self.forward_scorer = forward_scoring_method
        self.tree_scorer = tree_scoring_method

        if not parallel:
            for tree in tree_list:
                self.scored_trees.append(self.evaluate_tree(tree, context_recommender, context_scoring_method, forward_scoring_method, tree_scoring_method,
                                                            rank_threshold, prob_threshold, mincount, nproc, batch_size, n, is_target=True))
        else:
            self.spin_up_workers(nproc_t)
            self.populate_queue(tree_list)
            self.get_scored_trees()

        return self.scored_trees

    def evaluate_tree(self, tree, context_recommender='', context_scoring_method='', forward_scoring_method='',
                      tree_scoring_method='', rank_threshold=5, prob_threshold=0.2, mincount=25, nproc=1,
                      batch_size=500, n=10, is_target=False, reset=False , worker_no = 0):
        if is_target and reset:
            self.reset()
            self.get_context_prioritizer(context_scoring_method)
            self.rank_threshold = rank_threshold
            self.prob_threshold = prob_threshold
            self.mincount = mincount
            self.recommender = context_recommender
            self.nproc = nproc
            self.batch_size = batch_size
            self.forward_scorer = forward_scoring_method
            self.tree_scorer = tree_scoring_method

        if not tree['children']:
            # Reached the end of the synthesis tree -> Stop
            if is_target:
                return {'tree': tree, 'plausible': True, 'score': 1.0}
            else:
                return True, 1.0
        else:
            with allow_join_result():
                target = tree['smiles']
                reaction = tree['children'][0]
                reactants = [child['smiles'] for child in reaction['children']]
                reaction_smiles = reaction['smiles']
                necessary_reagent = reaction['necessary_reagent']
                ###############################################################
                # If reaction encountered before: get data from dict.
                ###############################################################
                if reaction_smiles in self.evaluation_dict:
                    evaluation = self.evaluation_dict[reaction_smiles]
                ###############################################################
                # Otherwise create data
                ###############################################################
                else:
                    if necessary_reagent:
                        contexts = self.get_contexts(reaction_smiles, 1)
                        reactants.extend(contexts[0][2].split('.')) # add reagents
                    else:
                        contexts = self.get_contexts(reaction_smiles, n)
                    evaluation = self.evaluate_reaction(
                        '.'.join(reactants), target, contexts, worker_no = worker_no)
                    self.evaluation_dict[reaction_smiles] = evaluation
                ###############################################################
                # Process data
                ###############################################################
                if len(evaluation) == 1:
                    top_result = evaluation[0]
                else:
                    top_result = self.get_top_context(evaluation)
                # Add evaluation information to the reaction
                MyLogger.print_and_log('Evaluated reaction: {} - ranked {} with a {}% probability.'
                                       .format(reaction_smiles, top_result['target']['rank'], top_result['target']['prob']*100.0),
                                       treeEvaluator_loc)

                score = self.score_step(
                    reaction['template_score'], top_result['target']['prob'])

                plausible = self.is_plausible(top_result)
                print reaction_smiles, plausible
                all_children_plausible = True
                for child in reaction['children']:
                    child_plausible, score_child = self.evaluate_tree(child)
                    score *= score_child
                    if not child_plausible:
                        all_children_plausible = False

                if all_children_plausible and plausible and is_target:
                    MyLogger.print_and_log(
                        'Found a fully plausible tree!', treeEvaluator_loc)
                elif is_target:
                    MyLogger.print_and_log(
                        'Evaluated tree has unfeasible children.', treeEvaluator_loc)

                reaction['top_product'] = {'smiles': top_result['top_product']['smiles'],
                                           'score': top_result['top_product']['score'],
                                           'prob': top_result['top_product']['prob'],
                                           }
                reaction['forward_score'] = top_result['target']['prob']
                reaction['cumul_score'] = score
                reaction['rank'] = top_result['target']['rank']
                reaction['templates'] = top_result['target']['template_ids']
                reaction['context'] = top_result['context']
                # overwrite
                tree['children'] = [reaction]
                if is_target:
                    return {'tree': tree, 'plausible': plausible and all_children_plausible, 'score': score}
                else:
                    return plausible and all_children_plausible, score

    #############################################################
    # MULTIPROCESSING CODE
    #############################################################
    def spin_up_workers(self, nproc_t):
        self.running = True
        MyLogger.print_and_log('Tree evaluator spinning off {} child processes'.format(
            nproc_t), treeEvaluator_loc)
        for i in range(nproc_t):
            self.idle.append(False)
            p = Process(target=self.work, args=(i,))
            self.workers.append(p)
            p.start()

    def populate_queue(self, tree_list):
        for tree in tree_list:
            self.evaluation_queue.put(tree)

    def waiting_for_results(self):
        time.sleep(0.05)
        waiting = [self.evaluation_queue.empty()]
        waiting.append(self.results_queue.empty())
        waiting.extend(self.idle)
        return (not all(waiting))

    def get_scored_trees(self):
        while self.waiting_for_results():
            try:
                tree, plausible, score = self.results_queue.get(0.2)
                self.scored_trees.append(
                    {'tree': tree, 'plausible': plausible, 'score': score})
            except VanillaQueue.Empty:
                pass
        self.terminate()

    def terminate(self):
        if not self.running:
            return
        self.done.value = 1
        MyLogger.print_and_log(
            'Terminating tree evaluation process.', treeEvaluator_loc)
        time.sleep(1)
        for p in self.workers:
            if p and p.is_alive():
                p.terminate()
        MyLogger.print_and_log(
            'All tree evaluation processes done.', treeEvaluator_loc)
        self.running = False
    ###############################################################
    # MULTIPROCESSING CODE END
    ###############################################################

if __name__ == '__main__':
    from synthetic.context.nearestneighbor import NNContextRecommender
    MyLogger.initialize_logFile()

    ev = TreeEvaluator(context_recommender=gc.nearest_neighbor, celery=False)
    trees = [{'is_chemical': True, 'smiles': 'CN1C2CCC1CC(C2)OC(=O)C(CO)c3ccccc3', 'ppg': 0.0, 'id': 1, 'children': [{'info': '', 'smiles': 'CN1C2CCC1CC(O)C2.O=C(O)C(CO)c1ccccc1>>CN1C2CCC1CC(C2)OC(=O)C(CO)c3ccccc3', 'is_reaction': True, 'num_examples': 19578, 'template_score': 0.017628178000450134, 'children': [
        {'is_chemical': True, 'smiles': 'CN1C2CCC1CC(O)C2', 'ppg': 1.0, 'id': 3, 'children': []}, {'is_chemical': True, 'smiles': 'O=C(O)C(CO)c1ccccc1', 'ppg': 1.0, 'id': 4, 'children': []}], 'id': 2, 'necessary_reagent': u''}]}]
    tree = trees[0]
    res = ev.evaluate_tree(tree, gc.nearest_neighbor, gc.probability,
                           gc.templatebased, gc.product, is_target=True, reset=True, nproc=16)
    #res = ev.evaluate_trees(trees, gc.probability, gc.templatebased, gc.product, nproc = 8, parallel=True, nproc_t=3)
    print res
