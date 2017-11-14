import global_config as gc
from retro_transformer import RetroTransformer
from heuristic_prioritizer import HeuristicPrioritizer
from utilities.buyable.pricer import Pricer
from utilities.i_o.logging import MyLogger
treebuilder_loc = 'tree_builder'

class TreeBuilder:
    '''
    Class for retro-synthetic tree expansion. The expansion is performed based on the retroTransformer that is provided.
    '''
    def __init__(self, retroTransformer = None, pricer = None, max_branching = 20, tree_dict = None, current_id = 2):
        if retroTransformer:
            self.retroTransformer = retroTransformer
        else:
            MyLogger.print_and_log('Cannot build a tree without a transformer. Exiting...', treebuilder_loc, level = 3)
    
        if pricer:
            self.pricer = pricer
        else:
            MyLogger.print_and_log('Cannot build a buyable tree without a pricer. Exiting...', treebuilder_loc, level = 3)
            
        if not tree_dict == None:
            self.tree_dict = tree_dict
        else:
            MyLogger.print_and_log('Cannot build a tree without a tree dictionary. Exiting...', treebuilder_loc, level = 3)
        
        self.max_branching = max_branching
        self.current_id = current_id
        
    def get_children(self, smiles):
        '''
        Find children for the desired molecule (in smiles format). Should be used in the "Worker" in parallel processing
        '''
        result = self.retroTransformer.perform_transformations(smiles)
        precursors = result.return_top(n=self.max_branching)
        children = []
        for precursor in precursors:
            children.append((
                {'template': precursor['tforms'][0],
                'necessary_reagent': precursor['necessary_reagent'],
                'num_examples': precursor['num_examples'],
                'score': precursor['score'],
                },
                precursor['smiles_split']
                ))

        return children

    def add_children(self, children, unique_id):
        '''
        Add the precursors to the dictionary the builder was initiated with. Should be used in the "Coordinator" in
        in parallel processing.
        To ensure correct working of the indexation, only one treebuilder object should be instantiated.
        '''
        parent_chem_doc = self.tree_dict[_id] # copy to overwrite later
        parent_chem_prod_of = parent_chem_doc['prod_of']

        # Assign unique number
        for (rxn, mols) in children:
            rxn_id = 0
            #depending on whether current_id was given as 'Manager.Value' type or 'Integer':
            try:
                rxn_id = self.current_id.value
                self.current_id.value += 1 # this is only okay because there is/should be only ONE treebuilder
            except AttributeError:
                rxn_id = self.current_id
                self.current_id += 1
                
                
            # For the parent molecule, record child reactions
            parent_chem_prod_of.append(rxn_id)
        
            # For the reaction, keep track of children IDs
            chem_ids = []
            for mol in mols:
                
                # New chemical?
                if mol not in chem_to_id:

                    try:
                        chem_id = self.current_id.value
                        self.current_id.value += 1 # this is only okay because there is/should be only ONE treebuilder
                    except AttributeError:
                        chem_id = self.current_id
                        self.current_id += 1

                    # Check if buyable
                    ppg = pricer.lookup_smiles(mol, alreadyCanonical = True)

                    tree_dict[chem_id] = {
                        'smiles': mol,
                        'prod_of': [],
                        'rct_of': [rxn_id],
                        'depth': parent_chem_doc['depth'] + 1,
                        'ppg': ppg
                    }
                    chem_to_id[mol] = chem_id

                    if ppg:
                        #print('{} buyable!'.format(mol))
                        buyable_leaves.append(chem_id)

                    else:
                        # Add to queue to get expanded
                        try:
                            expansion_queues[parent_chem_doc['depth']].put((chem_id, mol))
                            added_to_queue += 1
                        except IndexError:
                            # Maximum depth reached
                            if gc.DEBUG:
                                MyLogger.print_and_log('Reached maximum depth, so will not expand around {}'.format(tree_dict[chem_id]), treebuilder_loc)
                            pass
                else:
                    chem_id = chem_to_id[mol]

                    # Overwrite this chemical node to record it is a reactant of this rxn
                    chem_doc = tree_dict[chem_id]
                    chem_doc['rct_of'] += [rxn_id]
                    tree_dict[chem_id] = chem_doc

                # Save ID
                chem_ids.append(chem_id)

            # Record by overwriting the whole dict value
            rxn['rcts'] = chem_ids
            rxn['prod'] = _id
            rxn['depth'] = parent_chem_doc['depth'] + 0.5
            tree_dict[rxn_id] = rxn

        # Overwrite dictionary entry for the parent
        parent_chem_doc['prod_of'] = parent_chem_prod_of
        tree_dict[_id] = parent_chem_doc
        
if __name__ == '__main__':
    
    MyLogger.initialize_logFile()
    prioritizer = HeuristicPrioritizer()
    retroTrans = RetroTransformer(prioritizer = prioritizer, mincount = 4)
    retroTrans.load()
    pricer = Pricer()
    pricer.load()
    treedict = []
    treeBuilder = TreeBuilder(retroTransformer=retroTrans, pricer=  pricer, tree_dict=treedict)
    print(treeBuilder.get_children('c1ccccc1C(=O)OCCN'))
    
        