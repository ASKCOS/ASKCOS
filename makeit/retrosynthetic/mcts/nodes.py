class Chemical(object):
    """Represents a chemical compound."""

    def __init__(self, chem_smi):
        # Initialize lists for incoming and outgoing reactions.
        self.smiles = chem_smi
        # self.reactions = {} # TODO: deprecate, since it doesn't allow for branching
        self.template_idx_results = {} # key is template_idx, value is a CTA
        self.purchase_price = -1
        self.as_reactant = -1
        self.as_product = -1
        self.visit_count = 0
        self.terminal = False 

        # Counter param used for the DFS search. 
        self.estimate_price = -1         # estimated min cost
        self.estimate_price_sum = 0.
        self.estimate_price_cnt = 0

        self.best_template = -1

        self.prob = None
        self.value = None         # output of vaue network, not currently used
        self.sorted_id = None

        self.price = -1           # valid min cost - means not buyable
        self.done = False

        self.pathway_count = 0

    def __repr__(self):
        return "%s(%r)" % (self.__class__, self.smiles)

    def __str__(self):
        return "%s" % self.smiles

    def set_price(self, value):
        try:
            ppg = float(value)
            # if ppg > 0 and all_cost_1:
            #     ppg = 1.
            self.price = ppg
            self.estimate_price = ppg
            # self.done = True # done isn't determined by price
        except:
            pass

    def set_template_relevance_probs(self, top_probs, top_indeces, value):
        '''Information about how we might expand this chemical later'''
        self.top_probs = top_probs 
        self.top_indeces = top_indeces
        self.prob = {top_indeces[i]: top_probs[i] for i in range(len(top_probs))} # dict
        self.sorted_id = top_indeces # copy
        self.value = value 
        self.update_estimate_price(value)

    def update_estimate_price(self, value):
        self.estimate_price_sum += value
        self.estimate_price_cnt += 1
        self.estimate_price = self.estimate_price_sum / self.estimate_price_cnt

    def reset(self):
        return

    # def uniquify_templates(self):
    #     '''Looks through all single-step precursors and merges the CTAs that
    #     give identical results'''
    #     seen_reactants = {}
    #     template_idx_static = list(self.template_idx_results.keys())
    #     for template_idx in template_idx_static:
    #         CTA = self.template_idx_results[template_idx]
    #         for reactant_smiles in CTA.reactions:
    #             if reactant_smiles not in seen_reactants:
    #                 seen_reactants[reactant_smiles] = template_idx 
    #             else:
    #                 # do the merge
    #                 prev_CTA = self.template_idx_results[seen_reactants[reactant_smiles]]
    #                 prev_CTA.tforms.append(CTA.template_idx)
    #                 prev_CTA.template_score = max(CTA.template_score, prev_CTA.template_score)
    #                 assert prev_CTA.plausibility == CTA.plausibility
    #                 del self.template_idx_results[template_idx]

class ChemicalTemplateApplication(object):
    """Represents the application of a template to a chemical. This is essentially
    a list-wrapper for the Reaction class because applying one template to one
    product can lead to many distinct product sets."""

    def __init__(self, smiles, template_idx):
        self.smiles = smiles.strip()
        self.template_idx = template_idx

        self.waiting = True 
        self.valid = True

        self.reactions = {} # key is reactant SMILES string, value is a Reaction


class Reaction(object):
    """Represents a reaction."""

    def __init__(self, smiles, template_idx):
        """Initialize entry."""
        self.smiles = smiles.strip()
        self.template_idx = template_idx
        # self.depth  = depth 
        self.valid = True
        self.reactant_smiles = []
        self.visit_count = 0

        # self.waiting = True
        self.done = False

        self.estimate_price = -1
        self.estimate_price_sum = 0.
        self.estimate_price_cnt = 0

        self.price = -1

        self.pathway_count = 0

        self.filter_score = 1.0


        # Attributes that will need to be merged if there are two templates that can lead
        # to the same reactant SMILES - will need to check if exists
        self.tforms = [template_idx] # will be list if multiple template_idx possible
        # when merging, the redundant template_idx will be considered INVALID
        self.template_score = 0.0 # template relevance score
        self.plausibility = 0.0 # fast filter score
        # note: no necessary_reagent or num_examples considered yet

    def __repr__(self):
        return "%s(%r)" % (self.__class__, self.smiles)

    def __str__(self):
        return "%s" % self.smiles

    def update_estimate_price(self, value):
        self.estimate_price_sum += value
        self.estimate_price_cnt += 1
        self.estimate_price = self.estimate_price_sum / self.estimate_price_cnt

    def reset(self):
        self.price = -1 

