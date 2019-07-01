class Chemical(object):
    """Represents a chemical compound.

    Attributes:
        smiles (str): SMILES string of the chemical.
        template_idx_results (dict):
        purchase_price (float or int): Price per gram of the chemical.
            ?? Seems redundant with price ??
        as_reactant (int): Number of times the chemical is seen as a reactant in
            the reaction database.
        as_product (int): Number of times the chemical is seen as a product in
            the reaction database.
        visit_count (int):
        terminal (bool): Whether this node is a terminal node of the tree.
        estimate_price (float or int): Estimated minimum cost.
        estimate_price_sum (float): Sum of estimated costs. Used to calculate
            mean estimate price.
        estimate_price_cnt (int): Count of cost estimates. Used to calculate
            mean estimate price.
        best_template (int): Best template found so far. Gives the lowest price.
        prob (dict or None): Template relevance probabilities to be used for
            later expansion.
        value (None or int): Output of value network, not currently used.
        sorted_id (None or list): Sorted template IDs.
        price (float or int): Valid minimum cost. Negative means not buyable.
        done (bool): Whether this node is done being expanded.
        pathway_count (int): Number of pathways that can be used to make this
            chemical.
    """

    def __init__(self, chem_smi):
        """Initializes Chemical.

        Args:
            chem_smi (str): SMILES string of the chemical.
        """
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
        self.value = None         # output of value network, not currently used
        self.sorted_id = None

        self.price = -1           # valid min cost - means not buyable
        self.done = False

        self.pathway_count = 0

    def __repr__(self):
        return "%s(%r)" % (self.__class__, self.smiles)

    def __str__(self):
        return "%s" % self.smiles

    def set_price(self, value):
        """Sets the price of the chemical.

        Args:
            value (float): Price per gram to set the cost of the chemical to.
        """
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
        """Information about how we might expand this chemical later.

        Args:
            top_probs (list): Highest template relevance probabilities.
            top_indeces (list): Indeces for the corresponding templates.
            value (int): Value of precursor. Currently always 1.
        """
        self.top_probs = top_probs
        self.top_indeces = top_indeces
        self.prob = {top_indeces[i]: top_probs[i] for i in range(len(top_probs))} # dict
        self.sorted_id = top_indeces # copy
        self.value = value
        self.update_estimate_price(value)

    def update_estimate_price(self, value):
        """Updates the current estimated price for the chemical.

        Args:
            value (float or int): The next price estimate to average into
                the current estimate.
        """
        self.estimate_price_sum += value
        self.estimate_price_cnt += 1
        self.estimate_price = self.estimate_price_sum / self.estimate_price_cnt

    def reset(self):
        """Resets the chemical.

        Does nothing.
        """
        return

    # def uniquify_templates(self):
    #     """Looks through all single-step precursors and merges the CTAs that
    #     give identical results"""
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
    """Represents the application of a template to a chemical.

    This is essentially a list-wrapper for the Reaction class because applying
    one template to one product can lead to many distinct product sets.

    Attributes:
        smiles (string): SMILES string of the chemical.
        template_idx (): Index of template to apply.
        waiting (bool):
        vaild (bool):
        reactions (dict of {str: Reaction}):
    """

    def __init__(self, smiles, template_idx):
        self.smiles = smiles.strip()
        self.template_idx = template_idx

        self.waiting = True
        self.valid = True

        self.reactions = {} # key is reactant SMILES string, value is a Reaction


class Reaction(object):
    """Represents a reaction.

    Note: No necessary_reagent or num_examples considered yet.

    Attributes:
        smiles (str): SMILES string of reaction.
        template_idx (): Index of reaction template.
        valid (bool): Whether the reaction is considered valid.
        reactant_smiles (list of str): SMILES strings of reactants.
        visit_count (int):
        done (bool): Whether this node is done being expanded.
        estimate_price (float or int): Estimated minimum cost.
        estimate_price_sum (float): Sum of estimated costs. Used to calculate
            mean estimate price.
        estimate_price_cnt (int): Count of cost estimates. Used to calculate
            mean estimate price.
        price (float or int): Valid minimum cost. Negative means not buyable.
        pathway_count (int): Number of pathways to this reaction.
        filter_score (float):
        tforms (list of ??): Will be list if multiple template_idx possible.
            When merging, the redundant template_idx will be considered INVALID.
        template_score (float): Template relevance score.
        plausibility (float): Fast filter score.
    """

    def __init__(self, smiles, template_idx):
        """Initializes Reaction.

        Args:
            smiles (str): SMILES string of reaction.
            template_idx (): Index of reaction template.
        """
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
        """Updates the current estimated price for the reaction.

        Args:
            value (float or int): The next price estimate to average into
                the current estimate.
        """
        self.estimate_price_sum += value
        self.estimate_price_cnt += 1
        self.estimate_price = self.estimate_price_sum / self.estimate_price_cnt

    def reset(self):
        """Resets the cost of the reaction."""
        self.price = -1
