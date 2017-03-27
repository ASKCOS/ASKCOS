from django.db import models

class SeparationInput(models.Model):
    # Molecules
    tgt_smiles = models.CharField("Enter SMILES string of the target molecule", max_length=200)
    other_smiles = models.CharField('Enter SMILES strings of other solute molecules separated with comma',
                                   max_length=800)

    # Separation conditions
    org1 = models.CharField("Enter SMILES string of the organic solvent", max_length=200)
    org2 = models.CharField("Enter SMILES string of a second organic solvent and its molar fraction "
                            "separated with comma, otherwise 'None", max_length=200, blank=True)
    pH = models.FloatField("Enter the pH of the aqueous phase or 'None'", blank=True)
    ionic_strength = models.FloatField("Enter the ionic strength of the aqueous phase or 'None", blank=True)

    temp = models.FloatField("Enter the separation temperature in Celsius degree or 'None'", blank=True)
    flow_ratio = models.FloatField("Enter the org/aqu flow rate ratio or 'None'", blank=True)

    # Optimization parameters
    pH_step = models.FloatField('Enter the pH increment for optimization or None', blank=True)
    flow_ratio_min = models.FloatField("The min org/aqu flow rate ratio or 'None'", blank=True)
    flow_ratio_max = models.FloatField("The max org/aqu flow rate ratio or 'None'", blank=True)

    PHASE_CHOICES = [('a', 'aqueous'), ('o', 'organic'), ('e','either')]
    max_phase = models.CharField("In which phase the target molecule should be concentrated:",
                                 max_length=1, choices=PHASE_CHOICES)
    min_phase = models.CharField("In which phase the target molecule should be minimized:",
                                 max_length=1, choices=PHASE_CHOICES, blank=True)

    OPT_TGT_CHOICES = [('y', 'yield'), ('p', 'purity'), ('b', 'both')]
    opt_tgt_choice = models.CharField('Select the optimization target(s) of interest',
                                      max_length=1, choices=OPT_TGT_CHOICES)
    if opt_tgt_choice is ('p' or 'b'):
        aqu_conc = models.CharField("Enter the concentration for each molecule or 'None' "
                                         "in the aqueous phase and separate with comma", max_length=200, blank=True)
        org_conc = models.CharField("Enter the concentration for each molecule or 'None' "
                                         "in the organic phase and separate with comma", max_length=200, blank=True)

    ##### TO DO:
        # solvent selection option
        # Score

    # def is_valid_pH(self):
    #     # Need to check whether this limit does exist
    #     if self.pH <= 14.0 and self.pH >= 0.0:
    #         return True
    #
    # def is_valid_pH_step(self):
    #     if self.pH_step < 14.0 and self.pH_step > 0.0:
    #         return True


class ConditionPredictorSetup(models.Model):

    # Optimization parameters
    N_conditions = [('1', 'Top 1'), ('2', 'Top 1 + one more')]
    num_cond = models.CharField("Number of recommendations required, no more than 2:", max_length=1, choices=N_conditions)
    dist_limit = models.FloatField('Cosine distance limit cutoff for Nearest Neighbor in [0, 2], default = 0.3', blank=True)
    outputString = models.BooleanField("Don't output detailed reaction conditions, default")

    ##### Potentially the NN model and the list of corresponding reaction IDs can be varied as well