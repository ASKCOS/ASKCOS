from django.db import models
from django.contrib.auth.models import User

class SavedResults(models.Model):
    user = models.ForeignKey(User)
    description = models.CharField(max_length=200)
    created = models.DateTimeField('saved on')
    fpath = models.CharField(max_length=500)


class SeparationInput(models.Model):
    tgt_smiles = models.CharField("Enter SMILES string of the target molecule", max_length=200)
    other_smiles = models.CharField('Enter SMILES strings of other solute molecules separated with comma',
                                   max_length=800)

    # Separation conditions
    org1 = models.CharField("Enter SMILES string of the organic solvent", max_length=200)
    org2 = models.CharField("Enter SMILES string of a second organic solvent and its molar fraction "
                            "separated with comma, otherwise 'None", max_length=200, blank=True)
    pH = models.FloatField("Enter the pH of the aqueous phase or blank", blank=True)
    ionic_strength = models.FloatField("Enter the ionic strength of the aqueous phase or blank", blank=True)

    temp = models.FloatField("Enter the separation temperature in Celsius degree or blank", blank=True)
    flow_ratio = models.FloatField("Enter the org/aqu flow rate ratio or blank", blank=True)

    # Optimization parameters
    OPT_CHOICES = (('h', 'pH_only'), ('f', 'flow_only'), ('b', 'both'))
    opt_option = models.CharField("What to vary: pH, flow or both:", max_length=1, choices=OPT_CHOICES)
    pH_step = models.FloatField('If to vary pH, enter the pH increment for optimization or blank', blank=True)
    flow_ratio_min = models.FloatField("If to vary org/aqu flow rate ratio, enter the minimum or blank", blank=True)
    flow_ratio_max = models.FloatField("If to vary org/aqu flow rate ratio, enter the maximum or blank", blank=True)
    num_flowratio = models.IntegerField("The number of flowratios in the range or blank", blank=True)

    PHASE_CHOICES = (('aqu', 'aqueous'), ('org', 'organic'), ('e', 'either'))
    max_phase = models.CharField("In which phase the target molecule should be concentrated:",
                                 max_length=3, choices=PHASE_CHOICES)
    min_phase = models.CharField("In which phase the target molecule should be minimized or leave blank:",
                                 max_length=3, choices=PHASE_CHOICES, blank=True)

    OPT_TGT_CHOICES = (('y', 'yield'), ('p', 'purity'), ('b', 'both'))
    opt_tgt_choice = models.CharField('Select the optimization target(s) of interest',
                                      max_length=1, choices=OPT_TGT_CHOICES)
    aqu_conc = models.CharField("If purity needed, enter the concentration for each molecule in the aqueous phase "
                                "and separate with comma", max_length=200)
    org_conc = models.CharField("If purity needed, enter the concentration for each molecule in the organic phase "
                                "and separate with comma", max_length=200)

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