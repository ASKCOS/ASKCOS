import numpy as np
import requests

class TFServingAPIModel(object):
    """Template relevance prioritization model served using TF serving.
    
    Attributes:
        hostname (str): hostname of service serving tf model.
        model_name (str): Name of model provided to tf serving.
        fp_length (int): Fingerprint length.
        fp_radius (int): Fingerprint radius.

    """
    def __init__(self, hostname, model_name, version=1):
        self.url = 'http://{}:8501/v{}/models/{}:predict'.format(
            hostname, version, model_name
        )

    def load(self, model_path=None):
        """Override load method, no model to load"""
        pass

    def transform_input(self, *args):
        if len(args) == 0:
            return args[0]
        return args

    def transform_output(self, *args):
        if len(args) == 0:
            return args[0]
        return args

    def predict(self, *args, **kwargs):
        """Makes prediction using TF Serving API.

        Args:
            smiles (str): SMILES string of input molecule
            max_num_templates (int): Maximum number of template scores
                and indices to return
            max_cum_prob (float): Maximum cumulative probability of template
                scores to return. Scores and indices will be returned up until
                max_cum_prob is exceeded.

        Returns:
            (scores, indices): np.ndarrays of scores and indices for 
                prioritized templates
        """
        x = self.transform_input(*args, **kwargs)
        resp = requests.post(self.url, json={'instances': x})
        pred = np.array(resp.json()['predictions']).reshape(-1)
        pred = self.transform_output(pred, **kwargs)
        return pred