import numpy as np
import requests

class TFServingAPIModel(object):
    """Base tensorflow serving API Model class.
    
    Attributes:
        hostname (str): hostname of service serving tf model.
        model_name (str): Name of model provided to tf serving.
        version (int): version of the model to use when serving

    """
    def __init__(self, hostname, model_name, version=1):
        self.url = 'http://{}:8501/v{}/models/{}:predict'.format(
            hostname, version, model_name
        )

    def load_model(self, model_path=None):
        """Override load method, no model to load"""
        pass

    def transform_input(self, *args):
        """Identity transformation. Return the arguments as they were passed. If a single argument was passed, return it as a single argument and not a list."""
        if len(args) == 1:
            return args[0]
        return args

    def transform_output(self, *args):
        """Identity transformation. Return the arguments as they were passed. If a single argument was passed, return it as a single argument and not a list."""
        if len(args) == 1:
            return args[0]
        return args

    def predict(self, *args, **kwargs):
        """Makes a prediction using TF Serving API.
        Calls a transformation function before and after actually calling the API endpoint to allow for customizable pipelines.
        """
        x = self.transform_input(*args, **kwargs)
        resp = requests.post(self.url, json={'instances': x})
        pred = np.array(resp.json()['predictions']).reshape(-1)
        pred = self.transform_output(pred, **kwargs)
        return pred