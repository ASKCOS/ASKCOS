class Prioritizer(object):
    
    def __init__(self):
        raise NotImplementedError
    
    def get_priority(self, object_to_prioritize, **kwargs):
        '''
        Get the priority for an object that can be prioritized. This is either a retro-synthetic precursor or a tuple of
        a template and a target.
        '''
        raise NotImplementedError
    
    def load_model(self):
        '''
        If a neural network model is used to determine the priority, load it!
        '''
        raise NotImplementedError
    
    def set_max_templates(self, max):
        self.template_count = max
    
