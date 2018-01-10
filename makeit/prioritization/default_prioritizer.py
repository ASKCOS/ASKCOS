from makeit.prioritization.prioritizer import Prioritizer
from makeit.utilities.i_o.logging import MyLogger
defaultPrioritizer_loc = 'default_rioritizer'

class DefaultPrioritizer(Prioritizer):
    '''
    Will not prioritize the objects in any way. Returns a priority of 1 for any object. Can therefore be used both for (non)-
    prioritization of templates and/or retro-synthetic precursors.
    '''
    
    def __init__(self):
        MyLogger.print_and_log('No specific prioritization will be carried out.', defaultPrioritizer_loc)
        
    def get_priority(self, object_to_prioritize):
        try:
            (templates, target) = object_to_prioritize
            return templates
        #if not a tuple: prioritization of a retro-precursor element.
        except TypeError:
            return 1.0
    
    def load_model(self):
        pass