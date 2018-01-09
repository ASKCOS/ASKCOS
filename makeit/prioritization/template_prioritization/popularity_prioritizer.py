from prioritization.prioritizer import Prioritizer
from utilities.i_o.logging import MyLogger
popularityPrioritizer_loc = 'popularity_prioritizer'

class PopularityPrioritizer(Prioritizer):
    '''
    Allows to prioritize a template based on the number of times it appears, or its reported popularity
    '''
    def __init__(self, no_log = False):
        if not no_log:
            MyLogger.print_and_log('Prioritizing templates based on popularity in the literature.', popularityPrioritizer_loc)
        #only do 'reorder' once
        self.sorted = False
        self.reordered_templates = None
        
    def get_priority(self, input_tuple):
        (templates, target) = input_tuple
        return self.reorder(templates)
    
    def reorder(self, templates):
        '''
        Re-orders the list of templates (self.templates) according to 
        field 'count' in descending order. This means we will apply the
        most popular templates first
        '''
        if self.sorted:
            return self.reordered_templates
        else:
            templates[:] = [x for x in sorted(templates, key = lambda z: z['count'], reverse = True)]
            self.sorted = True
            for template in templates:
                template['score'] = 1
            self.reordered_templates = templates
            return templates
    
    def load_model(self):
        pass