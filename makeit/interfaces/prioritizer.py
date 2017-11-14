
class Prioritizer(object):
    
    def __init__(self):
        raise NotImplementedError
    
    def give_priority_score(self, retroPrecursor):
        '''
        Priority score should always be calculated in some way based on the RetroPrecursor object
        '''
        raise NotImplementedError