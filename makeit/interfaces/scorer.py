

class Scorer():
    
    def __init__(self, args, kwargs={}):
        raise NotImplementedError
    
    def evaluate(self, args, kwargs={}):
        '''
        Implemented method should return:
        A list of results (one for each context)
        Each result should contain:
        A list of outcomes, of which each outcome is a dictionnary containing:
        a rank, a forward result, a score and a probability.
        '''
        raise NotImplementedError
    def stop_expansion(self):
        '''
        Method to kill all spun up workers
        '''
        raise NotImplementedError