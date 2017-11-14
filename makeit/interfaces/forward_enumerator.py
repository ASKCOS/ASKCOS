

class ForwardEnumerator(object):
    '''
    Interface for forward enumeration classes. At least an initialization method and a 'get outcomes' of the enumeration 
    should be present.
    '''
    def __init__(self):
        raise NotImplementedError
    
    def get_outcomes(self, smiles):
        '''
        Use the enumerator's methods to get the possible products for a certain smiles
        '''
        raise NotImplementedError
    