

class ForwardEnumerator(object):
    """Interface for forward enumeration classes.

    At least an initialization method and a ``get_outcomes`` of the enumeration
    should be present.
    """
    def __init__(self):
        raise NotImplementedError

    def get_outcomes(self, smiles):
        """Enumerates the possible products for a given smiles."""
        raise NotImplementedError
