

class Scorer():
    """Interface for scorer classes."""
    def __init__(self, args, kwargs={}):
        raise NotImplementedError

    def evaluate(self, args, kwargs={}):
        """Scores reactions with different contexts.

        Implemented method should return:
        
        * A list of results (one for each context).

          Each result should contain:

          * A list of outcomes, of which each outcome is a dictionnary
            containing:

            * rank
            * forward result
            * score
            * probability
        """
        raise NotImplementedError
    def stop_expansion(self):
        """Method to kill all spun up workers."""
        raise NotImplementedError
