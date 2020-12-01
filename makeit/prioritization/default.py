from makeit.prioritization.prioritizer import Prioritizer
from makeit.utilities.io.logger import MyLogger
default_prioritizer_loc = 'default_prioritizer'


class DefaultPrioritizer(Prioritizer):
    """A default Prioritizer that assigns the same priority to everything.

    Will not prioritize the objects in any way. Returns a priority of 1 for any
    object. Can therefore be used both for (non)-prioritization of templates
    and/or retro-synthetic precursors.
    """

    def __init__(self):
        """Initializes DefaultPrioritizer."""
        pass
    def get_priority(self, object_to_prioritize, **kwargs):
        """Returns priority of given object.

        If object is a tuple (for prioritization of a template), then return the
        first element (the list of templates).
        Otherwise returns 1.

        Args:
            object_to_prioritize (2-tuple or ??): Object to prioritize.
            **kwargs: Unused.

        Returns:
            list or float: Unsorted input list of templates or 1.0
        """
        try:
            (templates, target) = object_to_prioritize
            return templates
        # if not a tuple: prioritization of a retro-precursor element.
        except TypeError:
            return 1.0

    def load_model(self):
        """Loads default model.

        DefaultPrioritizer does not use a neural network, so this does nothing.
        """
        pass
