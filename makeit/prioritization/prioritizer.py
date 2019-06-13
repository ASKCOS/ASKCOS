class Prioritizer(object):
    """Base class for prioritizers."""

    def __init__(self):
        """Initializes Prioritizer."""
        raise NotImplementedError

    def get_priority(self, object_to_prioritize, **kwargs):
        """Gets the priority for an object that can be prioritized.

        The object is either a retro-synthetic precursor or a tuple of
        a template and a target.
        """
        raise NotImplementedError

    def load_model(self):
        """
        If a neural network model is used to determine the priority, load it!
        """
        raise NotImplementedError

    def set_max_templates(self, max):
        """Sets maximum number of templates to be included in output."""
        self.template_count = max

    def set_max_cum_prob(self, max):
        """Sets maximum cumulative probability of output."""
        self.max_cum_prob = max
