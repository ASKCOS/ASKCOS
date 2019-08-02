from makeit.prioritization.prioritizer import Prioritizer
from makeit.utilities.io.logger import MyLogger
popularity_template_prioritizer_loc = 'popularity_template_prioritizer'


class PopularityTemplatePrioritizer(Prioritizer):
    """A template Prioritizer ordering by popularity.

    Allows to prioritize a template based on the number of times it appears,
    or its reported popularity.

    Attributes:
        sorted (bool): Whether the templates have been sorted yet.
        reordered_templates (list of ??): Templates ordered by popularity.
        template_count (int): Number of templates to return the priority of.
        max_cum_prob (float): Maximum cumulative probability of returned
            tempates. Unused.
    """

    def __init__(self, no_log=False):
        """Initializes PopularityTemplatePrioritizer.

        Args:
            no_log (bool, optional): Whether to not log. Unused.
                (default: {False})
        """
        # only do 'reorder' once
        self.sorted = False
        self.reordered_templates = None
        self.template_count = 1e9
        self.max_cum_prob = 1

    def get_priority(self, input_tuple, **kwargs):
        """Returns list of templates ordered by popularity.

        Args:
            input_tuple (2-tuple of (list of ??, ??)): Templates to get the
                priority of.
            **kwargs: Additional optional parameters. Used for template_count.
        """
        (templates, target) = input_tuple
        return self.reorder(templates)[:min(len(templates), kwargs.get('template_count', 1e10))]

    def reorder(self, templates):
        """Reorders templates by popularity.

        Re-orders the list of templates according to field 'count' in descending
        order. This means we will apply the most popular templates first.

        Args:
            templates (list of ??): Unordered templates to be reordered by
                popularity.

        Returns:
            list of ??: Templates sorted by popularity.
        """
        if self.sorted:
            return self.reordered_templates
        else:
            templates[:] = [x for x in sorted(
                templates, key=lambda z: z['count'], reverse=True)]
            self.sorted = True
            for template in templates:
                template['score'] = 1
            self.reordered_templates = templates
            return templates

    def load_model(self):
        """Loads popularity model.

        PopularityTemplatePrioritizer does not use a neural network, so this
        does nothing.
        """
        pass
