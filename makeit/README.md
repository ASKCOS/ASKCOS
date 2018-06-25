
# MAKE-IT:
Software package for the prediction of feasible synthetic routes towards a desired compound and associated tasks. Is interdependent with ASKCOS_WEBSITE [https://github.com/connorcoley/ASKCOS_Website/].

## Dependencies
The code has primarily been developed for Python 2.7.6 on Ubuntu 16.04. However, we have made an effort to make it work on Python 3.6.1 as well (tested on macOS 10.13.3). There are likely some lingering bugs, so please let us know if you find any.

Included in this repository is a conda environment ```askcos.yml```. Note that this is a fairly messy file and there are some extraneous packages listed. Additionally, some packages are slightly out of date and could be updated without any issues. We are woring to clean this script up and streamline the deployment process. 

## How to install
A coarse installation guide can be found in ```install_cli.sh``` (i.e., "install command line interface"). Note that this shell script is  _not_ meant to actually be run as a shell script. 

We also have an installation guide for the Django web interface, which uses Celery for asynchronous task management (```install_webapp.sh```). This relies on RabbitMQ and Redis servers and uWSGI/NGINX or Apache for deployment.

## Options 
For the integrated synthesis planning tool at ```makeit/application/run.py```, there are several options available. The currently enabled options for the command-line tool can be found at ```makeit/utilities/io/arg_parser.py```. There are some options that are only available for the website and some that are only available for the command-line version. As an example of the former, the consideration of popular but non-buyable chemicals as suitable "leaf nodes" in the search.

### Model choices: the following options influence which models are used to carry out the different
tasks within the algorithm.

- Context recommendation: via '--context_recommender', currently has the following options:
	-'Nearest_Neighbor': Uses a nearest neighbor based database search (memory intensive, ~30GB, and slow; relies on external data file)
	-'Neural_Network': Uses a pretrained neural network (highly recommended!!)

- Context prioritization: via '--context_prioritization', specifies how we should determine the ``best'' context for a proposed reaction. It currently has the following options:
	-'Probability': uses the likelihood of success for the reaction under that condition
	-'Rank': uses the rank of the reaction under that condition relative to all other outcomes

- Forward evaluation: via '--forward_scoring', is used to evaluate the likelihood of success of a reaction. It currently has the following options:
	-'Template_Based': uses the original forward evaluation method enumerating all possible outcomes by applying templates and then predicting the most likely main product [https://pubs.acs.org/doi/abs/10.1021/acscentsci.7b00064]
    -'Template_Free': uses the higher-performing and faster template-free method based on graph convolutional neural networks [https://arxiv.org/abs/1709.04555]
    -'Fast_Filter': uses a binary classifier to distinguish good and bad reaction suggestions. It is imperfect, but very fast. Based on the ``in-scope filter'' suggested by Marwin Segler [https://www.nature.com/articles/nature25978]

- Retrosynthetic template prioritization: via '--template_prioritization', is used to minimize the number of reaction templates that must be applied to the target compound at each iteration. It currently has the following options:
	-'Relevance': Quantifies how relevant a given template is for the considered reactants, based on the approach suggested by Marwin Segler [https://onlinelibrary.wiley.com/doi/abs/10.1002/chem.201605499]
	-'Popularity': Ranking based on number of references in literature, independent of the product species

- Precursor prioritization: via '--precusor_prioritization', is used to determine which precursor is the most promising branch to pursue. It currently has the following options:
	-'Heuristic': Simple heuristic, with decreasing score as number of atoms, rings and chiral centers increases
	-'SCScore': Synthetic Complexity Score - learned quantity indicating how complex a molecule is. Tries to interpret molecules with a protection/deprotection group as less complex than their non-protected counterparts. [https://pubs.acs.org/doi/abs/10.1021/acs.jcim.7b00622]



- Tree scoring: via '--tree_scoring', determines how final synthesis trees should be sorted/ranked. It currently has the following options:
	-'Product': uses the product of template score and forward prediction score
	-'Forward_only': uses only the forward prediction score
	-'Template_only': uses only the template score

### Limits and thresholds: the following options will set limits for the different parts of the program

- Expansion time via '--expansion_time': limit the amount of time the program spends expanding the retro synthetic tree. Default value is 60 seconds.
  
- Maximum search depth via '--max_depth': limit the search depth in the retro synthetic expansion. Default value is 4.
  
- Maximum degree of branching via '--max_branching': limit the number of branches generated in each layer of the retro synthetic tree. Default value is 20
  
- Maximum number of buyable trees via '--max_trees': limit the number of buyable trees the program  should search for. Default value is 500.

- Maximum number of templates to be applied via '--template_count': limit the number of templates that are appied for each expansion when using the popularity prioritizer. Default value is 10000.

- Minimal number of templates to be considered in retro synthetic direction for non-chiral reactions via '--mincount_retro'. Default value is 25.
  
- Minimal number of templates to be considered in retro synthetic direction for chiral reactions via '--mincoun_retro_c'. Default value is 10.
 
- Minimal number of templates to be considered in synthetic direction via '--synth_mincount'. Default value is 25.
  
- Minimal target rank for considering a target feasible via '--rank_threshold'. Default value is 10.

- Minimal probability for considering a target feasible via '--prob_threshold'. Default value is 0.01.

- Maximum number of contexts to be proposed for each reaction via '--max_contexts'. Default value is 10

- Maximum price per gram for a component to be considered buyable via '--max_ppg'. Default value is 100
  
- Precursor filtering: via '--apply_fast_filter' and '--filter_threshold', is used to impose rapid filtering of low-quality retrosynthetic suggestions. Default is True and 0.75
