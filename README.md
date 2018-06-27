
# Make-It:
Software package for the prediction of feasible synthetic routes towards a desired compound and associated tasks. Is interdependent with ASKCOS_WEBSITE [https://github.com/connorcoley/ASKCOS_Website/].

## Dependencies
The code has primarily been developed for Python 2.7.6 on Ubuntu 16.04. However, we have made an effort to make it work on Python 3.6.1 as well (tested on macOS 10.13.3). There are likely some lingering bugs, so please let us know if you find any.

Included in this repository is a conda environment ```askcos.yml```. Note that this is a fairly messy file and there are some extraneous packages listed. Additionally, some packages are slightly out of date and could be updated without any issues. We are woring to clean this script up and streamline the deployment process. 

## How to install
A coarse installation guide can be found in ```install_cli.sh``` (i.e., "install command line interface"). Note that this shell script is  _not_ meant to actually be run as a shell script. 

We also have an installation guide for the Django web interface, which uses Celery for asynchronous task management (```install_webapp.sh```). This relies on RabbitMQ and Redis servers and uWSGI/NGINX or Apache for deployment.

Please note that this code relies on either (1) additional data files not contained in this repo, but available from ccoley@mit.edu or (2) connection to a MongoDB with specific expectations for databases/collections/etc.

#### Make-It installation instructions (Ubuntu 16.04, Python 2.7.6)
1. Place the Make-It and ASKCOS_Website repositories inside an ```ASKCOS``` folder in your home folder

1. Download miniconda2 (or miniconda3 for Python 3) and add it to the system path
	```
	wget https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh
	bash Miniconda2-latest-Linux-x86_64.sh -b -p $HOME/miniconda
	export PATH=~/miniconda/bin:$PATH
	echo 'export PATH=~/miniconda/bin:$PATH' >> ~/.bashrc
	```

1. Prepare system-wide dependencies before setting up python packages

	```
	sudo apt install heimdal-dev
	sudo apt install libkrb5-dev 
	sudo apt-get install build-essential
	```

1. Create the "askcos" conda environment using the provided ```askcos.yml``` as a starting point

	```
	conda-env create -f askcos.yml -n askcos
	```

1. Add Make-It and ASKCOS_Website folders to your ```PYTHONPATH``` environment variable

	```
	export PYTHONPATH=~/ASKCOS/Make-It:~/ASKCOS/ASKCOS_Website:$PYTHONPATH
	echo 'export PYTHONPATH=~/ASKCOS/Make-It:~/ASKCOS/ASKCOS_Website:$PYTHONPATH' >> ~/.bashrc 
	```

1. [OPTIONAL] create a link between ```Make-It/makeit/data``` and wherever you actually want to store the data

1. Install RDKit and Pymongo

	```
	source activate askcos
	conda install rdkit 
	pip install pymongo
	sudo apt install libxrender-dev 

	export PYTHONPATH=~/miniconda/envs/askcos/lib/python2.7/site-packages:$PYTHONPATH
	export LD_LIBRARY_PATH=~/miniconda/envs/askcos/lib:$LD_LIBRARY_PATH

	echo 'export PYTHONPATH=~/miniconda/envs/askcos/lib/python2.7/site-packages:$PYTHONPATH' >> ~/.bashrc
	echo 'export LD_LIBRARY_PATH=~/miniconda/envs/askcos/lib:$LD_LIBRARY_PATH' >> ~/.bashrc
	```

1. Install a more recent verison of Tensorflow

	```
	pip install tensorflow==1.4.1
	```

1. Set the default Keras backend to Theano by editing ```~/.keras/keras.json```

1. Test Make-It by running any of the individual modules below and/or the full planning script.

#### ASKCOS_Website installation instructions (Ubuntu 16.04, Python 2.7.6)

1. Install Redis as the results backend

	```
	sudo add-apt-repository -y ppa:chris-lea/redis-server 
	sudo apt-get update
	sudo apt-get install redis-server
	sudo service redis start 
	```

1. Install RabbitMQ as the message broker

	```
	echo "deb https://dl.bintray.com/rabbitmq/debian xenial main" | sudo tee /etc/apt/sources.list.d/bintray.rabbitmq.list 
	wget -O- https://dl.bintray.com/rabbitmq/Keys/rabbitmq-release-signing-key.asc |
	     sudo apt-key add - 
	sudo apt-get update 

	wget https://packages.erlang-solutions.com/erlang-solutions_1.0_all.deb 
	sudo dpkg -i erlang-solutions_1.0_all.deb
	sudo apt-get update
	sudo apt-get install erlang 
	sudo apt-get install rabbitmq-server # on OS X, brew install rabbitmq

	sudo service redis start 
	sudo service rabbitmq-server start
	```
	
1. Install uWSGI

	```
	pip install uwsgi
	sudo apt-get install nginx
	sudo /etc/init.d/nginx start 
	```
	
1. Test that uWSGI is working by temporarily running the webserver

	```
	uwsgi --http :8000 --wsgi-file wsgi.py
	```

	and in a new terminal window...

	```
	curl http://localhost:8000
	```

1. Inside the ASKCOS_Website folder, get uwsgi_params

	```
	wget https://raw.githubusercontent.com/nginx/nginx/master/conf/uwsgi_params
	```

1. Edit /etc/nginx/nginx.conf as needed, according to http://uwsgi-docs.readthedocs.io/en/latest/tutorials/Django_and_nginx.html

1. Start the uWSGI server (long-term, this should be run as a service)

	```
	uwsgi --socket :8000 --wsgi-file wsgi.py
	```

1. Start the NGINX server

	```
	sudo service nginx start
	```
	
1. Test basic website functionality
	- SCScore
	- Buyable page
	- Drawing
	- Saving/restoring pages
	
1. Test Celery functionality by starting up context recommender workers
	1. 

		```
		celery -A askcos_site worker -c 1 -Q cr_coordinator -n "cr_coordinator@$(hostname)" --max-tasks-per-child 1000 --loglevel=$CELERY_LOG_LEVEL  --logfile=celery_logs/%p.log &
		celery -A askcos_site worker -c 1 -Q cr_network_worker -n "cr_network_worker@$(hostname)" --max-tasks-per-child 5000 --loglevel=$CELERY_LOG_LEVEL  --logfile=celery_logs/%p.log &
		```
	1. Go to the /context/ page on the website and make a request for the neurla network context recommender model
	
1. Edit ```spawn_workers.sh``` to reflect the anticipated server load (i.e., number of workers needed to support each task) before running the script to spawn them.

	_note: if you want the workers to run on a separate server, you will need to edit ```ASKCOS_Website/askcos_site/celery.py``` to give an explicit SERVERHOST instead of localhost. You will also need to open up the ports for the message broker (5672) and redis server (6379) to the other server._


#### Setting up the Celery workers on a different server from the Webserver

For scalability, it is possible to have the webserver and compute servers completely separately. The webserver will host the website, RabbitMQ, and Redis servers. The compute server will be where Celery workers are spun up.

1. If using AWS, open up connections between instances in the same security group (following the slightly-outdated https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/security-group-rules-reference.html#sg-rules-other-instances)

1. [On the webserver] Allow remote connections to RabbitMQ. It is possible to create specific users for the default vhost to limit access, but it is easiest to open up all remote connections to the "guest" user (https://www.rabbitmq.com/access-control.html). This relies on non-RabbitMQ security measures to prevent unwanted connections to the server's port. Create ```/etc/rabbitmq/rabbitmq.conf``` if it does not exits and enter one line:

	```
	loopback_users = none
	```

1. [On the webserver] Edit ```/etc/redis/redis.conf``` to allow remote connections from worker machine. Ideally, you can bind to a specific IP by editing the ```bind 127.0.0.1``` line. This seems to work fine on our local servers but did not work well on AWS. The alternative is to completely open up connections by changing the ```bind``` line to ```bind 0.0.0.0``` and edit the following block to say ```protected-mode off```. As with the RabbitMQ settings, this means that Redis will rely on other firewalls/security measures.


1. [On the compute server] Change the ```SERVERHOST``` in ```ASKCOS_Website/askcos_site/celery.py``` to the IP of the webserver, so the Celery workers know where to look for the RabbitMQ/Redis servers. If using AWS, this should be the _private_ IP of the webserver instance.

1. [On the compute server] Spin up as many workers as desired using ```spawn_workers.sh``` from inside the ```ASKCOS_Website``` folder. The most computationally intensive task is the TreeBuilder, for which we generally allocate 4-12 processes in each worker pool. Each process requires 3-4 GB of RAM. The other Celery workers (aside from the nearest neighbor context recommender, which should be deprecated) are relatively small in memory footprint.


## How to run individual modules
Many of the individual modules -- at least the ones that are the most interesting -- can be run "standalone". Examples of how to use them are often found in the ```if __name__ == '__main__'``` statement at the bottom of the script definitions. For example...

#### Using the learned synthetic complexity metric (SCScore)
```makeit/prioritization/precursors/scscore.py```

#### Obtaining a single-step retrosynthetic suggestion with consideration of chirality
```makeit/retrosynthetic/transformer.py```

#### Finding recommended reaction conditions based on a trained neural network model
```makeit/synthetic/context/neuralnetwork.py```

#### Using the template-free forward predictor
```makeit/synthetic/evaluation/template_free.py```

#### Using the coarse "fast filter" (binary classifier) for evaluating reaction plausibility
```makeit/synthetic/evaluation/fast_filter.py```

## Integrated CASP tool
For the integrated synthesis planning tool at ```makeit/application/run.py```, there are several options available. The currently enabled options for the command-line tool can be found at ```makeit/utilities/io/arg_parser.py```. There are some options that are only available for the website and some that are only available for the command-line version. As an example of the former, the consideration of popular but non-buyable chemicals as suitable "leaf nodes" in the search. An example of how to use this module is:

```python ASKCOS/Make-It/makeit/application/run.py --TARGET atropine```

### Model choices: the following options influence which models are used to carry out the different
tasks within the algorithm.

- Context recommendation: via '--context_recommender', currently has the following options:

	-'Nearest_Neighbor': Uses a nearest neighbor based database search (memory intensive, ~30GB, and slow; relies on external data file)
	
	-'Neural_Network': Uses a pretrained neural network (highly recommended!!)

- Context prioritization: via '--context_prioritization', specifies how we should determine the "best" context for a proposed reaction. It currently has the following options:

	-'Probability': uses the likelihood of success for the reaction under that condition
	
	-'Rank': uses the rank of the reaction under that condition relative to all other outcomes

- Forward evaluation: via '--forward_scoring', is used to evaluate the likelihood of success of a reaction. It currently has the following options:

	-'Template_Based': uses the original forward evaluation method enumerating all possible outcomes by applying templates and then predicting the most likely main product [https://pubs.acs.org/doi/abs/10.1021/acscentsci.7b00064]
	
    -'Template_Free': uses the higher-performing and faster template-free method based on graph convolutional neural networks [https://arxiv.org/abs/1709.04555]
    
    -'Fast_Filter': uses a binary classifier to distinguish good and bad reaction suggestions. It is imperfect, but very fast. Based on the "in-scope filter" suggested by Marwin Segler [https://www.nature.com/articles/nature25978]

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
