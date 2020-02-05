# ASKCOS:
Software package for the prediction of feasible synthetic routes towards a desired compound and associated tasks related to synthesis planning. Originally developed under the DARPA Make-It program and now being developed under the [MLPDS Consortium](http://mlpds.mit.edu).

# v0.4.1 Release

Release Notes  
User notes:  
* New impurity predictor module
* New reaction atom-mapping module
* Upgrade to rdkit version 2019.03.3
* Migration of rdchiral to standalone pypi package. rdchiral development can now be found at https://github.com/connorcoley/rdchiral
* Improved buyables lookup consistency
* Canonicalizes SMILES strings before lookup in buyables module
* Improved granularity of feedback after buyable upload
* Improved handling of reaction templates in rdchiral for "hypervalent" nitrogens. Significant improvement for nitration reactions

Developer notes:
* makeit data migrated to a separate repository (https://gitlab.com/mlpds_mit/ASKCOS/makeit-data)
* Data copied from data-only docker container (registry.gitlab.com/mlpds_mit/askcos/makeit-data:0.4.1)
* Docker builds are now much, much, much faster
* chemhistorian data migrated to mongodb. This increases initialization mongodb seeding time, but decreases memory footprint
* All dependencies, including third-party docker images, are now pinned to specific versions
* Seeding of mongo db now occurs using the backend "app" service
* Template relevance and fast filter models moved to tensorflow serving API endpoints

Bug fixes:
* Buyables page bugfixes
* nginx service restarts like the rest of the services now

### Using GitLab Deploy Tokens

ASKCOS can also be downloaded using deploy tokens, these provide __read-only__ access to the source code and our container registry in GitLab. Below is a complete example showing how to deploy the ASKCOS application using deploy tokens (omitted in this example). The deploy tokens can be found on the [MLPDS Member Resources ASKCOS Versions Page](https://mlpds.mit.edu/member-resources-releases-versions/). The only software prerequisites are git, docker, and docker-compose.

```bash
$ export DEPLOY_TOKEN_USERNAME=
$ export DEPLOY_TOKEN_PASSWORD=
$ git clone https://$DEPLOY_TOKEN_USERNAME:$DEPLOY_TOKEN_PASSWORD@gitlab.com/mlpds_mit/askcos/askcos.git
$ docker login registry.gitlab.com -u $DEPLOY_TOKEN_USERNAME -p $DEPLOY_TOKEN_PASSWORD
$ cd askcos/deploy
$ git checkout v0.4.1
$ bash deploy.sh deploy
```

__NOTE:__ Starting with version 0.4.1, the chemhistorian data has been migrated to mongodb, which may take up to ~5 minutes to initially seed for the first time upgrade/deployment. Subsequent upgrades should not require the re-seeding of the chemhistorian information.

### Upgrade information

The easiest way to upgrade to a new version of ASKCOS is using Docker and docker-compose.
To get started, make sure both docker and docker-compose are installed on your machine.
We have a pre-built docker image of ASKCOS hosted on GitLab.
It is a private repository; if you do not have access to pull the image, please [contact us](mailto:mlpds_support@mit.edu).
In addition, you need to have the deploy/ folder from the ASKCOS code repository for the specific version for which you would like to upgrade to. Due to backend changes introduced with v0.3.1, the upgrade information is different is older versions, and the steps are summarized below:

#### From v0.3.1 or v0.4.0
```
$ git checkout v0.4.1
$ bash deploy.sh start                      # updates services to v0.4.1
$ bash deploy.sh set-db-defaults seed-db    # this may take ~5 minutes to load "default chemicals data" (new in 0.4.1)
```

#### From v0.2.x or v0.3.0
```
$ git checkout v0.4.1
$ bash backup.sh
$ bash deploy.sh start                      # updates services to v0.4.1
$ bash deploy.sh set-db-defaults seed-db    # this may take ~5 minutes to load "default chemicals data" (new in 0.4.1)
$ bash restore.sh
```

# First Time Deployment with Docker

### Prerequisites

 - If you're buidling the image from scratch, make sure git (and git lfs) is installed on your machine
 - Install Docker [OS specific instructions](https://docs.docker.com/install/)
 - Install docker-compose [installation instructions](https://docs.docker.com/compose/install/#install-compose)

### Deploying the web application

Deployment is initiated by a bash script that runs a few docker-compose commands in a specific order.
Several database services need to be started first, and more importantly seeded with data, before other services 
(which rely on the availability of data in the database) can start. The bash script can be found and should be run 
from the deploy folder as follows:

```bash
$ bash deploy.sh command [optional arguments]
```

There are a number of available commands, including the following for common deploy tasks:
* `deploy`: runs standard first-time deployment tasks, including `seed-db`
* `update`: pulls new docker image from GitLab repository and restarts all services
* `seed-db`: seed the database with default or custom data files
* `start`: start a deployment without performing first-time tasks
* `stop`: stop a running deployment
* `clean`: stop a running deployment and remove all docker containers

For a running deployment, new data can be seeded into the database using the `seed-db` command along with arguments
indicating the types of data to be seeded. Note that this will replace the existing data in the database.
The available arguments are as follows:
* `-b, --buyables`: specify buyables data to seed, either `default` or path to data file
* `-c, --chemicals`: specify chemicals data to seed, either `default` or path to data file
* `-x, --reactions`: specify reactions data to seed, either `default` or path to data file
* `-r, --retro-templates`: specify retrosynthetic templates to seed, either `default` or path to data file
* `-f, --forward-templates`: specify forward templates to seed, either `default` or path to data file

For example, to seed default buyables data and custom retrosynthetic pathways, run the following from the deploy folder:

```bash
$ bash deploy.sh seed-db --buyables default --retro-templates /path/to/my.retro.templates.json.gz
```

To update a deployment, run the following from the deploy folder:

```bash
$ bash deploy.sh update --version x.y.z
```

To stop a currently running application, run the following from the deploy folder:

```bash
$ bash deploy.sh stop
```

If you would like to clean up and remove everything from a previous deployment (__NOTE: you will lose user data__), run the following from the deploy folder:

```bash
$ bash deploy.sh clean
```

### Backing up user data

If you are upgrading from v0.3.1 or later, the backup/restore process is no longer needed unless you are moving deployments to a new machine.

If you are upgrading the deployment from a previous version (prior to v0.3.1), or moving the application to a different server, you may want to retain user accounts and user-saved data/results.
The provided `backup.sh` and `restore.sh` scripts are capable of handling the backup and restoring process. Please read the following carefully so as to not lose any user data:

1) Start by making sure the previous version you would like to backup is __currently up and running__ with `docker-compose ps`.
2) Checkout the newest version of the source code `git checkout v0.4.1`
3) Run `$ bash backup.sh`
4) Make sure that the `deploy/backup` folder is present, and there is a folder with a long string of numbers (year+month+date+time) that corresponds to the time you just ran the backup command
5) If the backup was successful (`db.json` and `user_saves` (\<v0.3.1) or `results.mongo` (\>=0.3.1) should be present), you can safely tear down the old application with `docker-compose down [-v]`
6) Deploy the new application with `bash deploy.sh deploy` or update with `bash deploy.sh update -v x.y.z`
7) Restore user data with `bash restore.sh`

Note: For versions >=0.3.1, user data persists in docker volumes and is not tied to the lifecycle of the container services. If the [-v] flag is not used with `docker-compose down`, volumes do not get removed, and user data is safe. In this case, the backup/restore procedure is not necessary as the containers that get created upon an install/upgrade will continue to use the docker volumes that contain all the important data. If the [-v] flag is used, all data will be removed and a restore will be required to recover user data.

### (Optional) Building the ASKCOS Image

The askcos image itself can be built using the Dockerfile in this repository.

```bash
$ git clone https://gitlab.com/mlpds_mit/askcos/askcos  
$ cd askcos/  
$ git lfs pull   
$ docker build -t askcos .
```

__NOTE:__ For application deployment, double check the image tag used in the `docker-compose.yml` file and be sure to tag your newly built image with the same image name. Otherwise, the image tag sed in `docker-compose.yml` will be pullled and deployed instead of the image that was just built.

### Add customization

There are a few parts of the application that you can customize:
* Header sub-title next to ASKCOS (to designate this as a local deployment at your organization)

This is handled as an environment variable that can change upon deployment (and are therefore not tied into the image directly). This can be found in `deploy/customization`. Please let us know what other degrees of customization you would like.

### Managing Django

If you'd like to manage the Django app (i.e. - run python manage.py ...), for example, to create an admin superuser, you can run commands in the _running_ app service (do this _after_ `docker-compose up`) as follows:

`docker-compose exec app bash -c "python /usr/local/ASKCOS/askcos/manage.py createsuperuser"`

In this case you'll be presented an interactive prompt to create a superuser with your desired credentials.

## Important Notes

#### Scaling workers

Only 1 worker per queue is deployed by default with limited concurrency. This is not ideal for many-user demand. You can easily scale the number of celery workers you'd like to use with `docker-compose up -d --scale tb_c_worker=N` where N is the number of workers you want, for example. The above note applies to each worker you start, however, and each worker will consume RAM.


# How to run individual modules
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

#### Integrated CASP tool
For the integrated synthesis planning tool at ```makeit/application/run.py```, there are several options available. The currently enabled options for the command-line tool can be found at ```makeit/utilities/io/arg_parser.py```. There are some options that are only available for the website and some that are only available for the command-line version. As an example of the former, the consideration of popular but non-buyable chemicals as suitable "leaf nodes" in the search. An example of how to use this module is:

```python ASKCOS/Make-It/makeit/application/run.py --TARGET atropine```

##### Model choices.
The following options influence which models are used to carry out the different tasks within the algorithm.

- Context recommendation: via '--context_recommender', currently has the following options:

	-'Nearest_Neighbor': Uses a nearest neighbor based database search (memory intensive, ~30GB, and slow; relies on external data file)
	
	-'Neural_Network': Uses a pretrained neural network (highly recommended!!)

- Context prioritization: via '--context_prioritization', specifies how we should determine the "best" context for a proposed reaction. It currently has the following options:

	-'Probability': uses the likelihood of success for the reaction under that condition
	
	-'Rank': uses the rank of the reaction under that condition relative to all other outcomes

- Forward evaluation: via '--forward_scoring', is used to evaluate the likelihood of success of a reaction. It currently has the following options:

	-'Template_Based': uses the original forward evaluation method enumerating all possible outcomes by applying templates and then predicting the most likely main product [https://pubs.acs.org/doi/abs/10.1021/acscentsci.7b00064] (NOTE: the template-based forward predictor requires a custom built version of RDKit from https://github.com/connorcoley/rdkit - we highly recommend using the template-free approach)
	
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

#### Limits and thresholds: the following options will set limits for the different parts of the program

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
