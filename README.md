# ASKCOS:
Software package for the prediction of feasible synthetic routes towards a desired compound and associated tasks related to synthesis planning. Originally developed under the DARPA Make-It program and now being developed under the [MLPDS Consortium](http://mlpds.mit.edu).

Please note that the MPL 2.0 license for this repository does not apply to the data and trained models. The data and trained models are released under CC BY-NC-SA (i.e., are for noncommercial use only).

Contributors include Connor Coley, Mike Fortunato, Hanyu Gao, Pieter Plehiers, Matthew Cameron, Max Liu, Yuran Wang, Thomas Struble, Jiannan Liu, and Yiming Mo.

# Where to find newer ASKCOS releases

This repository contains the original ASKCOS project and will no longer be updated, with the final release being v0.4.1. However, the ASKCOS project is still under very active development! The project has been split into multiple smaller repositories for many reasons including improved modularity. The askcos, deploy, and makeit directories have been migrated to [askcos-site](https://github.com/ASKCOS/askcos-site), [askcos-deploy](https://github.com/ASKCOS/askcos-deploy), and [askcos-core](https://github.com/ASKCOS/askcos-core), respectively. In addition, version numbers have changed to a date based system, with the first new open-source release being version [2020.07](https://askcos.github.io/#/2020-07/). We highly encourage you to check out and update to the new version! For any questions or comments, please contact askcos_support@mit.edu.

# Quick start using Google Cloud

```
# (1) Create a Google Cloud instance
#     - recommended specs: 8 vCPUs, 64 GB memory
#     - select Ubuntu 18.04 LTS Minimal
#     - upgrade to a 100 GB disk
#     - allow HTTP and HTTPS traffic

# (2) Install docker
#     - https://docs.docker.com/engine/install/ubuntu/
sudo apt-get update
sudo apt-get install \
    apt-transport-https \
    ca-certificates \
    curl \
    gnupg-agent \
    software-properties-common -y
curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo apt-key add -
sudo add-apt-repository \
   "deb [arch=amd64] https://download.docker.com/linux/ubuntu \
   $(lsb_release -cs) \
   stable"
sudo apt-get update
sudo apt-get install docker-ce docker-ce-cli containerd.io -y
sudo groupadd docker
sudo usermod -aG docker $USER
newgrp docker

# (3) Install docker-compose
#     - https://docs.docker.com/compose/install/
sudo curl -L "https://github.com/docker/compose/releases/download/1.27.4/docker-compose-$(uname -s)-$(uname -m)" -o /usr/local/bin/docker-compose
sudo chmod +x /usr/local/bin/docker-compose
sudo ln -s /usr/local/bin/docker-compose /usr/bin/docker-compose

# (4) Install git lfs
#     - https://github.com/git-lfs/git-lfs/wiki/Installation
sudo apt-get install software-properties-common -y
sudo add-apt-repository ppa:git-core/ppa -y
curl -s https://packagecloud.io/install/repositories/github/git-lfs/script.deb.sh | sudo bash
sudo apt-get install git-lfs -y
git lfs install

# (5) Pull
git clone https://github.com/ASKCOS/ASKCOS
cd ASKCOS
git lfs pull

# (6) Build & run
docker build -t askcos/askcos . # build
cd deploy
bash deploy.sh deploy     # start containers (detached) and run other initialization tasks
docker-compose logs -f    # start tailing logs (can CTRL+C to exit)

# (7) Navigate to your instance's external IP
#     - note that it may take ~5 minutes for the retro transformer workers to start up
#     - you can check the status of their startup by looking at "server status"
#     - the first request to a website process may take ~10 seconds
#     - the first request to a retro transform worker may take ~5-10 seconds
#     - the first request to the forward predictor may take ~60 seconds
```

# First Time Deployment with Docker

## Prerequisites

 - If you're building the image from scratch, make sure git (and git lfs) is installed on your machine
 - Install Docker [OS specific instructions](https://docs.docker.com/install/)
 - Install docker-compose [installation instructions](https://docs.docker.com/compose/install/#install-compose)

## (Optional) Building the ASKCOS Image

The askcos image itself can be built using the Dockerfile in this repository.

```bash
$ git clone https://github.com/ASKCOS/ASKCOS
$ cd ASKCOS
$ git lfs pull
$ docker build -t askcos/askcos .
```

__NOTE:__ For application deployment, double check the image tag used in the `docker-compose.yml` file and be sure to tag your newly built image with the same image name. Otherwise, the image tag used in `docker-compose.yml` will be pulled and deployed instead of the image that was just built.

## Deploying the Web Application

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

# Important Notes

## Recommended hardware

We recommend running this code on a machine with at least 8 compute cores (16 preferred) and 64 GB RAM (128 GB preferred).

## First startup

The celery worker will take a few minutes to start up (possibly up to 5 minutes; it reads a lot of data into memory from disk). The web app itself will be ready before this, however upon the first get request (only the first for each process) a few files will be read from disk, so expect a 10-15 second delay.

## Scaling workers

Only 1 worker per queue is deployed by default with limited concurrency. This is not ideal for many-user demand.
You can easily scale the number of celery workers you'd like to use with

```
$ docker-compose up -d --scale tb_c_worker=N tb_c_worker
```

where N is the number of workers you want, for example. The above note applies to each worker you start, however, and
each worker will consume RAM. You can also adjust the default number of workers defined by the variables at the top of
the `deploy.sh` script.

## Managing Django

If you'd like to manage the Django app (i.e. - run python manage.py ...), for example, to create an admin superuser, you can run commands in the _running_ app service (do this _after_ `docker-compose up`) as follows:

```
$ docker-compose exec app bash -c "python /usr/local/ASKCOS/askcos/manage.py createsuperuser"
```

In this case you'll be presented an interactive prompt to create a superuser with your desired credentials.

# Data migration to askcos-data
In the v0.4.1 release of ASKCOS, data and models have been migrated to a separate repository at https://github.com/ASKCOS/askcos-data.
The pre-built ASKCOS Docker image available from Docker Hub already contains the data and models.
For local use, you will need to clone the askcos-data repository separately:

```
$ cd ASKCOS/makeit
$ git clone https://github.com/ASKCOS/askcos-data data
```

# How to run individual modules
Many of the individual modules -- at least the ones that are the most interesting -- can be run "standalone". Examples of how to use them are often found in the ```if __name__ == '__main__'``` statement at the bottom of the script definitions. For example...

## Using the learned synthetic complexity metric (SCScore)
```makeit/prioritization/precursors/scscore.py```

## Obtaining a single-step retrosynthetic suggestion with consideration of chirality
```makeit/retrosynthetic/transformer.py```

## Finding recommended reaction conditions based on a trained neural network model
```makeit/synthetic/context/neuralnetwork.py```

## Using the template-free forward predictor
```makeit/synthetic/evaluation/template_free.py```

## Using the coarse "fast filter" (binary classifier) for evaluating reaction plausibility
```makeit/synthetic/evaluation/fast_filter.py```

## Integrated CASP tool
For the integrated synthesis planning tool at ```makeit/application/run.py```, there are several options available. The currently enabled options for the command-line tool can be found at ```makeit/utilities/io/arg_parser.py```. There are some options that are only available for the website and some that are only available for the command-line version. As an example of the former, the consideration of popular but non-buyable chemicals as suitable "leaf nodes" in the search. It is highly recommended to use the web interface when possible.
