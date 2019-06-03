
# ASKCOS:
Software package for the prediction of feasible synthetic routes towards a desired compound and associated tasks related to synthesis planning. Originally developed under the DARPA Make-It program and now being developed under the [MLPDS Consortium](http://mlpds.mit.edu). The tools are deployed as a Django webapp. A version similar but not identical to this code is currently hosted at [askcos.mit.edu](http://askcos.mit.edu).

Please note that the MPL 2.0 license for this repository does not apply to the data and trained models. The data and trained models are released under CC BY-NC-SA (i.e., are for noncommercial use only).

Contributors include Connor Coley, Mike Fortunato, Hanyu Gao, and Pieter Plehiers.



# Quick start using Google Cloud

```
# (1) Create a Google Cloud instance 
#     - 8 vCPUs, 52 GB memory is the maximum for their free trial
#     - select Ubuntu 18.04 LTS Minimal
#     - upgrade to a 100 GB disk
#     - allow HTTP traffic

# (2) Install docker
#     - https://docs.docker.com/install/linux/docker-ce/ubuntu/
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
sudo curl -L "https://github.com/docker/compose/releases/download/1.24.0/docker-compose-$(uname -s)-$(uname -m)" -o /usr/local/bin/docker-compose
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
git clone https://github.com/connorcoley/ASKCOS ASKCOS
cd ASKCOS
git lfs pull 

# (6) Build & run
cd deploy
docker build -t askcos .. # build
docker-compose up -d      # start containers (detached)
docker-compose logs -f    # start tailing logs (can CTRL+C to exit)

# (7) Navigate to your instance's external IP 
#     - note that it may take ~5 minutes for the retro transformer workers to start up
#     - you can check the status of their startup by looking at "server status"
#     - the first request to a website process may take ~10 seconds
#     - the first request to a retro transform worker may take ~5-10 seconds
#     - the first request to the forward predictor may take ~60 seconds
```

# Installation with Docker

### Prerequisites

 - If you're buidling the image from scratch, make sure git (and git lfs) is installed on your machine
 - Install Docker [OS specific instructions](https://docs.docker.com/install/)
 - Install docker-compose [installation instructions](https://docs.docker.com/compose/install/#install-compose)


### Building the ASKCOS Image

The askcos image itself can be built using the Dockerfile in this repository `ASKCOS/Dockerfile`.

```bash
$ git clone https://github.com/connorcoley/ASKCOS  
$ cd ASKCOS/makeit/data  
$ git lfs pull  
$ cd ../../  
$ docker build -t askcos .
```


### Deploy with docker-compose

The `Make-It/deploy/docker-compose.yml` file contains the configuration to deploy the askcos stack with docker-compose. This requires that the askcos image is built (see previous step), and a few environment variables are set in the .env file. The default ENV values will work, but it is better to set `CURRENT_HOST` to the IP address of the machine you are deploying on, and to set the MongoDB credentials if you have access.

```bash
$ cd deploy  
$ docker-compose up -d
```

The services will start in a detached state. You can view logs with `docker-compose logs [-f]`.

To stop the containers use `docker-compose stop`. To restart the containers use `docker-compose start`. To completely delete the containers and volumes use `docker-compose down -v` (this deletes user database and saves; read section about backing up data first).

### Managing Django

If you'd like to manage the Django app (i.e. - run python manage.py ...), for example, to create an admin superuser, you can run commands in the _running_ app service (do this _after_ `docker-compose up`) as follows:

`docker-compose exec app bash -c "python /usr/local/ASKCOS/askcos/manage.py createsuperuser"`

In this case you'll be presented an interactive prompt to create a superuser with your desired credentials.

## Important Notes

#### Recommended hardward

We recommend running this code on a machine with at least 8 compute cores (16 preferred) and 64 GB RAM (128 GB preferred)

#### First startup

The celery worker will take a few minutes to start up (possibly up to 5 minutes; it reads a lot of data into memory from disk). The web app itself will be ready before this, however upon the first get request (only the first for each process) a few files will be read from disk, so expect a 10-15 second delay.

#### Scaling workers

Only 1 worker per queue is deployed by default with limited concurrency. This is not ideal for many-user demand. You can easily scale the number of celery workers you'd like to use with `docker-compose up -d --scale tb_c_worker=N` where N is the number of workers you want, for example. The above note applies to each worker you start, however, and each worker will consume RAM.


### Dependencies
The code has primarily been developed for Python 2.7.6 on Ubuntu 16.04. However, we have made an effort to make it work on Python 3.6.1 as well (tested on macOS 10.13.3). It heavily relies on RDKit.

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

#### Using the tree builder to find full pathways
```makeit/retrosynthetic/mcts/tree_builder.py```

#### Integrated CASP tool
For the integrated synthesis planning tool at ```makeit/application/run.py```, there are several options available. The currently enabled options for the command-line tool can be found at ```makeit/utilities/io/arg_parser.py```. There are some options that are only available for the website and some that are only available for the command-line version. As an example of the former, the consideration of popular but non-buyable chemicals as suitable "leaf nodes" in the search. It is highly recommended to use the web interface when possible.
