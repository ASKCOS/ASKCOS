# ASKCOS Deployment - Containerized with Docker, deployed with docker-compose

### Prerequisites

 - If you're buidling the image from scratch, make sure git (and git lfs) is installed on your machine
 - Install Docker [OS specific instructions](https://docs.docker.com/install/)
 - Install docker-compose [installation instructions](https://docs.docker.com/compose/install/#install-compose)

### ASKCOS Image

The askcos image itself can be built using the Dockerfile in this repository `Make-It/Dockerfile`.

```bash
$ git clone https://github.com/connorcoley/Make-It  
$ cd Make-It/makeit/data  
$ git lfs pull  
$ cd ../../  
$ docker build -t askcos .
```

Alternatively pull from docker hub

```bash
$ docker login  
# login with credentials  
$ docker pull mefortunato/askcos  
$ docker tag mefortunato/askcos askcos # docker-compose expects there to be an image names askcos  
```

### Deploy with docker-compose

The `Make-It/deploy/docker-compose.yml` file contains the configuration to deploy the askcos stack with docker-compose. This requires that the askcos image is built (see previous step), and a few environment variables are set in the .env file.

```bash
$ cd deploy  
$ vi .env # fill out missing environment variables  
$ docker-compose up -d
```

The services will start in a detached state. You can view logs with `docker-compose logs [-f]`.

To stop the containers use `docker-compose down`

## Important Notes

#### First startup

The celery worker will take a few minutes to start up (possibly up to 5 minutes; it reads a lot of data into memory from disk). The web app itself will be ready before this, however upon the first get request (only the first) a few files will be read from disk, so expect a 10-15 second delay. You will know the celery workers are ready when the logs show messages such as `celery_worker_1  | ### TREE BUILDER WORKER STARTED UP ###` and `celery_worker_1  | ### NEURAL NETWORK CONTEXT RECOMMENDER STARTED UP ###`

#### Scaling workers

Only 2 worker are deployed by default. This may limit the number of concurrent tasks that can be performed. This is not ideal for multi-user demand. You can easily scale the number of celery workers you'd like to use with `docker-compose up -d --scale tb_worker=N` where N is te number of workers you want, for example. The above note applies to each worker you start, however, and each worker will comsume RAM.

#### Celery worker configuration

The current celery architecture may not be ideal for optimal concurrency. Currently two worker containers are spun up and they listen to separate sets of queues. This can be tinkered with in the docker-compose.yml file, to split more workers into more containers that listen to separate queues.

## Future Improvements

 - Optimize celery configuration (separate workers per queue(s) but reduce replicated setup, maybe ENV on deploy?)
 - Container orchestration with Kubernetes (will allow for distributed celery workers on multiple machines)
