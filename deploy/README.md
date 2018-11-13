# ASKCOS - Containerized with Docker, deployed with docker-compose

### Prerequisites

 - Make sure git is installed on your machine
 - You may need git-lfs to acquire Make-It/makeit/data
 - Install Docker [OS specific instructions](https://docs.docker.com/install/)
 - Install docker-compose [installation instructions](https://docs.docker.com/compose/install/#install-compose)

### ASKCOS Image

The askcos image itself can be built using the Dockerfile in this repository `Make-It/Dockerfile`.

```bash
# from directory with Dockerfile
$ docker build -t askcos .
```


__*NOTE*__: The docker image is built excluding the Make-It/makeit/data directory. This data directory is intended to be mounted into the container when it is run:

```bash
$ docker run -it --name=askcos -v /full/path/to/Make-It/makeit/data:/home/askcos/ASKCOS/makeit/data askcos
```


### Deploy with docker-compose

The `Make-It/deploy/docker-compose.yml` file contains the configuration to deploy the askcos stack with docker-compose. This requires that the askcos image is built, the makeit/data directory is intact (must use git lfs to pull large data files), and a few environment variables are set in the .env file.

```bash
$ cd deploy

$ ln -s ../makeit/data # this allows the data directory to be mounted into container during deployment

$ vi .env # fill out missing environment variables

$ docker-compose up -d
```

The services will start in a detached state. You can view logs with `docker-compose logs [-f]`.

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
