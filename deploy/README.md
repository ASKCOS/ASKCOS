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
# login with credentials; must have access to private repo  
$ docker pull mefortunato/askcos  
$ docker tag mefortunato/askcos askcos # docker-compose expects there to be an image names askcos  
```

### Deploy with docker-compose

The `Make-It/deploy/docker-compose.yml` file contains the configuration to deploy the askcos stack with docker-compose. This requires that the askcos image is built (see previous step), and a few environment variables are set in the .env file. The default ENV values will work, but it is better to set `CURRENT_HOST` to the IP address of the machine you are deploying on, and to set the MongoDB credentials if you have access.

```bash
$ cd deploy  
$ docker-compose up -d
```

The services will start in a detached state. You can view logs with `docker-compose logs [-f]`.

To stop the containers use `docker-compose stop`. To restart the containers use `docker-compose start`. To completely delete the containers and volumes use `docker-compose down -v` (this deletes user database and saves).

## Important Notes

#### First startup

The celery worker will take a few minutes to start up (possibly up to 5 minutes; it reads a lot of data into memory from disk). The web app itself will be ready before this, however upon the first get request (only the first for each process) a few files will be read from disk, so expect a 10-15 second delay. You will know the celery workers are ready when the logs show messages such as `tb_c_worker_1  | ### TREE BUILDER WORKER STARTED UP ###`.

#### Scaling workers

Only 1 worker per queue is deployed by default with limited concurrency. This is not ideal for many-user demand. You can easily scale the number of celery workers you'd like to use with `docker-compose up -d --scale tb_c_worker=N` where N is te number of workers you want, for example. The above note applies to each worker you start, however, and each worker will comsume RAM.

## Future Improvements

 - Container orchestration with Kubernetes (will allow for distributed celery workers on multiple machines)
