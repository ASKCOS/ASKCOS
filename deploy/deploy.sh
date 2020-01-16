#!/usr/bin/env bash

################################################################################
#
#   ASKCOS Deployment Utilities
#
#    ~ To streamline deployment commands ~
#
################################################################################

set -e  # exit with nonzero exit code if anything fails

usage() {
  echo
  echo "Deployment Utilities for ASKCOS"
  echo
  echo "Specify a task to perform, along with any desired options."
  echo
  echo "Valid commands:"
  echo "    deploy:            performs initial deployment steps using https"
  echo "    deploy-http:       performs initial deployment steps using http"
  echo "    update:            update an existing deployment"
  echo "    seed-db:           seed mongo database with data"
  echo "    migrate:           perform user database migrations"
  echo "    start:             (re)start an existing deployment"
  echo "    stop:              stop a currently running deployment"
  echo "    clean:             stop and remove a currently running deployment"
  echo
  echo "Optional arguments:"
  echo "    -f,--compose-file  specify docker-compose file(s) for deployment"
  echo "    -v,--version       specify desired version for updating a deployment"
#  echo "    -t,--templates     template data for reseeding mongo database"
#  echo "    -b,--buyables      buyables data for reseeding mongo database"
  echo "    -d,--dev           use docker-compose configuration for development (fewer workers)"
  echo
  echo "Examples:"
  echo "    ./deploy.sh deploy -f docker-compose.yml"
  echo "    ./deploy.sh update -v x.y.z"
#  echo "    ./deploy.sh seed-db -t templates.json -b buyables.json"
  echo "    ./deploy.sh clean"
  echo
}

# Default argument values
COMPOSE_FILE=""
VERSION="0.4.1"
TEMPLATES=""
BUYABLES=""

COMMANDS=""
while (( "$#" )); do
  case "$1" in
    -h|--help|help)
      usage
      exit
      ;;
    -f|--compose-file)
      if [ -z "$COMPOSE_FILE" ]; then
        COMPOSE_FILE=$2
      else
        COMPOSE_FILE=$COMPOSE_FILE:$2
      fi
      shift 2
      ;;
    -d|--dev)
      COMPOSE_FILE="docker-compose.yml:docker-compose.dev.yml"
      shift 1
      ;;
    -v|--version)
      VERSION=$2
      shift 2
      ;;
    -t|--templates)
      TEMPLATES=$2  # TODO: This is not used anywhere
      shift 2
      ;;
    -b|--buyables)
      BUYABLES=$2  # TODO: This is not used anywhere
      shift 2
      ;;
    --) # end argument parsing
      shift
      break
      ;;
    -*) # any other flag
      echo "Error: Unsupported flag $1" >&2  # print to stderr
      exit 1
      ;;
    *) # preserve positional arguments
      COMMANDS="$COMMANDS $1"
      shift
      ;;
  esac
done

# Set positional arguments in their proper place
eval set -- "$COMMANDS"

# Export VERSION and COMPOSE_FILE so they're available to docker-compose
export VERSION
export COMPOSE_FILE

# Define various functions
clean-static() {
  echo "Cleaning up old static file volume..."
  docker-compose stop app nginx
  docker-compose rm -f app nginx
  docker volume rm deploy_staticdata
  echo "Clean up complete."
  echo
}

start-db-services() {
  echo "Starting database services..."
  docker-compose up -d mysql mongo redis rabbit
  sleep 1
  echo "Start up complete."
  echo
}

seed-db() {
  echo "Seeding mongo database..."
  docker-compose exec app python -c "from askcos_site.main.db import seed_mongo_db;seed_mongo_db(reactions=False, chemicals=False)"
  echo "Seeding complete."
  echo
}

copy-http-conf() {
  echo "Using http nginx configuration."
  cp nginx.http.conf nginx.conf
}

copy-https-conf() {
  echo "Using https nginx configuration."
  cp nginx.https.conf nginx.conf
  echo
}

create-ssl() {
  if [ ! -f "askcos.ssl.cert" ]; then
    echo "Creating SSL certificates."
    openssl req -new -newkey rsa:4096 -days 3650 -nodes -x509 -subj "/C=US/ST=MA/L=BOS/O=askcos/CN=askcos.$RANDOM.com" -keyout askcos.ssl.key -out askcos.ssl.cert
    echo
  fi
}

start-web-services() {
  echo "Starting web services..."
  docker-compose up -d nginx app
  echo "Start up complete."
  echo
}

start-tf-server() {
  echo "Starting tensorflow serving worker..."
  docker-compose up -d template-relevance-reaxys
  echo "Start up complete."
  echo
}

start-celery-workers() {
  echo "Starting celery workers..."
  docker-compose up -d te_coordinator sc_coordinator ft_worker cr_coordinator cr_network_worker tb_coordinator_mcts \
                       tb_c_worker tb_c_worker_preload sites_worker impurity_worker atom_mapping_worker
  echo "Start up complete."
  echo
}

migrate() {
  echo "Migrating user database..."
  docker-compose exec app bash -c "python /usr/local/ASKCOS/askcos/manage.py makemigrations main"
  docker-compose exec app bash -c "python /usr/local/ASKCOS/askcos/manage.py migrate"
  echo "Migration complete."
  echo
}

# Handle positional arguments, which should be commands
if [ $# -eq 0 ]; then
  # No arguments
  echo "Must provide a valid task, e.g. deploy|update|migrate."
  echo "See 'deploy.sh help' for more options."
  exit 1;
else
  for arg in "$@"
  do
    case "$arg" in
      clean-static | start-db-services | seed-db | copy-http-conf | copy-https-conf | create-ssl | \
      start-web-services | start-tf-server | start-celery-workers | migrate)
        # This is a defined function, so execute it
        $arg
        ;;
      deploy)
        # Normal first deployment, do everything
        copy-https-conf
        create-ssl
        start-db-services
        start-web-services
        seed-db  # Must occur after starting app
        start-tf-server
        start-celery-workers
        migrate
        ;;
      deploy-http)
        # Deploy with http, only difference is ssl cert creation and nginx conf
        copy-http-conf
        start-db-services
        start-web-services
        seed-db  # Must occur after starting app
        start-tf-server
        start-celery-workers
        migrate
        ;;
      update)
        # Update an existing configuration, database seeding is not performed
        docker pull registry.gitlab.com/mlpds_mit/askcos/askcos:"$VERSION"
        clean-static
        start-db-services
        start-web-services
        start-tf-server
        start-celery-workers
        migrate
        ;;
      start)
        # (Re)start existing deployment
        start-db-services
        start-web-services
        start-tf-server
        start-celery-workers
        ;;
      stop)
        # Stop currently running containers
        docker-compose stop
        ;;
      clean)
        # Clean up current deployment
        docker-compose down -v  # make sure there's a prompt for confirmation
        ;;
      *)
        echo "Error: Unsupported command $1" >&2  # print to stderr
        exit 1;
    esac
  done
fi
