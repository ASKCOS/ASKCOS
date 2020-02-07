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
  echo "    deploy:                   performs initial deployment steps using https"
  echo "    deploy-http:              performs initial deployment steps using http"
  echo "    update:                   update an existing deployment"
  echo "    seed-db:                  seed mongo database with data"
  echo "    migrate:                  perform user database migrations"
  echo "    start:                    (re)start an existing deployment"
  echo "    stop:                     stop a currently running deployment"
  echo "    clean:                    stop and remove a currently running deployment"
  echo
  echo "Optional arguments:"
  echo "    -f,--compose-file         specify docker-compose file(s) for deployment"
  echo "    -v,--version              specify desired version for updating a deployment"
  echo "    -b,--buyables             buyables data for reseeding mongo database"
  echo "    -c,--chemicals            chemicals data for reseeding mongo database"
  echo "    -x,--reactions            reactions data for reseeding mongo database"
  echo "    -r,--retro-templates      retrosynthetic template data for reseeding mongo database"
  echo "    -t,--forward-templates    forward template data for reseeding mongo database"
  echo "    -d,--dev                  use docker-compose configuration for development (fewer workers)"
  echo
  echo "Examples:"
  echo "    bash deploy.sh deploy -f docker-compose.yml"
  echo "    bash deploy.sh update -v x.y.z"
  echo "    bash deploy.sh seed-db -r retro-templates.json.gz -b buyables.json.gz"
  echo "    bash deploy.sh clean"
  echo
}

# Worker scales (i.e. number of celery workers)
n_te_coordinator=1       # Tree evaluation coordinator
n_sc_coordinator=1       # Scoring coordinator
n_ft_worker=1            # Forward transformer worker
n_cr_network_worker=1    # Context recommender neural network worker
n_tb_coordinator_mcts=2  # Tree builder coordinator
n_tb_c_worker=12         # Tree builder chiral worker
n_tb_c_worker_preload=1  # Tree builder chiral worker with template preloading
n_sites_worker=1         # Site selectivity worker
n_impurity_worker=1      # Impurity worker
n_atom_mapping_worker=1  # Atom mapping worker

# Default argument values
COMPOSE_FILE=""
VERSION=""
BUYABLES=""
CHEMICALS=""
REACTIONS=""
RETRO_TEMPLATES=""
FORWARD_TEMPLATES=""

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
      n_tb_coordinator_mcts=1  # Tree builder coordinator
      n_tb_c_worker=1          # Tree builder chiral worker
      shift 1
      ;;
    -v|--version)
      VERSION=$2
      shift 2
      ;;
    -b|--buyables)
      BUYABLES=$2
      shift 2
      ;;
    -c|--chemicals)
      CHEMICALS=$2
      shift 2
      ;;
    -x|--reactions)
      REACTIONS=$2
      shift 2
      ;;
    -r|--retro-templates)
      RETRO_TEMPLATES=$2
      shift 2
      ;;
    -t|--forward-templates)
      FORWARD_TEMPLATES=$2
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

set-db-defaults() {
  # Set default values for seeding database if values are not already defined
  BUYABLES=${BUYABLES:-default}
  RETRO_TEMPLATES=${RETRO_TEMPLATES:-default}
  FORWARD_TEMPLATES=${FORWARD_TEMPLATES:-default}
  CHEMICALS=${CHEMICALS:-default}
}

run-mongo-js() {
  # arg 1 is js command
  docker-compose exec mongo bash -c 'mongo --username ${MONGO_USER} --password ${MONGO_PW} --authenticationDatabase admin ${MONGO_HOST}/askcos --quiet --eval '"'$1'"
}

seed-db-collection() {
  # arg 1 is collection name
  # arg 2 is file path
  # arg 3 is a flag to pass to docker-compose exec, e.g. -d to detach
  docker-compose exec $3 mongo bash -c 'gunzip -c '$2' | mongoimport --host ${MONGO_HOST} --username ${MONGO_USER} --password ${MONGO_PW} --authenticationDatabase admin --db askcos --collection '$1' --type json --jsonArray --drop'
}

seed-db() {
  if [ -z "$BUYABLES" ] && [ -z "$CHEMICALS" ] && [ -z "$REACTIONS" ] && [ -z "$RETRO_TEMPLATES" ] && [ -z "$FORWARD_TEMPLATES" ]; then
    echo "Nothing to seed!"
    echo "Example usages:"
    echo "    bash deploy.sh seed-db -r default                  seed only the default retro templates"
    echo "    bash deploy.sh seed-db -r <templates.json.gz>      seed retro templates from local file <templates.json.gz>"
    echo "    bash deploy.sh set-db-defaults seed-db             seed all default collections"
    return
  fi

  echo "Seeding mongo database..."

  if [ "$BUYABLES" = "default" ]; then
    echo "Loading default buyables data in background..."
    buyables_file="/data/app/buyables/buyables.json.gz"
    seed-db-collection buyables "$buyables_file" -d
    run-mongo-js 'db.buyables.createIndex({smiles: "text"})'
  elif [ -n "$BUYABLES" ]; then
    echo "Loading buyables data from $BUYABLES in background..."
    buyables_file="/data/app/buyables/$(basename $BUYABLES)"
    docker cp "$BUYABLES" deploy_mongo_1:"$buyables_file"
    run-mongo-js "db.buyables.remove({})"
    seed-db-collection buyables "$buyables_file" -d
    run-mongo-js 'db.buyables.createIndex({smiles: "text"})'
  fi

  if [ "$CHEMICALS" = "default" ]; then
    echo "Loading default chemicals data in background..."
    chemicals_file="/data/app/historian/chemicals.json.gz"
    seed-db-collection chemicals "$chemicals_file" -d
    run-mongo-js 'db.chemicals.createIndex({smiles: "hashed"})'
  elif [ -n "$CHEMICALS" ]; then
    echo "Loading chemicals data from $CHEMICALS in background..."
    chemicals_file="/data/app/historian/$(basename $CHEMICALS)"
    docker cp "$CHEMICALS" deploy_mongo_1:"$chemicals_file"
    seed-db-collection chemicals "$chemicals_file" -d
    run-mongo-js 'db.chemicals.createIndex({smiles: "hashed"})'
  fi

  if [ "$REACTIONS" = "default" ]; then
    echo "Loading default reactions data in background..."
    reactions_file="/data/app/historian/reactions.json.gz"
    seed-db-collection reactions "$reactions_file" -d
  elif [ -n "$REACTIONS" ]; then
    echo "Loading reactions data from $REACTIONS in background..."
    reactions_file="/data/app/historian/$(basename $REACTIONS)"
    docker cp "$REACTIONS" deploy_mongo_1:"$reactions_file"
    seed-db-collection reactions "$reactions_file" -d
  fi

  if [ "$RETRO_TEMPLATES" = "default" ]; then
    echo "Loading default retrosynthetic templates..."
    retro_file="/data/app/templates/retro.templates.json.gz"
    seed-db-collection retro_templates "$retro_file"
    run-mongo-js 'db.retro_templates.createIndex({index: 1})'
  elif [ -n "$RETRO_TEMPLATES" ]; then
    echo "Loading retrosynthetic templates from $RETRO_TEMPLATES ..."
    retro_file="/data/app/templates/$(basename $RETRO_TEMPLATES)"
    docker cp "$RETRO_TEMPLATES" deploy_mongo_1:"$retro_file"
    seed-db-collection retro_templates "$retro_file"
    run-mongo-js 'db.retro_templates.createIndex({index: 1})'
  fi

  if [ "$FORWARD_TEMPLATES" = "default" ]; then
    echo "Loading default forward templates..."
    forward_file="/data/app/templates/forward.templates.json.gz"
    seed-db-collection forward_templates "$forward_file"
  elif [ -n "$FORWARD_TEMPLATES" ]; then
    echo "Loading forward templates from $FORWARD_TEMPLATES ..."
    forward_file="/data/app/templates/$(basename $FORWARD_TEMPLATES)"
    docker cp "$FORWARD_TEMPLATES" deploy_mongo_1:"$forward_file"
    seed-db-collection forward_templates "$forward_file"
  fi

  echo "Seeding complete."
  echo
}

count-mongo-docs() {
  echo "Buyables collection:          $(run-mongo-js "db.buyables.countDocuments({})" | tr -d '\r') / 106750 expected (default)"
  echo "Chemicals collection:         $(run-mongo-js "db.chemicals.countDocuments({})" | tr -d '\r') / 17562038 expected (default)"
  echo "Reactions collection:         $(run-mongo-js "db.reactions.countDocuments({})" | tr -d '\r') / 0 expected (default)"
  echo "Retro template collection:    $(run-mongo-js "db.retro_templates.countDocuments({})" | tr -d '\r') / 163723 expected (default)"
  echo "Forward template collection:  $(run-mongo-js "db.forward_templates.countDocuments({})" | tr -d '\r') / 17089 expected (default)"
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
  docker-compose up -d --remove-orphans nginx app
  echo "Start up complete."
  echo
}

start-tf-server() {
  echo "Starting tensorflow serving worker..."
  docker-compose up -d --remove-orphans template-relevance-reaxys fast-filter
  echo "Start up complete."
  echo
}

start-celery-workers() {
  echo "Starting celery workers..."
  docker-compose up -d --scale te_coordinator=$n_te_coordinator \
                       --scale sc_coordinator=$n_sc_coordinator \
                       --scale ft_worker=$n_ft_worker \
                       --scale cr_network_worker=$n_cr_network_worker \
                       --scale tb_coordinator_mcts=$n_tb_coordinator_mcts \
                       --scale tb_c_worker=$n_tb_c_worker \
                       --scale tb_c_worker_preload=$n_tb_c_worker_preload \
                       --scale sites_worker=$n_sites_worker \
                       --scale impurity_worker=$n_impurity_worker \
                       --scale atom_mapping_worker=$n_atom_mapping_worker \
                       --remove-orphans \
                       te_coordinator sc_coordinator ft_worker cr_network_worker tb_coordinator_mcts \
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
      start-web-services | start-tf-server | start-celery-workers | migrate | set-db-defaults | count-mongo-docs)
        # This is a defined function, so execute it
        $arg
        ;;
      deploy)
        # Normal first deployment, do everything
        copy-https-conf
        create-ssl
        start-db-services
        start-web-services
        set-db-defaults
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
        set-db-defaults
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
        echo "This will stop and remove all containers and also remove all data volumes. Are you sure you want to continue?"
        read -rp 'Continue (y/N): ' response
        case "$response" in
          [Yy] | [Yy][Ee][Ss])
            echo "Cleaning deployment."
            docker-compose down -v --remove-orphans
            ;;
          *)
            echo "Doing nothing."
            ;;
        esac
        ;;
      *)
        echo "Error: Unsupported command $1" >&2  # print to stderr
        exit 1;
    esac
  done
fi
