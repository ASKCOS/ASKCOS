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
  echo "    ./deploy.sh deploy -f docker-compose.yml"
  echo "    ./deploy.sh update -v x.y.z"
  echo "    ./deploy.sh seed-db -r retro-templates.json -b buyables.json"
  echo "    ./deploy.sh clean"
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
VERSION="0.4.1"
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

set_db_defaults() {
  # Set default values for seeding database if values are not already defined
  BUYABLES=${BUYABLES:-default}
  RETRO_TEMPLATES=${RETRO_TEMPLATES:-default}
  FORWARD_TEMPLATES=${FORWARD_TEMPLATES:-default}
  CHEMICALS=${CHEMICALS:-default}
}

seed-db() {
  echo "Seeding mongo database..."
  MAKEIT_PATH=$(docker-compose exec app bash -c "python -c 'import makeit; print(makeit.__file__.split(\"/__\")[0])'" | tr -d '\r')

  if [ "$BUYABLES" = "default" ]; then
    echo "Loading default buyables data..."
    docker-compose exec app python -c "from askcos_site.main.db import seed_mongo_db; seed_mongo_db(buyables=True, chemicals=False, reactions=False, retro_templates=False, forward_templates=False)"
  elif [ -n "$BUYABLES" ]; then
    echo "Loading buyables data from $BUYABLES ..."
    docker cp "$BUYABLES" deploy_app_1:"$MAKEIT_PATH/data/buyables/$(basename $BUYABLES)"
    docker-compose exec app python -c "from askcos_site.main.db import seed_buyables; seed_buyables('$MAKEIT_PATH/data/buyables/$(basename $BUYABLES)')"
  fi

  if [ "$CHEMICALS" = "default" ]; then
    echo "Loading default chemicals data..."
    docker-compose exec app python -c "from askcos_site.main.db import seed_mongo_db; seed_mongo_db(buyables=False, chemicals=True, reactions=False, retro_templates=False, forward_templates=False)"
  elif [ -n "$CHEMICALS" ]; then
    echo "Loading chemicals data from $CHEMICALS ..."
    docker cp "$CHEMICALS" deploy_app_1:"$MAKEIT_PATH/data/historian/$(basename $CHEMICALS)"
    docker-compose exec app python -c "from askcos_site.main.db import seed_chemicals; seed_chemicals('$MAKEIT_PATH/data/historian/$(basename $CHEMICALS)')"
  fi

  if [ "$REACTIONS" = "default" ]; then
    echo "Loading default reactions data..."
    docker-compose exec app python -c "from askcos_site.main.db import seed_mongo_db; seed_mongo_db(buyables=False, chemicals=False, reactions=True, retro_templates=False, forward_templates=False)"
  elif [ -n "$REACTIONS" ]; then
    echo "Loading reactions data from $REACTIONS ..."
    docker cp "$REACTIONS" deploy_app_1:"$MAKEIT_PATH/data/historian/$(basename $REACTIONS)"
    docker-compose exec app python -c "from askcos_site.main.db import seed_reactions; seed_reactions('$MAKEIT_PATH/data/historian/$(basename $REACTIONS)')"
  fi

  if [ "$RETRO_TEMPLATES" = "default" ]; then
    echo "Loading default retrosynthetic templates..."
    docker-compose exec app python -c "from askcos_site.main.db import seed_mongo_db; seed_mongo_db(buyables=False, chemicals=False, reactions=False, retro_templates=True, forward_templates=False)"
  elif [ -n "$RETRO_TEMPLATES" ]; then
    echo "Loading retrosynthetic templates from $RETRO_TEMPLATES ..."
    docker cp "$RETRO_TEMPLATES" deploy_app_1:"$MAKEIT_PATH/data/templates/$(basename $RETRO_TEMPLATES)"
    docker-compose exec app python -c "from askcos_site.main.db import seed_retro_templates; seed_retro_templates('$MAKEIT_PATH/data/templates/$(basename $RETRO_TEMPLATES)')"
  fi

  if [ "$FORWARD_TEMPLATES" = "default" ]; then
    echo "Loading default forward templates..."
    docker-compose exec app python -c "from askcos_site.main.db import seed_mongo_db; seed_mongo_db(buyables=False, chemicals=False, reactions=False, retro_templates=False, forward_templates=True)"
  elif [ -n "$FORWARD_TEMPLATES" ]; then
    echo "Loading forward templates from $FORWARD_TEMPLATES ..."
    docker cp "$FORWARD_TEMPLATES" deploy_app_1:"$MAKEIT_PATH/data/templates/$(basename $FORWARD_TEMPLATES)"
    docker-compose exec app python -c "from askcos_site.main.db import seed_forward_templates; seed_forward_templates('$MAKEIT_PATH/data/templates/$(basename $FORWARD_TEMPLATES)')"
  fi

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
  docker-compose up -d template-relevance-reaxys fast-filter
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
        set_db_defaults
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
        set_db_defaults
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
