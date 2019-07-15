#!/bin/bash

SKIP_SEED=false
SKIP_SSL=false
SKIP_MIGRATION=false

RANDOM_STRING=`cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 8 | head -n 1`

while [[ $# -gt 0 ]]
do
key="$1"

case $key in
    --skip-seed)
    SKIP_SEED=true
    shift # past argument
    ;;
    --skip-ssl)
    SKIP_SSL=true
    shift # past argument
    ;;
    --skip-migration)
    SKIP_MIGRATION=true
    shift # past argument
    ;;
    *)    # unknown option
    shift # past argument
    ;;
esac
done

echo "##########################"
echo "starting database services"
echo "##########################"
docker-compose up -d mysql mongo redis rabbit

if [ "$SKIP_SEED" = false ]; then
  echo "###############################"
  echo "seeding mongo database"
  echo "###############################"
  docker-compose up mongoseed
fi

if [ "$SKIP_SSL" = false ]; then
  echo "###############################"
  echo "creating SSL certificates"
  echo "###############################"
  openssl req   -new   -newkey rsa:4096   -days 3650   -nodes   -x509   -subj "/C=US/ST=MA/L=BOS/O=askcos/CN=askcos.$RANDOM_STRING.com"   -keyout askcos.ssl.key -out askcos.ssl.cert
fi

echo "#################################"
echo "starting web application services"
echo "#################################"
docker-compose up -d nginx app

echo "#######################"
echo "starting celery workers"
echo "#######################"
docker-compose up -d te_coordinator sc_coordinator ft_worker cr_coordinator cr_network_worker tb_coordinator_mcts tb_c_worker

if [ "$SKIP_MIGRATION" = false ]; then
  echo "#################"
  echo "migrating user db"
  echo "#################"
  docker-compose exec app bash -c "python /usr/local/ASKCOS/askcos/manage.py makemigrations main"
  docker-compose exec app bash -c "python /usr/local/ASKCOS/askcos/manage.py migrate"
fi
