#!/bin/bash

SKIP_SEED=false
SKIP_SSL=false
SKIP_MIGRATION=false
SKIP_HTTPS=false

while [[ $# -gt 0 ]]
do
key="$1"

case $key in
    --skip-seed)
    SKIP_SEED=true
    shift # past argument
    ;;
    --skip-https)
    SKIP_HTTPS=true
    SKIP_SSL=true
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

echo ""
echo "##################################"
echo "cleaning up old static file volume"
echo "##################################"
docker-compose stop app nginx
docker-compose rm -f app nginx
docker volume rm deploy_staticdata

echo ""
echo "##########################"
echo "starting database services"
echo "##########################"
docker-compose up -d mysql mongo redis rabbit
sleep 1

if [ "$SKIP_SEED" = false ]; then
  echo ""
  echo "###############################"
  echo "seeding mongo database"
  echo "###############################"
  cd mongo
  docker run -it --rm --network deploy_default --env-file ../.env -v ${PWD}:/init mongo bash /init/init.sh
  cd ../

fi

if [ "$SKIP_HTTPS" = true ]; then
  echo ""
  echo "###############################"
  echo "not using https"
  echo "###############################"
  cp nginx.http.conf nginx.conf
else
  echo ""
  echo "###############################"
  echo "using https"
  echo "###############################"
  cp nginx.https.conf nginx.conf
fi

if [ "$SKIP_SSL" = false ]; then
  echo ""
  echo "###############################"
  echo "creating SSL certificates"
  echo "###############################"
  openssl req   -new   -newkey rsa:4096   -days 3650   -nodes   -x509   -subj "/C=US/ST=MA/L=BOS/O=askcos/CN=askcos.$RANDOM.com"   -keyout askcos.ssl.key -out askcos.ssl.cert
fi

echo ""
echo "#################################"
echo "starting web application services"
echo "#################################"
docker-compose up -d nginx app
docker-compose exec app bash -c "python /usr/local/ASKCOS/askcos/manage.py collectstatic --noinput"

echo ""
echo "#######################"
echo "starting celery workers"
echo "#######################"
docker-compose up -d te_coordinator sc_coordinator ft_worker cr_coordinator cr_network_worker tb_coordinator_mcts tb_c_worker sites_worker
docker-compose up -d --scale tb_coordinator_mcts=2

if [ "$SKIP_MIGRATION" = false ]; then
  echo ""
  echo "#################"
  echo "migrating user db"
  echo "#################"
  docker-compose exec app bash -c "python /usr/local/ASKCOS/askcos/manage.py makemigrations main"
  docker-compose exec app bash -c "python /usr/local/ASKCOS/askcos/manage.py migrate"
fi
