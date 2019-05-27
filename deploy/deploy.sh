#!/bin/bash

echo "##########################"
echo "starting database services"
echo "##########################"
docker-compose up -d mysql mongo redis rabbit

echo "###############################"
echo "migrating and seeding databases"
echo "###############################"
docker-compose up -d mongoseed

echo "#################################"
echo "starting web application services"
echo "#################################"
docker-compose up -d nginx app

echo "#######################"
echo "starting celery workers"
echo "#######################"
docker-compose up -d te_coordinator sc_coordinator ft_worker cr_coordinator cr_network_worker tb_coordinator_mcts tb_c_worker

echo "#################"
echo "migrating user db"
echo "#################"
docker-compose exec app bash -c "python /usr/local/ASKCOS/askcos/manage.py makemigrations main"
docker-compose exec app bash -c "python /usr/local/ASKCOS/askcos/manage.py migrate"
