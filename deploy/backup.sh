#!/usr/bin/env bash

BACKUPFOLDER="backup/$(date +%Y%m%d%s)"
mkdir -p $BACKUPFOLDER

RES=$(docker-compose exec app bash -c "python -c 'import makeit; print makeit'")
ASKCOSPATH=$(echo $RES | grep -o "/.*ASKCOS")

docker-compose exec app bash -c "cd $ASKCOSPATH/askcos && python manage.py dumpdata > db.json"

docker cp deploy_app_1:$ASKCOSPATH/askcos/db.json $BACKUPFOLDER
docker cp deploy_app_1:$ASKCOSPATH/makeit/data/user_saves $BACKUPFOLDER

docker-compose exec mongo bash -c "mongoexport -u askcos -p askcos --authenticationDatabase admin -d results -c results -o results.mongo"
docker cp deploy_mongo_1:/results.mongo $BACKUPFOLDER
