#!/usr/bin/env bash

BACKUPFOLDER="backup/$(ls -t backup | head -1)"

RES=$(docker-compose exec app bash -c "python -c 'import makeit; print(makeit)'")
ASKCOSPATH=$(echo $RES | grep -o "/.*ASKCOS")

docker cp $BACKUPFOLDER/db.json deploy_app_1:$ASKCOSPATH/askcos/db.json
docker cp $BACKUPFOLDER/user_saves deploy_app_1:$ASKCOSPATH/makeit/data/

docker-compose exec --user=root app bash -c "chown -R askcos:askcos $ASKCOSPATH/askcos/db.json"
docker-compose exec --user=root app bash -c "chown -R askcos:askcos $ASKCOSPATH/makeit/data/user_saves"

docker-compose exec app bash -c "cd $ASKCOSPATH/askcos && echo 'from django.contrib.contenttypes.models import ContentType; ContentType.objects.all().delete()' | python manage.py shell"

docker cp transfer_results_to_mongo.py deploy_app_1:$ASKCOSPATH/askcos/
docker-compose exec app bash -c "cd $ASKCOSPATH/askcos && python transfer_results_to_mongo.py"

docker-compose exec app bash -c "cd $ASKCOSPATH/askcos && python manage.py loaddata db.json"

docker cp $BACKUPFOLDER/results.mongo deploy_mongo_1:/
docker-compose exec mongo bash -c "mongoimport -u askcos -p askcos --authenticationDatabase admin -d results -c results results.mongo"
