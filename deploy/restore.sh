#!/usr/bin/env bash

BACKUPFOLDER="backup/$(ls -t backup | head -1)"

RES=$(docker-compose exec app bash -c "python -c 'import makeit; print makeit'")
ASKCOSPATH=$(echo $RES | grep -o "/.*ASKCOS")

docker cp $BACKUPFOLDER/db.sqlite3 deploy_app_1:$ASKCOSPATH/askcos/db.sqlite3
docker cp $BACKUPFOLDER/user_saves deploy_app_1:$ASKCOSPATH/makeit/data/

docker-compose exec --user=root app bash -c "chown -R askcos:askcos $ASKCOSPATH/askcos/db.sqlite3"
docker-compose exec --user=root app bash -c "chown -R askcos:askcos $ASKCOSPATH/makeit/data/user_saves"
