#!/usr/bin/env bash

BACKUPFOLDER="backup/$(date +%Y%m%d%s)"
mkdir -p $BACKUPFOLDER

RES=$(docker-compose exec app bash -c "python -c 'import makeit; print makeit'")
ASKCOSPATH=$(echo $RES | grep -o "/.*ASKCOS")

docker cp deploy_app_1:$ASKCOSPATH/askcos/db.sqlite3 $BACKUPFOLDER
docker cp deploy_app_1:$ASKCOSPATH/makeit/data/user_saves $BACKUPFOLDER
