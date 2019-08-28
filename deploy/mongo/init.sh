#!/bin/bash

echo "clearing"
mongo --username ${MONGO_USER} --password ${MONGO_PW} --authenticationDatabase admin ${MONGO_HOST}/askcos /init/clearDB.js

echo "importing"
gunzip -c /init/buyables.json.gz | mongoimport --host ${MONGO_HOST} --username ${MONGO_USER} --password ${MONGO_PW} --authenticationDatabase admin --db askcos --collection buyables --type json --jsonArray
gunzip -c /init/retro.templates.json.gz | mongoimport --host ${MONGO_HOST} --username ${MONGO_USER} --password ${MONGO_PW} --authenticationDatabase admin --db askcos --collection retro_templates --type json --jsonArray

echo "creating indexes"
mongo --username ${MONGO_USER} --password ${MONGO_PW} --authenticationDatabase admin ${MONGO_HOST}/askcos /init/createSmilesIndex.js
