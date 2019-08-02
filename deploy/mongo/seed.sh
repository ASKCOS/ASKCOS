#!/bin/bash

if ls /init/*.json.gz &>/dev/null
then
  echo "gunzipping"
  gunzip /init/*.json.gz
else
  echo "clearing"
  mongo --username ${MONGO_USER} --password ${MONGO_PW} --authenticationDatabase admin ${MONGO_HOST}/askcos /init/clearDB.js
fi

echo "importing"
mongoimport --host ${MONGO_HOST} --username ${MONGO_USER} --password ${MONGO_PW} --authenticationDatabase admin --db askcos --collection buyables --type json --jsonArray --file /init/buyables.json
mongoimport --host ${MONGO_HOST} --username ${MONGO_USER} --password ${MONGO_PW} --authenticationDatabase admin --db askcos --collection retro_templates --type json --jsonArray --file /init/retro.templates.json

echo "creating indexes"
mongo --username ${MONGO_USER} --password ${MONGO_PW} --authenticationDatabase admin ${MONGO_HOST}/askcos /init/createSmilesIndex.js
