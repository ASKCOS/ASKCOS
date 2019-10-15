#!/bin/bash

POD=$(kubectl get pod -l pod=mongo -o jsonpath="{.items[0].metadata.name}")
kubectl exec -it $POD -- bash -c "mkdir /init"
kubectl cp buyables.json.gz $POD:/init
kubectl cp retro.templates.json.gz $POD:/init
kubectl cp clearDB.js $POD:/init
kubectl cp createSmilesIndex.js $POD:/init
kubectl cp init.sh $POD:/init
kubectl exec -it $POD -- bash /init/init.sh