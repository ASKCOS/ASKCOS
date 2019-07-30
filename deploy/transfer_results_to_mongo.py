import sys
import json
from pymongo import MongoClient
from makeit import global_config as gc

if len(sys.argv) == 2:
    db_json = sys.argv[1]
else:
    db_json = 'db.json'

client = MongoClient(
    gc.MONGO['path'],
    gc.MONGO['id'],
    connect=gc.MONGO['connect']
)
results_db = client['results']
results_collection = results_db['results']

with open(db_json) as f:
    j = json.load(f)

saved_results = [thing for thing in j if thing['model'] == 'main.savedresults']

for result in saved_results:
    if not result['fields'].get('fpath'): continue
    path = result['fields'].pop('fpath')
    result_id = path.split('/')[-1].split('.')[0]
    result['fields']['result_id'] = result_id
    result['fields']['result_type'] = 'html'
    with open(path) as f:
        html = f.read()
    results_collection.update(
        {'_id': result_id},
        {'result': html},
        upsert=True
    )

with open(db_json, 'w') as f:
    json.dump(j, f)
