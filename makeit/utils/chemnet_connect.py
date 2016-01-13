# Import relevant packages
from pymongo import MongoClient    # mongodb plugin

# Connect to database
client = MongoClient('localhost', 27017) # default port
db = client['chemnet']
chemicals = db['chemicals']
reactions = db['reactions']