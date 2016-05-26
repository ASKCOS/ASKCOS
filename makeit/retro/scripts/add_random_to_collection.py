from __future__ import print_function
from pymongo import MongoClient    # mongodb plugin
import numpy as np
import argparse

def main(collection, overwrite = False):
	'''
	Given a MongoDB collection, this function creates (or overwrites, 
	if specified) the 'random' field of each document to contain a 
	random float between 0 and 1.
	'''

	if overwrite: 
		find_filter = {}
		find_filter_string = '{}'
	else:
		find_filter = {'random': { '$exists': False }}
		find_filter_string = "{'random': { '$exists': false }}"

	cursor = collection.find(find_filter)
	print('Found {} matching entries'.format(cursor.count()))

	collection_name = str(collection._Collection__name)

	mongo_string = 'db.getCollection("{}").find({})'.format(collection_name, find_filter_string)
	mongo_string += '.forEach(function(mydoc) {}db.getCollection("{}")'.format('{', collection_name)
	mongo_string += '.update({}_id: mydoc._id{}, {}$set: {}random: Math.random() {}) {})'.format('{', '}', '{', '{', '}}', '}')
	
	print('SENDING COMMAND')
	print(mongo_string)
	print('')

	try:
		db = collection._Collection__database
		result = db.eval(mongo_string)
	except:
		print('Not authorized - copy and paste string into Mongo shell instead')

	
	
if __name__ == '__main__':

	parser = argparse.ArgumentParser()
	parser.add_argument('database', type = str,
						help = 'Database hosting collection')
	parser.add_argument('collection', type = str,
						help = 'Collection to add random field to')
	parser.add_argument('-o', '--overwrite', type = bool, default = False,
						help = 'Overwrite existing random field (if exists)')
	parser.add_argument('-v', type = bool, default = False,
						help = 'Verbose printing (incl. saving images); defaults to False')
	args = parser.parse_args()


	client = MongoClient('mongodb://guest:guest@rmg.mit.edu/admin', 27017)
	db = client[args.database]
	example_collection = db[args.collection]

	main(example_collection, overwrite = bool(args.overwrite))