class Collection_Manager:
	'''Collection_Manager is an auxilliary class to help manage large
	PyMongo collection objects. It returns a small number of results
	for the current collection based on the specified batch_size.'''

	def __init__(self, collection, batch_size = 10):
		'''Create a new Collection_Manager object'''

		# PyMongo collection object
		self.collection = collection # PyMongo collection

		# Batch size (number of results to return)
		self.batch_size = batch_size
		
		# Current index
		self.index = 0

		# Total number of documents
		self.num_documents = self.collection.find().count()

		# Have reported all records?
		self.done = False

	def pull(self):
		'''Pull the next set of documents'''

		# Get index numbers
		start = self.index
		end   = self.index + self.batch_size
		self.index = end # increment

		# Check dimension
		if end < self.num_documents:
			return self.collection.find()[start:end]
		else:
			self.done = True
			self.index = self.num_documents
			return self.collection.find()[start:]


def get_all_ids(collection, db_filter = {}):
	'''Finds all of the '_id's in the collection (chemicals or reactions) to 
	allow easier randomization of database subsets. Allows for filters to be
	used as well (db_filter == dictionary, first argument of find()).'''

	matched_docs = collection.find(db_filter, {"_id": 1})
	return [x['_id'] for x in matched_docs]
