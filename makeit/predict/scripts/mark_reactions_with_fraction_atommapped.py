# Import relevant packages
from __future__ import print_function
from global_config import USE_STEREOCHEMISTRY
import argparse
import numpy as np     	      	   # for simple calculations
import os                          # for saving
import sys
import rdkit.Chem as Chem
import rdkit.Chem.AllChem as AllChem
from pymongo import MongoClient    # mongodb plugin
import time
from tqdm import tqdm
from multiprocessing import Process, Lock, Queue
from Queue import Empty as QueueEmpty

'''
This script goes through the reactions we have scraped from Reaxys and 
adds additional fields for the FRAC_PROD_MAPPED so we know that some
reactions, even with RX_SKW: mapped reaction, aren't actually atom
mapped to a satisfactory degree.
'''


def process_one(rx, REACTION_DB, lock_pymongo):

	rxn_smiles = rx['RXN_SMILES']
	products = rxn_smiles.split('>')[2]

	products = Chem.MolFromSmiles(products)

	if not products:
		print('Could not load products for ID {}'.format(rx['_id']))
		return
	AllChem.RemoveHs(products)

	mapped_atoms = sum([a.HasProp('molAtomMapNumber') for a in products.GetAtoms()])
	total_atoms = products.GetNumAtoms()
	frac = mapped_atoms / float(total_atoms)

	#lock_pymongo.acquire()
	REACTION_DB.update_one(
		{'_id': rx['_id']},
		{'$set': 
			{'FRAC_PROD_MAPPED': frac,
			 'PROD_NUMATOMS': total_atoms},
		}
	)
	#lock_pymongo.release()

	#print('[RX {}] Processed, {} mapped'.format(rx['_id'], frac))

	return

def process_forever(queue, lock_pymongo):
	client = MongoClient('mongodb://guest:guest@rmg.mit.edu/admin', 27017)
	db = client['reaxys']
	REACTION_DB = db['reactions']
	while True:
		try:
			rx = queue.get()
			process_one(rx, REACTION_DB, lock_pymongo)
		except QueueEmpty:
			time.sleep(1)

def rx_generator():
	'''Return (rx, rxd) tuples that have not been processed'''

	for rx in REACTION_DB.find({'RX_SKW': 'mapped reaction', 'FRAC_PROD_MAPPED': {'$exists': False}}, ['_id', 'RXN_SMILES'], 
			no_cursor_timeout = True).sort('_id', 1):
		if not rx: continue 
		yield rx

if __name__ == '__main__':

	parser = argparse.ArgumentParser()
	parser.add_argument('-n', '--num', type = int, default = 50,
						help = 'Maximum number of records to examine; defaults to 50')
	parser.add_argument('--workers', type = int, default = 10,
						help = 'Number of parallel workers, default 10')
	args = parser.parse_args()

	n_max = int(args.num), 
	complete_only = True

	from rdkit import RDLogger
	lg = RDLogger.logger()
	lg.setLevel(4)

	client = MongoClient('mongodb://guest:guest@rmg.mit.edu/admin', 27017)
	db = client['reaxys']
	REACTION_DB = db['reactions']

	# Run
	try:
		queue = Queue() 
		processes = []
		lock = Lock()
		lock_pymongo = Lock()
		for i in range(int(args.workers)):
			processes.append(
				Process(target=process_forever, args=(queue, lock_pymongo))
			)
		[p.start() for p in processes]
		print('Created {} processes'.format(len(processes)))

		# Keep queue loaded
		generator = rx_generator()
		j = 0
		while True:
			if j % 5000 == 0:
				print('Total reactions queued up: {}'.format(j))
			if j == n_max: break
			if queue.qsize() < 5000:
				queue.put(generator.next())
				j += 1
			else:
				time.sleep(5)
				for i in range(len(processes)):
					process = processes[i]
					if not process.is_alive():
						new_process = Process(target=process_forever, args=(queue, lock_pymongo))
						new_process.start()
						processes[i] = new_process 
						del process
						print('### Restarted process {} ###'.format(i))

				continue # wait

		# Wait for queue to empty before calling .join()
		while not queue.empty():
			print('Stopped adding to queue, just waiting to empty it out...')
			time.sleep(3)

		[p.join() for p in processes]

	except KeyboardInterrupt:
		print('Stopped early, leaving pool')
