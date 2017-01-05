# Import relevant packages
from makeit.utils.chemnet_connect import *      # mongodb connection, gets 'chemicals', 'reactions'
import numpy as np     	      	   # for simple calculations
import matplotlib.pyplot as plt    # for visualization
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})
import os                          # for saving
from tqdm import tqdm
import cPickle as pickle


# DATABASE
from pymongo import MongoClient    # mongodb plugin
client = MongoClient('mongodb://guest:guest@rmg.mit.edu/admin', 27017)
db = client['reaxys']
INSTANCE_DB = db['instances']


if os.path.isfile(os.path.join(os.path.join(os.path.dirname(__file__), 'output'), 'yields_complete.pickle')):
	with open(os.path.join(os.path.join(os.path.dirname(__file__), 'output'), 'yields_complete.pickle'), 'rb') as fid:
		yields = pickle.load(fid)
else:
	# Get distribution of (median) yields for all reactions
	yields = np.zeros((INSTANCE_DB.count({'complete': True})))
	for i, doc in tqdm(enumerate(INSTANCE_DB.find({'complete': True}, ['RXD_NYD']).limit(len(yields)))):
		yields[i] = float(doc['RXD_NYD'])

yields = yields[yields > 0.0]

# Visualize yields in histogram
weights = np.ones_like(yields) / len(yields)
fig = plt.figure(figsize=(6,4), dpi = 300)
n, bins, patches = plt.hist(yields, range(101), facecolor = 'blue', alpha = 1, weights = weights)
plt.xlabel('Reaction yield [%]')
plt.ylabel('Normalized frequency [-]')
plt.title('Yields of {} Reaxys reaction examples'.format(len(yields)))
ymin, ymax = plt.ylim()
plt.axis([0, 100, 0, ymax])
plt.grid(True)
plt.savefig(os.path.join(os.path.join(os.path.dirname(__file__), 'output'), 'yields_complete.png'), bbox_inches = 'tight')
with open(os.path.join(os.path.join(os.path.dirname(__file__), 'output'), 'yields_complete.pickle'), 'wb') as fid:
	pickle.dump(yields, fid, pickle.HIGHEST_PROTOCOL)
