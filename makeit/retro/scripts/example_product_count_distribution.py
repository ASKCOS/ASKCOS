# Import relevant packages
import numpy as np     	      	   # for simple calculations
import matplotlib
import matplotlib.pyplot as plt    # for visualization
import os                          # for saving
matplotlib.rc('font', **{'size': 18})

# DATABASE
from pymongo import MongoClient    # mongodb plugin
client = MongoClient('mongodb://guest:guest@rmg.mit.edu/admin', 27017)
db = client['reaction_examples']
reactions = db['lowe_1976-2013_USPTOgrants_reactions']

# Get distribution of (median) yields for all reactions
N = reactions.count()
all_prodcounts = np.nan * np.ones((N, 1))
for i, reaction in enumerate(reactions.find()):
	all_prodcounts[i] = len(reaction['products'])

# Filter out NaN values (which we get if we take np.median of empty array)
print 'filtering ' + str(len(all_prodcounts[np.isnan(all_prodcounts)])) + ' entries with num prods == NaN'
all_prodcounts = all_prodcounts[~np.isnan(all_prodcounts)]
print str(len(all_prodcounts)) + ' entries remaining in final dataset'

# Visualize yields in histogram
weights = np.ones_like(all_prodcounts) / len(all_prodcounts)
n, bins, patches = plt.hist(all_prodcounts, 50, facecolor = 'blue', alpha = 0.5, weights = weights)
plt.xlabel('Number of products')
plt.ylabel('Normalized frequency')
plt.title('Histogram of product counts')
ymin, ymax = plt.ylim()
xmin, xmax = plt.xlim()
plt.axis([1, 4, 0, ymax])
plt.grid(True)
plt.savefig(os.path.join(os.path.join(os.path.dirname(__file__), 'output'), \
	'reaction_prod_counts.png'), bbox_inches = 'tight')

# Number with one
print('Number with only one: {}'.format(sum(all_prodcounts == 1)))