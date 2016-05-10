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
N = reactions.count({'yield':{'$exists':True}})
all_yields = np.nan * np.ones((N, 1))
for i, reaction in enumerate(reactions.find({'yield':{'$exists':True}})):
	all_yields[i] = reaction['yield']

# Filter out NaN values (which we get if we take np.median of empty array)
print 'filtering ' + str(len(all_yields[np.isnan(all_yields)])) + ' entries with yield == NaN'
all_yields = all_yields[~np.isnan(all_yields)]
# Filter out > 1 values (which we get from poor data quality)
print 'filtering ' + str(len(all_yields[all_yields > 100])) + ' entries with yield > 100'
all_yields = all_yields[all_yields <= 100]
print str(len(all_yields)) + ' entries remaining in final dataset'

# Visualize yields in histogram
weights = np.ones_like(all_yields) / len(all_yields)
n, bins, patches = plt.hist(all_yields, 50, facecolor = 'blue', alpha = 0.5, weights = weights)
plt.xlabel('Recorded yield (%)')
plt.ylabel('Normalized frequency')
plt.title('Histogram of yields')
ymin, ymax = plt.ylim()
xmin, xmax = plt.xlim()
plt.axis([0, 100, 0, ymax])
plt.grid(True)
plt.savefig(os.path.join(os.path.join(os.path.dirname(__file__), 'output'), \
	'reaction_yields.png'), bbox_inches = 'tight')
