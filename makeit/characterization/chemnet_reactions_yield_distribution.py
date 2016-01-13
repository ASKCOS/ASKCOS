# Import relevant packages
from makeit.utils.chemnet_connect import *      # mongodb connection, gets 'chemicals', 'reactions'
import numpy as np     	      	   # for simple calculations
import matplotlib.pyplot as plt    # for visualization
import os                          # for saving

# Get distribution of (median) yields for all reactions
all_yields = np.nan * np.ones(1000000)
for i, reaction in enumerate(reactions.find()):

	# Get rxid tag
	rxid = reaction['rxid'][0]
	
	# Try reading reference list
	try: 
		references_tag = 'rx' + str(rxid)
		references = reaction[references_tag]

		# Get list of yields from all references
		this_yields = np.zeros(0)
		for reference in references:
			try:
				this_yields = np.append(this_yields, reference['yield'][0])
			except: # yield not reported
				pass 

		# Get median yield
		if len(this_yields) != 0:
			all_yields[i] = np.median(this_yields)

	# No reference list provided, only one yield entry
	except: 
		all_yields[i] = reaction['yield']

	# # Terminate early (for testing)
	# if i > 1000: 
	# 	all_yields = all_yields[:1000]
	# 	break

# Filter out NaN values (which we get if we take np.median of empty array)
print 'filtering ' + str(len(all_yields[np.isnan(all_yields)])) + ' entries with yield == NaN'
all_yields = all_yields[~np.isnan(all_yields)]
# Filter out > 1 values (which we get from poor data quality)
print 'filtering ' + str(len(all_yields[all_yields > 1])) + ' entries with yield > 1'
all_yields = all_yields[all_yields <= 1]
print str(len(all_yields)) + ' entries remaining in final dataset'
# Determine how many are just reporting 0.75
print str(len(all_yields[all_yields == 0.75])) + ' entries with median = 0.75'

# Visualize yields in histogram
weights = np.ones_like(all_yields) / len(all_yields)
n, bins, patches = plt.hist(all_yields, 50, facecolor = 'blue', alpha = 0.5, weights = weights)
plt.xlabel('Median reaction yield from all sources')
plt.ylabel('Normalized frequency')
plt.title('Histogram of yields in Bishop\'s 1M reaction dataset')
ymin, ymax = plt.ylim()
plt.axis([0, 1, 0, ymax])
plt.grid(True)
plt.savefig(os.path.dirname(__file__) + '/Histogram of reaction yields.png', bbox_inches = 'tight')

# Filter out yields equal to 0.75
all_yields = all_yields[all_yields != 0.75]

# Visualize yields in histogram
plt.clf()
weights = np.ones_like(all_yields) / len(all_yields)
n, bins, patches = plt.hist(all_yields, 50, facecolor = 'blue', alpha = 0.5, weights = weights)
plt.xlabel('Median reaction yield from all sources')
plt.ylabel('Normalized frequency')
plt.title('Histogram of yields in Bishop\'s 1M reaction dataset')
ymin, ymax = plt.ylim()
plt.axis([0, 1, 0, ymax])
plt.grid(True)
plt.savefig(os.path.dirname(__file__) + '/Histogram of reaction yields filter_0d75.png', bbox_inches = 'tight')