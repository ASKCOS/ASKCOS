# Import relevant packages
from makeit.utils.chemnet_connect import *      # mongodb connection, gets 'chemicals', 'reactions'
import numpy as np     	      	   # for simple calculations
import matplotlib.pyplot as plt    # for visualization
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})
import os                          # for saving
from tqdm import tqdm
import cPickle as pickle


def string_or_range_to_float(text):
	try:
		return float(text)
	except Exception as e:
		x = [z for z in text.strip().split('-') if z not in [u'', u' ']]
		if text.count('-') == 1: # 20 - 30
			try:
				return (float(x[0]) + float(x[1])) / 2.0
			except Exception as e:
				print('Could not convert {}, {}'.format(text, x))
				#print(e)
		elif text.count('-') == 2: # -20 - 0
			try:
				return (-float(x[0]) + float(x[1])) / 2.0
			except Exception as e:
				print('Could not convert {}, {}'.format(text, x))
				#print(e)
		elif text.count('-') == 3: # -20 - -10
			try:
				return (-float(x[0]) - float(x[1])) / 2.0
			except Exception as e:
				print('Could not convert {}, {}'.format(text, x))
				#print(e)
		else:
			print('Could not convert {}'.format(text))
			print(e)
	return np.nan


# DATABASE
from pymongo import MongoClient    # mongodb plugin
client = MongoClient('mongodb://guest:guest@rmg.mit.edu/admin', 27017)
db = client['reaxys']
INSTANCE_DB = db['instances']

if os.path.isfile(os.path.join(os.path.join(os.path.dirname(__file__), 'output'), 'temps_complete.pickle')):
	with open(os.path.join(os.path.join(os.path.dirname(__file__), 'output'), 'temps_complete.pickle'), 'rb') as fid:
		temps = pickle.load(fid)
else:
	# Get distribution of (median) temps for all reactions
	temps = np.nan * np.ones((INSTANCE_DB.count({'complete': True})))
	for i, doc in tqdm(enumerate(INSTANCE_DB.find({'complete': True}, ['RXD_T']).limit(len(temps)))):
		try:
			temps[i] = string_or_range_to_float(doc['RXD_T'])
		except Exception as e:
			print(e)

temps = temps[~np.isnan(temps)]

# Visualize temps in histogram
weights = np.ones_like(temps) / len(temps)
fig = plt.figure(figsize=(6,4), dpi = 300)
n, bins, patches = plt.hist(temps, range(-200, 400, 5), facecolor = 'blue', alpha = 1, weights = weights)
plt.xlabel('Temperature [C]')
plt.ylabel('Normalized frequency [-]')
plt.title('Temps of {} Reaxys reaction examples'.format(len(temps)))
ymin, ymax = plt.ylim()
# plt.axis([0, 100, 0, ymax])
plt.grid(True)
plt.savefig(os.path.join(os.path.join(os.path.dirname(__file__), 'output'), 'temps_complete.png'), bbox_inches = 'tight')
with open(os.path.join(os.path.join(os.path.dirname(__file__), 'output'), 'temps_complete.pickle'), 'wb') as fid:
	pickle.dump(temps, fid, pickle.HIGHEST_PROTOCOL)
