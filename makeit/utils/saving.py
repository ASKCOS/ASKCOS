from __future__ import print_function
import datetime

def save_model_history(hist, fpath):
	''''This function saves the history returned by model.fit to a tab-
	delimited file, where model is a keras model'''

	# Open file
	fid = open(fpath, 'a')
	print('trained at {}'.format(datetime.datetime.utcnow()))
	print('iteration\tnum_batches\tbatch_size\tloss\taccuracy', file = fid)

	# Iterate through
	for i in range(len(hist.history['batch'])):
		print('{}\t{}\t{}\t{}'.format(i + 1, 
			                          hist.history['batch'][i], 
			                          hist.history['size'][i],
			                          hist.history['loss'][i]), file = fid)

	# Close file
	fid.close()
	return True