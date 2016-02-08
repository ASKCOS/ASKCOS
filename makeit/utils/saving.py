from __future__ import print_function
import datetime

def save_model_history(hist, fpath):
	''''This function saves the history returned by model.fit to a tab-
	delimited file, where model is a keras model'''

	# Open file
	fid = open(fpath, 'a')
	print('trained at {}'.format(datetime.datetime.utcnow()))
	print('iteration\tnum_batches\tbatch_size\tloss\tval_loss', file = fid)

	try:
		# Iterate through
		for i in range(len(hist.history['batch'])):
			print('{}\t{}\t{}\t{}\t{}'.format(i + 1, 
							hist.history['batch'][i],
							hist.history['size'][i],
							hist.history['loss'][i], 
							hist.history['val_loss'][i] if 'val_loss' in hist.history else 'n/a'),
							file = fid)
	except KeyError:
		print('<no history found>', file = fid)

	# Close file
	fid.close()
	return True

def save_model_history_manual(loss, val_loss, fpath):
	''''This function saves the history returned by model.fit to a tab-
	delimited file, where model is a keras model'''

	# Open file
	fid = open(fpath, 'a')
	print('trained at {}'.format(datetime.datetime.utcnow()))
	print('iteration\tloss\tval_loss', file = fid)

	try:
		# Iterate through
		for i in range(len(loss)):
			print('{}\t{}\t{}'.format(i + 1, 
							loss[i], val_loss[i]),
							file = fid)
	except KeyError:
		print('<no history found>', file = fid)

	# Close file
	fid.close()
	return True
