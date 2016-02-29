## Generate input files for CV
# Assumes '{}' are in order of (# folds, fold #, fold #, # folds, seed)
import numpy.random as rand
import sys
import os

if __name__ == '__main__':
	# Read number of folds from command line
	if len(sys.argv) < 2:
		N_folds = 5
	else:
		N_folds = int(sys.argv[1])

	## Get random seed to use for all folds
	seed = rand.randint(100000)

	## Read/write
	run_string = ''
	run_root = 'python makeit/main/neural_fp_opensol.py '
	directory = os.path.dirname(os.path.realpath(__file__))
	print(directory)
	with open(os.path.join(directory, 'input.cfg'), 'r') as template_fid:
		template = template_fid.read()

		# Write each fold file
		for i in range(N_folds):
			fold_num = i + 1
			with open(os.path.join(directory, 'fold{}.cfg'.format(fold_num)), 'w') as fold_fid:
				fold_fid.write(template.format(N_folds, fold_num, fold_num, N_folds, seed))
				run_string += run_root + '"' + os.path.join(directory, 'fold{}.cfg'.format(fold_num)) + '"\n'

	with open(os.path.join(directory, 'run.sh'), 'w') as run_fid:
		run_fid.write(run_string)