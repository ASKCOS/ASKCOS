## Generate input files for CV
import numpy.random as rand
import sys
import os

if __name__ == '__main__':
	# Read number of folds from command line
	if len(sys.argv) < 2:
		N_folds = 5
	else:
		N_folds = int(sys.argv[1])

	# List of hidden nodes to try
	hiddens = [50]
	lengths = [256]
	decays = [40]
	depths = [3, 4]
	solvents = ['methanol']

	# Get random tag for all files
	tag = rand.randint(1000)

	## Read/write
	run_string = ''
	run_root = 'python makeit/main/neural_fp_opensol.py '
	directory = os.path.dirname(os.path.realpath(__file__))
	print(directory)
	with open(os.path.join(directory, 'input.cfg'), 'r') as template_fid:
		template = template_fid.read()

		for solvent in solvents:
			num = 0
			for depth in depths:
				for decay in decays:
					for length in lengths:
						for hidden in hiddens:
							num = num + 1
							# Model directory
							model_directory = 'makeit/models/opensol_<solvent>_<tag>_v<num>'
							model_directory = model_directory.replace('<solvent>', solvent)
							model_directory = model_directory.replace('<tag>', str(tag))
							model_directory = model_directory.replace('<num>', str(num))
							try:
								os.makedirs(model_directory)
							except: # file exists
								pass

							# Get random seed to use for all folds
							seed = rand.randint(100000)

							# Customize template for these conditions
							this_cfg = template.replace('<solvent>', solvent)
							this_cfg = this_cfg.replace('<tag>', str(tag))
							this_cfg = this_cfg.replace('<num>', str(num))
							this_cfg = this_cfg.replace('<depth>', str(depth))
							this_cfg = this_cfg.replace('<decay>', str(float(decay)))
							this_cfg = this_cfg.replace('<length>', str(length))
							this_cfg = this_cfg.replace('<hidden>', str(hidden))
							this_cfg = this_cfg.replace('<seed>', str(seed))
							this_cfg = this_cfg.replace('<num_folds>', str(N_folds))

							flabel = 'solvent{}_tag{}_num{}_depth{}_decay{}_length{}_hidden{}'.format(solvent, tag, num, depth, decay, length, hidden)

							# Write each fold file
							for i in range(N_folds):
								fold_num = i + 1
								with open(os.path.join(directory, flabel + '_fold{}.cfg'.format(fold_num)), 'w') as fold_fid:
									fold_fid.write(this_cfg.replace('<this_fold>', str(fold_num)))
									run_string += run_root + '"' + os.path.join(directory, flabel + '_fold{}.cfg'.format(fold_num)) + '"\n'

	with open(os.path.join(directory, 'run.sh'), 'w') as run_fid:
		run_fid.write(run_string)