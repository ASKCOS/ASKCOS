# parse_cfg is used to read settings in a configuration file and check the validitiy of those
import ConfigParser
import sys
import os


def read_config(fpath):
	'''This function reads a configuration file and returns an equivalent dictionary'''

	# Create new parser
	config = ConfigParser.ConfigParser()
	with open(fpath, 'r') as fid:
		config.readfp(fid)

	# TODO: write defaults and check data

	return config._sections

if __name__ == '__main__':
	if len(sys.argv) < 2:
		quit('usage: {} "file.cfg"'.format(sys.argv[0]))

	print(read_config(sys.argv[1]))
