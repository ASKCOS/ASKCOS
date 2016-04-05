from makeit.utils.draw import *
import argparse
import shutil
import os

def main(tform_file, folder, N = 0, ext = 'jpg'):
	'''This function generates a webpage displaying all of the
	retrosyntheses in tform_file in a semi-nice format.'''

	# Prepare folders
	if not os.path.exists(folder):
		os.makedirs(folder)
	media_folder = os.path.join(folder, 'media')
	if not os.path.exists(media_folder):
		os.makedirs(media_folder)

	# Copy transform file into folder
	tform_file_copy = os.path.join(folder, 'source_transforms.txt')
	shutil.copy(tform_file, tform_file_copy)

	with open(os.path.join(folder, 'index.html'), 'w') as out_fid:
		# Begin webpage
		out_fid.write('<html>\n')
		out_fid.write('\t<body>\n')
		out_fid.write('\t\t<h1>List of extracted retrosynthetic transforms</h1>\n')
		out_fid.write('\t\t<h3><a href="source_transforms.txt">(see source)</a><h3>\n')
		out_fid.write('\t\t<table border="1" align="center">\n')

		# Headers
		out_fid.write('\t\t\t<tr>\n')
		out_fid.write('\t\t\t\t<th align="center">Index</th>\n')
		out_fid.write('\t\t\t\t<th align="center">\n')
		out_fid.write('\t\t\t\t\t<font face="Courier New, Courier monospace">Transform</font><br>\n')
		out_fid.write('\t\t\t\t</th>\n')
		out_fid.write('\t\t\t\t<th align="center">Score</th>\n')
		out_fid.write('\t\t\t</tr>\n')

		# Add all 
		with open(tform_file, 'r') as in_fid:
			for i, line in enumerate(in_fid):
				if N and i == N: break

				# Get image
				j, tform, score = line.strip().split('\t')
				img = TransformStringToImage(tform)
				if ext != 'png':
					img = StripAlphaFromImage(img)
				img.save(os.path.join(media_folder, 'rxn{}_score{}.{}'.format(j, score, ext)))

				# Add to table
				out_fid.write('\t\t\t<tr align="center">\n')
				out_fid.write('\t\t\t\t<td align="center">{}</td>\n'.format(j))
				out_fid.write('\t\t\t\t<td align="center">\n')
				out_fid.write('\t\t\t\t\t<font face="Courier New, Courier monospace">{}</font><br>\n'.format(tform))
				out_fid.write('\t\t\t\t\t<img src="media/rxn{}_score{}.{}">\n'.format(j, score, ext))
				out_fid.write('\t\t\t\t</td>\n')
				out_fid.write('\t\t\t\t<td align="center">{}</td>\n'.format(score))
				out_fid.write('\t\t\t</tr>\n')

		# Finish webpage
		out_fid.write('\t\t</table>\n')
		out_fid.write('\t</body>\n')
		out_fid.write('</html>\n')

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('tform_file', type = str, 
		 				help = 'Tab-delimited file where each line contains an index, '
		 				'SMARTS retrosynthetic transform, and score.')
	parser.add_argument('folder', type = str,
						help = 'Folder name to generate webpage inside.')
	parser.add_argument('-n', '--num', type = int, default = 0,
						help = 'Maximum number of records to examine; default is unlimited')
	parser.add_argument('-e', '--ext', type = str, default = 'jpg',
						help = 'Extension of image files; default is jpg')
	args = parser.parse_args()

	main(args.tform_file, args.folder, N = args.num, ext = args.ext)