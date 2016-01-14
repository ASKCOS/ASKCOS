import re # regular expression processing

def SplitChemicalName(name):
	'''This function takes a raw chemical name and replaces any common
		chemical delimeters with spaces for easier processing'''

	# (1) Split by delimeters: ; . , \W -
	# (2) Strip each string in list
	# (3) Remove empty strings from list
	return  filter(None, [str(token).strip() for token in re.split(r'[\(\);.,\W-]', name)])