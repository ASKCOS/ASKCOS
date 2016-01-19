import re # regular expression processing
import rdkit.Chem as Chem
import numpy as np

def SplitChemicalName(name):
	'''This function takes a raw chemical name and replaces any common
		chemical delimeters with spaces for easier processing'''

	# (1) Split by delimeters: ; . , \W -
	# (2) Strip each string in list
	# (3) Remove empty strings from list
	return  filter(None, [str(token).strip() for token in re.split(r'[\(\);.,\W-]', name)])

def input_to_bool(txt):
	'''This function converts a raw_input string to a booolean'''
	return txt.lower() in ['yes', 'true', 'y', 't', '1']

def sequences_to_texts(tokenizer, sequences):
	'''This function takes list of sequences (i.e., indeces) and converts each
	back to the original sentence using the tokenizer's indexing. The 
	sequences can be padded or not.'''

	# Convert tokenizer word_index to sorted list
	ordered_vocab = sorted(tokenizer.__dict__['word_index'], 
		                   key = tokenizer.__dict__['word_index'].get)

	# Custom evaluation
	texts = []
	for i in range(len(sequences)):
		words = []
		for index in sequences[i]:
			if index == 0:
				continue
			words.append(ordered_vocab[index - 1])
		texts.append(' '.join(words))
	return texts

def smiles_to_boolean_fps(all_smiles):
	'''This function takes a list of SMILES strings and returns a fingerprint 
	in the form of an explicit boolean vector'''

	BVs = [Chem.RDKFingerprint(Chem.MolFromSmiles(smiles)) for smiles in all_smiles]
	for i in range(len(BVs)):
		BVs[i] = np.array([bool(index) for index in BVs[i]])
	return np.vstack(BVs)
