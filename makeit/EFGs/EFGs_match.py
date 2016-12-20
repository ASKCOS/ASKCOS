import rdkit.Chem as Chem
import rdkit.Chem.AllChem as AllChem

def get_EFGS_matches(mol, library = []):
	'''Given a molecule and a library (list) of EFG documents, 
	this function returns a list of which ones match as a 
	boolean vector'''

	if not library:
		from pymongo import MongoClient
		client = MongoClient('mongodb://guest:guest@rmg.mit.edu/admin', 27017)
		db = client['askcos_transforms']
		EFG_DB = db['EFGs']
		library = [doc for doc in EFG_DB.find()]

	return [EFG_match(mol, doc['SMARTS']) for doc in library]

def EFG_match(mol, smarts):
	'''Given a molecule and an EFG smarts, this function 
	checks to see if there is a match and returns a bool'''

	match = True 
	for frag in str(smarts).split('AND'):
		negate = False
		if 'NOT' in frag:
			frag = frag.split('NOT')[1]
			negate = True
		frag = frag.strip()
		match *= mol.HasSubstructMatch(AllChem.MolFromSmarts(frag)) != negate
	return bool(match)