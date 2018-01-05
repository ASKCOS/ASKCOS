from __future__ import print_function
import rdkit.Chem as Chem          
from rdkit.Chem import AllChem
import numpy as np
from collections import defaultdict

def average_template_list(INSTANCE_DB, CHEMICAL_DB, id_list):
	'''
	Given the INSTANCE_DB with many reaction examples and
	a list of IDs to include, this function will look at 
	the reaction conditions and come up with an estimate
	of proposed conditions
	'''

	solvents = defaultdict(int)
	reagents = defaultdict(int)
	catalysts = defaultdict(int)
	pressures = []
	temps = []
	yields = []
	times = []

	def string_or_range_to_float(text):
		try:
			return float(text)
		except Exception as e:
			if '-' in text:
				try:
					return sum([float(x) for x in text.split('-')]) / len(text.split('-'))
				except Exception as e:
					print(e)
			else:
				print(e)
		return None

	N_id = float(len(id_list))
	for _id in id_list:
		doc = INSTANCE_DB.find_one({'_id': _id})
		if not doc: continue 
		for xrn in doc['RXD_SOLXRN']: solvents[xrn] += 1
		for xrn in doc['RXD_CATXRN']: reagents[xrn] += 1
		for xrn in doc['RXD_RGTXRN']: catalysts[xrn] += 1
		if doc['RXD_P'] != -1: 
			P = string_or_range_to_float(doc['RXD_P'])
			if P: pressures.append(P)
		if doc['RXD_T'] != -1: 
			T = string_or_range_to_float(doc['RXD_T'])
			if T: temps.append(T)
		if doc['RXD_TIM'] != -1: 
			t = string_or_range_to_float(doc['RXD_TIM'])
			if t: times.append(t)
		if doc['RXD_NYD'] != -1: yields.append(float(doc['RXD_NYD']))
	conditions = []

	# Solvents
	solvent_string = ''
	for (solxrn, count) in sorted(solvents.items(), key = lambda x: x[1], reverse = True)[:5]:
		doc = CHEMICAL_DB.find_one({'_id': solxrn}, ['SMILES', 'IDE_CN'])
		cn = doc['IDE_CN'] if doc else solxrn
		solvent_string += '{} ({:.0f}%); '.format(cn, count*100.0/N_id)
	conditions.append('SOLVENT: ' + solvent_string)

	# Reagents
	reagent_string = ''
	for (rgtxrn, count) in sorted(reagents.items(), key = lambda x: x[1], reverse = True)[:5]:
		doc = CHEMICAL_DB.find_one({'_id': rgtxrn}, ['SMILES', 'IDE_CN'])
		cn = doc['IDE_CN'] if doc else rgtxrn
		reagent_string += '{} ({:.0f}%); '.format(cn, count*100.0/N_id)
	conditions.append('REAGENT: ' + reagent_string)

	# Catalysts
	catalyst_string = ''
	for (catxrn, count) in sorted(catalysts.items(), key = lambda x: x[1], reverse = True)[:5]:
		doc = CHEMICAL_DB.find_one({'_id': catxrn}, ['SMILES', 'IDE_CN'])
		cn = doc['IDE_CN'] if doc else catxrn
		catalyst_string += '{} ({:.0f}%); '.format(cn, count*100.0/N_id)
	conditions.append('CATALYST: ' + catalyst_string)

	# Time
	if times:
		conditions.append(
			'TIME: {:.1f} +/- {:.1f} hours (min {:.1f}, max {:.1f}, N={})'.format(
				np.mean(times), np.std(times), min(times), max(times), len(times)
			)
		)
	else:
		conditions.append('TIME UNKNOWN')

	# Temp
	if temps:
		conditions.append(
			'TEMP: {:.1f} +/- {:.1f} Celsius (min {:.1f}, max {:.1f}, N={})'.format(
				np.mean(temps), np.std(temps), min(temps), max(temps), len(temps)
			)
		)
	else:
		conditions.append('TEMP UNKNOWN')

	# Pressure
	if pressures:
		conditions.append(
			'PRESSURE: {:.0f} +/- {:.0f} Torr (min {:.0f}, max {:.0f}, N={})'.format(
				np.mean(pressures), np.std(pressures), min(pressures), max(pressures), len(pressures)
			)
		)
	else:
		conditions.append('PRESSURE UNKNOWN')

	# Yields
	if yields:
		conditions.append(
			'To provide a yield of: {:.1f} +/- {:.1f} percent (min {:.1f}, max {:.1f}, N={})'.format(
				np.mean(yields), np.std(yields), min(yields), max(yields), len(yields)
			)
		)
	else:
		conditions.append('YIELD UNKNOWN')

	return conditions