import rdkit.Chem as Chem

def summarize_reaction_outcome(mols, outcome):

	h_lost = []
	h_gain = []
	bond_lost = []
	bond_gain = []
	
	conserved_maps = [a.GetProp('molAtomMapNumber') for a in outcome.GetAtoms() if a.HasProp('molAtomMapNumber')]
	changes = 0

	for atom_prev in mols.GetAtoms():
		atomMapNumber = atom_prev.GetProp('molAtomMapNumber')
		atom_new = [a for a in outcome.GetAtoms() if a.HasProp('molAtomMapNumber') and a.GetProp('molAtomMapNumber') == atomMapNumber]
		if not atom_new: continue
		atom_new = atom_new[0]
		
		Hs_prev = atom_prev.GetTotalNumHs()
		Hs_new  = atom_new.GetTotalNumHs()
		if Hs_prev < Hs_new:
			#print('    atom {} gained {} hydrogens'.format(atomMapNumber, Hs_new - Hs_prev))
			for i in range(Hs_prev, Hs_new):
				h_gain.append(atomMapNumber)
				changes += 1
		if Hs_prev > Hs_new:
			#print('    atom {} lost {} hydrogens'.format(atomMapNumber, Hs_prev - Hs_new))
			for i in range(Hs_new, Hs_prev): 
				h_lost.append(atomMapNumber)
				changes += 1

		# Charge_prev = atom_prev.GetFormalCharge()
		# Charge_new = atom_new.GetFormalCharge()
		# if Charge_prev != Charge_new:
		# 	#print('    atom {} changed charge ({} to {})'.format(atomMapNumber, Charge_prev, Charge_new))
		# 	changes += 1

	bonds_prev = {}
	for bond in mols.GetBonds():
		nums = sorted(
			[bond.GetBeginAtom().GetProp('molAtomMapNumber'),
			bond.GetEndAtom().GetProp('molAtomMapNumber')])
		if (nums[0] not in conserved_maps) and (nums[1] not in conserved_maps): continue
		bonds_prev['{}~{}'.format(nums[0], nums[1])] = bond.GetBondTypeAsDouble()
	bonds_new = {}
	for bond in outcome.GetBonds():
		nums = sorted(
			[bond.GetBeginAtom().GetProp('molAtomMapNumber'),
			bond.GetEndAtom().GetProp('molAtomMapNumber')])
		bonds_new['{}~{}'.format(nums[0], nums[1])] = bond.GetBondTypeAsDouble()
	
	# print('prev: {}'.format(Chem.MolToSmarts(mols)))
	# print('new: {}'.format(Chem.MolToSmarts(outcome)))
	# print('bonds_prev: {}'.format(bonds_prev))
	# print('bonds_new: {}'.format(bonds_new))

	for bond in bonds_prev:
		if bond not in bonds_new:
			#print('    lost bond {}, order {}'.format(bond, bonds_prev[bond]))
			bond_lost.append((bond.split('~')[0], bond.split('~')[1], bonds_prev[bond]))
			changes += 1
		else:
			if bonds_prev[bond] != bonds_new[bond]:
				#print('    changed bond {} from order {} to {}'.format(bond, bonds_prev[bond], bonds_new[bond]))
				bond_lost.append((bond.split('~')[0], bond.split('~')[1], bonds_prev[bond]))
				bond_gain.append((bond.split('~')[0], bond.split('~')[1], bonds_new[bond]))
				changes += 1
	for bond in bonds_new:
		if bond not in bonds_prev:
			#print('    new bond {}, order {}'.format(bond, bonds_new[bond]))
			bond_gain.append((bond.split('~')[0], bond.split('~')[1], bonds_new[bond]))
			changes += 1

	#print('{} total changes'.format(changes))
	return (sorted(h_lost), sorted(h_gain), sorted(bond_lost), sorted(bond_gain))