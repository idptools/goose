'''
User facing functionality for analysis of generated proteins.

Limited for now, planning on expanding on in the future
'''

from goose.backend.protein import Protein
from sparrow import Protein as pr

from goose.backend.predictors.predict_nes import predict_nes_seq as _predict_nes_seq
from goose.backend.predictors.predict_mito import predict_mitochondrial_targeting as _predict_mitochondrial_targeting
from goose.backend.predictors.predict_nls import predict_nls_seq as _predict_nls_seq
from goose.backend.predictors.predict_phosphosite import predict_phosphorylation as _predict_phosphorylation
from goose.backend.predictors.predict_tad import predict_tad_seq as _predict_tad_seq

def properties(sequence, fractions=True):
	'''
	analyzes various sequence properties and returns them to the user

	parameters
	-----------
	sequence : string
		the amino acid sequence as a string
	fractions : bool
		whether or not to include the specific fractions 
		of all amino acids for the sequence

	returns
	-------
	props_dict : dict
		returns a dictionary of the properties. 
	'''
	props_dict = Protein.calc_basic_properties(sequence)
	props_dict['kappa'] = pr(sequence).kappa
	# if fractions are wanted, add to dict
	if fractions == True:
		props_dict['fractions'] = Protein.calc_frac(sequence)
	# return the dict
	return props_dict

def phosphosites(sequence, raw_vals = False):
	'''
	for getting analysis of potential phosphosites
	'''
	seq_phopsho = _predict_phosphorylation(sequence)
	if raw_vals == True:
		return seq_phospho

	Yphospho = seq_phopsho['Y']
	Tphospho = seq_phopsho['T']
	Sphospho = seq_phopsho['S']

	potential_S = []
	potential_Y = []
	potential_T = []
	for i in range(0, len(sequence)):
		curY = Yphospho[i]
		curT = Tphospho[i]
		curS = Sphospho[i]
		cur_aa = sequence[i]
		if i ==0:
			curslice = sequence[i:i+4]
		elif i == 1:
			curslice = sequence[0:i+4]
		elif i == 2:
			curslice = sequence[0:i+4]
		elif i == 3:
			curslice = sequence[0:i+4]
		elif i >= 4:
			if i < len(sequence)+5:
				curslice = sequence[i-4:i+4]
			else: 
				curslice = sequence[i-4:]
		else:
			if i < len(sequence)+5:
				curslice = sequence[i:i+4]
			else:
				curslice = sequence[:i]

		if curY >= 0.6:
			if 'Y' in curslice:
				ylocation = sequence.index(curslice)+curslice.index('Y')
				if sequence[ylocation] == 'Y':
					if ylocation+1 not in potential_Y:
						potential_Y.append(ylocation+1)

		if curT >= 0.6:
			if 'T' in curslice:
				tlocation = sequence.index(curslice)+curslice.index('T')
				if sequence[tlocation] == 'T':
					if tlocation+1 not in potential_T:
						potential_T.append(tlocation+1)
		if curS >= 0.6:
			if 'S' in curslice:
				slocation = sequence.index(curslice)+curslice.index('S')
				if sequence[slocation] == 'S':
					if slocation+1 not in potential_S:
						potential_S.append(slocation+1)

		phosphosite_dict = {}
		if 	potential_S != []:
			phosphosite_dict['S']=potential_S
		if 	potential_Y != []:
			phosphosite_dict['Y']=potential_Y
		if 	potential_T != []:
			phosphosite_dict['T']=potential_T	
		if phosphosite_dict == {}:
			phosphosite_dict['S']='No S phosphorylation predicted'
			phosphosite_dict['Y']='No Y phosphorylation predicted'
			phosphosite_dict['T']='No T phosphorylation predicted'

	return phosphosite_dict				

def cellular_localization(sequence):
	'''
	returns whether the sequence may have any 
	cellular localization signals

	parameters
	---------
	sequence : string
		amino acid sequence as a string

	returns 
	-------
	all_locs : dict
		a dictionary of all predicted mitochondrial, NLS, 
		and NES sequences. The dict has three keys
		1. mitochondria
		2. NES
		3. NLS
		within each of those key pair values is another dictionary that has the 
		sequence of the predicted feature as well as the sequence coordinates.
	'''

	all_locs = {}
	mito_loc = {}
	nes_loc = {}
	nls_loc = {}

	mitochondrial_targeting = _predict_mitochondrial_targeting(sequence)
	if mitochondrial_targeting.count(1) > 10:
		cur_consecutive = 0
		best=0
		mito_seq = ''
		best_mito = ''
		for i in range(0, len(sequence)):
			cur_mito = mitochondrial_targeting[i]
			cur_aa = sequence[i]
			if cur_mito==1:
				cur_consecutive += 1
				mito_seq += cur_aa
			else:
				if cur_consecutive > best:
					best = cur_consecutive
					best_mito = mito_seq
				cur_consecutive = 0
				mito_seq=''
		sec_loc = sequence.index(best_mito)+1
		end_loc = sec_loc+len(best_mito)
		mito_loc = {best_mito:[sec_loc, end_loc]}
	
	nes_targeting = _predict_nes_seq(sequence)
	num_consec=0
	potential_nes = ''
	nes_seqs = []
	for prob in range(0, len(nes_targeting)):
		cur_prob = nes_targeting[prob]
		cur_aa = sequence[prob]
		if cur_prob > 0.5:
			num_consec += 1
			potential_nes += cur_aa
		else:
			if num_consec > 4:
				nes_seqs.append(potential_nes)
			potential_nes = ''
			num_consec = 0

	for seq in nes_seqs:
		cur_start = sequence.index(seq)+1
		cur_end = cur_start + len(seq)

		nes_loc[seq]=[cur_start, cur_end]

	nls_targeting = _predict_nls_seq(sequence)
	num_consec=0
	potential_nls = ''
	nls_seqs = []
	for prob in range(0, len(nls_targeting)):
		cur_prob = nls_targeting[prob]
		cur_aa = sequence[prob]
		if cur_prob > 0.5:
			num_consec += 1
			potential_nls += cur_aa
		else:
			if num_consec > 4:
				nls_seqs.append(potential_nls)
			potential_nls = ''
			num_consec = 0

	for seq in nls_seqs:
		cur_start = sequence.index(seq)+1
		cur_end = cur_start + len(seq)

		nls_loc[seq]=[cur_start, cur_end]

	if mito_loc != {}:
		all_locs['mitochondria'] = mito_loc
	else:
		all_locs['mitochondria'] = 'No mitochondrial targeting sequences predicted.'	

	if nes_loc != {}:
		all_locs['NES'] = nes_loc
	else:
		all_locs['NES'] = 'No NES sequences predicted.'

	if nls_loc != {}:
		all_locs['NLS']=nls_loc
	else:
		all_locs['NLS'] = 'No NLS targeting sequences predicted.'
	
	return all_locs


def transcriptional_activation(sequence):
	'''
	function to return any potential transcriptional activation domains

	parameters
	---------
	sequence : string
		amino acid sequence as a string

	returns 
	-------
	all_locs : dict
		a dictionary of the tad sequences and their coordinates in the seq.
	'''
	tad_locs = _predict_tad_seq(sequence)
	num_consec=0
	potential_tad = ''
	tad_seqs = []
	for prob in range(0, len(tad_locs)):
		cur_prob = tad_locs[prob]
		cur_aa = sequence[prob]
		if cur_prob > 0.5:
			num_consec += 1
			potential_tad += cur_aa
		else:
			if num_consec > 5:
				tad_seqs.append(potential_tad)
			potential_tad = ''
			num_consec = 0

	all_locs = {}
	if tad_seqs != []:
		for seq in tad_seqs:
			start_seq = sequence.index(seq)+1
			end_seq = start_seq + len(seq)
			all_locs[seq]=[start_seq, end_seq]
	else:
		all_locs = {'None': 'No predicted locations found for TAD.'}

	return all_locs


def everything(sequence):
	'''
	for when you want all the info in one sweet go.
	parameters
	---------
	sequence : string
		amino acid sequence as a string

	returns 
	-------
	all_info : dict
		a dictionary of all of the info.

	'''
	# first get the basics
	all_info = Protein.calc_all_properties(sequence)
	all_info['kappa'] = pr(sequence).kappa
	# now get the rest
	all_info['predicted phosphosites'] = phosphosites(sequence)
	all_info['predicted cellular localization'] = cellular_localization(sequence)
	all_info['predicted transcriptional activation'] = transcriptional_activation(sequence)
	all_info['fractions'] = Protein.calc_frac(sequence)
	all_info['sequence'] = sequence

	return all_info

