"""
This module implements STACCATO algorithm, which generates minimal hitting-sets,
given spectrum-matrix and error-vector.

STACCATO is described in
	Abreu, R., & van Gemund, A. J. (2009, June).
	Statistics-directed minimal hitting set algorithm.

This module does not implement STACCATO exactly as described in the paper.
It might not be stable, and suffer some flaws.
This code is not suitable for distribution.
"""

__author__ = "Shay Lapid"
__email__ = "lapidshay@gmail.com"


import numpy as np
from datetime import datetime
from itertools import combinations


######################################
# Ochiai
######################################

def compute_component_ochiai_arguments(c, mat, vec, comps, verbose=False):

	# extract component spectra (in which test it participated (1) and not (0))
	comp_spectra = mat[:, c].reshape(-1,1)

	# compute Ochiai arguments
	n11 = np.sum(comp_spectra * vec)
	n10 = np.sum(comp_spectra * (1-vec))
	n01 = np.sum((1-comp_spectra) * vec)

	if verbose:
		print(f'Component {comps[c]+1}: n11: {n11} ;   n10: {n10} ;   n01: {n01}')

	return n11, n10, n01


def compute_ochiai_similiarity(c, mat, vec, comps):
	# c is index of a component

	# compute Ochiai arguments
	n11, n10, n01 = compute_component_ochiai_arguments(c, mat, vec, comps)

	# compute Ochiai similarity
	if (n11 + n10) * (n11 + n01) != 0:
		return {'och': n11/np.sqrt((n11 + n10) * (n11 + n01)), 'n11': n11}
	return {'och': 0, 'n11': n11}


def not_superset(sut_sup: set, set_list: list):
	# return True if sut_sup is not a superset of any element in set_list, False otherwise
	return not any(sut_sup.issuperset(elem) for elem in set_list)


######################################
# STACCATO sub-functions
######################################

def strip_component(mat, comps, ranking):
	# remove components columns if they appeared in all tests (i.e)
	return mat[:, sorted([comps[i] for i in ranking])].copy()


def strip(mat, vec, comps, c, verbose=False):

	# extract component from set
	c = list(c)[0]

	# determine which rows and columns to not remove
	conflict_not_involved = np.where(mat[:, comps[c]] == 0)[0]
	other_comps = [comps[i] for i in comps if i != c]

	# reduce the matrix and vector
	mat_1 = mat[conflict_not_involved, :][:, other_comps].copy()
	vec_1 = vec[conflict_not_involved].copy()

	if verbose:
		print(f'\ninside stip() function:')
		print(f'\tcomponent to strip: {c}')
		print(f'\tconflict_not_involved: {conflict_not_involved}')
		print(f'\tcomps: {comps}')
		print(f'\tother_comps: {other_comps}')

	return mat_1, vec_1


def update_comp_indices_single(comps_dict, removed_comp):
	for c in comps_dict:
		if c >= removed_comp:
			comps_dict[c] -= 1


def update_comp_indices_all(comps_dict, removed_comps):
	for comp in removed_comps:
		c = list(comp)[0]
		update_comp_indices_single(comps_dict, c)
		del(comps_dict[c])


def is_not_superset(set_to_check: set, set_list: list):
	# return True if set_to_check is not a superset of any element in set_list, False otherwise
	return not any(set_to_check.issuperset(elem) for elem in set_list)


def remove_supersets(diags):
	output = []
	for diag in diags:
		cur_diags = diags.copy()
		cur_diags.remove(diag)
		if not_superset(set(diag), cur_diags):
			output.append(diag)
	return output


######################################
# original STACCATO algorithm
######################################

def staccato_reg(mat: np.array, vec: np.array, m=None, lam: float=1.0, l: int=100, comps=None):
	# create a local copy, to not make changes to original arrays
	mat_0, vec_0 = np.copy(mat), np.copy(vec)

	# extract number of observed conflicts and number of components
	num_conflicts = np.sum(vec_0)

	# extract original number of components
	if m is None:
		num_comps = mat_0.shape[1]
	else:
		num_comps = m

	# create an indexing map in 1st iteration, and use it in next iterations
	if comps is None:
		comps = {c: c for c in range(num_comps)}
	else:
		comps = comps.copy()

	# create a dictionary with och rank and n11 argument of each component
	och_rank = {c: compute_ochiai_similiarity(comps[c], mat_0, vec_0, comps) for c in comps}

	# sort components by Ochiai similarity rank
	ranking = sorted(och_rank, key=lambda r: och_rank[r]['och'], reverse=True)

	# instantiate limit counter and output
	seen_comps = 0
	output = []

	# add single-fault diagnoses to output and remove from ranked components list
	for c in comps:
		if och_rank[c]['n11'] == num_conflicts:
			output.append({c})
			ranking.remove(c)
			seen_comps += 1/num_comps

	# remove components' columns who appeared in all tests
	mat_1 = strip_component(mat_0, comps, ranking)

	# update the indexing map
	update_comp_indices_all(comps, output)

	comps_c = comps.copy()

	while len(ranking) > 0 and seen_comps < lam and len(output) < l:

		mat_2 = mat_1.copy()
		seen_comps += 1 / num_comps
		unite_comp = {ranking.pop(0)}

		mat_hat, vec_hat = strip(mat_2, vec_0, comps_c, unite_comp)

		comps_c2 = comps_c.copy()
		update_comp_indices_all(comps_c2, [unite_comp])

		output_hat = staccato_reg(mat_hat, vec_hat, m=num_comps, lam=lam-seen_comps, l=l-len(output), comps=comps_c2)

		while output_hat and len(output_hat) > 0:
			cur_comp = output_hat.pop(0)
			cur_comp |= unite_comp

			if is_not_superset(cur_comp, output):
				output.append(cur_comp)

	output = remove_supersets(output)  # substitute with subsumed()?

	return output


######################################
# STACCATO additions
######################################

def reduce_hs_to_mhs(mat, vec, hs):

	# if hitting-set is too large for Brute-Force, don't try to reduce, due to complexity
	if len(hs) > 15:
		return hs

	# reduce the matrix to only rows containing failed test and columns included in hitting-set
	partial_failed_spectra = (mat[np.where(vec)[0], :][:, hs])

	# map the components indices, by creating {ind:comp} dictionary
	comp_indices = {ind: comp for ind, comp in zip(np.arange(len(hs)), hs)}

	# create all sizes combination from hitting-set's components
	comps_combinations = []
	for i in range(len(hs)+1):
		comps_combinations += combinations(comp_indices, i)

	# for every combination, check if it's a hitting-set and is not longer than previous hitting-sets,
	# and map back to original components
	output = []
	min_len_hs = np.inf
	for comb in comps_combinations:
		if partial_mat_is_hs(comb, partial_failed_spectra):
			hit_set = tuple([comp_indices[c] for c in comb])
			if len(hit_set) <= min_len_hs:
				output.append(hit_set)
				min_len_hs = len(hit_set)

	return output


def partial_mat_is_hs(set_to_check, partial_matrix):
	"""
	:param set_to_check: set ot tuple of components (indices) to check
	:param partial_matrix: a matrix composed of only failed tests
	:return: True if the set is a hitting-set, False otherwise
	"""
	return True if all(np.sum(partial_matrix[:, set_to_check], axis=1)) >= 1 else False


def tuple_add_one(tup):
	tup = tuple(tup)
	# add 1 to each element of a tuple
	return tuple(tup[i]+1 for i in range(len(tup)))


def staccato(matrix, vector, lam=1, l=100):
	"""
	generates minimal hitting-sets, given spectrum-matrix and error-vector.
	:param matrix: spectra-matrix
	:param vector: error vector (fuzzy or crisp)
	:param lam, l: tools for early stopping
	:return: minimal-hititng-sets
	"""

	# perform regular STACCATO
	start_time = datetime.now()
	stac_diags = staccato_reg(matrix, vector, lam=lam, l=l)
	comp_time_1 = datetime.now()-start_time

	# reduce hitting-sets to minimal hitting-sets
	upgraded_stac_diags = []
	for diag in stac_diags:
		upgraded_stac_diags += reduce_hs_to_mhs(matrix, vector, tuple(diag))

	# sort each diagnosis in ascending order
	upgraded_stac_diags = [tuple(sorted(diag)) for diag in upgraded_stac_diags]

	# sort diagnoses list by cardinality in ascending order
	upgraded_stac_diags = sorted(upgraded_stac_diags, key=lambda x: x[0])

	return upgraded_stac_diags
