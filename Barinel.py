"""
This module implements BARINEL algorithm, which generates diagnoses and rank them
according to their probability of being correct, given spectrum-matrix and error-vector.

BARINEL is described in
	Abreu, R., Zoeteweij, P., & Van Gemund, A. J. (2009, November).
	Spectrum-based multiple fault localization.

The fuzzy logic addition is described in
	Cardoso, N., Abreu, R., Feldman, A., & Kleer, J. D. (2016, August).
	A framework for automatic debugging of functional and degradation failures.

This module does not implement BARINEL exactly as described in the paper.
It might not be stable, and suffer some flaws.
This code is not suitable for distribution.
"""

__author__ = "Shay Lapid"
__email__ = "lapidshay@gmail.com"

import numpy as np
from scipy.optimize import minimize
from Staccato import staccato, tuple_add_one

######################################
# global constant a-priori probability for a component failure
######################################

P = 0.05


######################################
# Utility functions
######################################

def rank_pretty_print(ranked_probs, k_top=20):
	"""
	visualization tool
	"""
	output = '\nBarinel ranked probabilities:'
	for i, tup in enumerate(ranked_probs):
		if i < k_top:
			tmp = f'\n\t{i+1}. Diagnosis {tuple_add_one(tup[0])}:'.ljust(30)
			tmp += f'Probability: {tup[1]:.4f}'
			output += tmp
	return output


def mat_pretty_print(mat, vec):
	"""
	visualization tool
	"""
	output = '\nSpectra and Error vector:'
	e_string = '|  e' if vec.dtype == np.dtype('float64') else '| e'
	output += '\n\t' + ' '.join([f'c{i+1}'.ljust(3) for i in range(mat.shape[1])]) + e_string
	for row, e in zip(mat, vec):
		output += '\n\t'+' '.join([f' {c} ' for c in row]) + f'| {e[0]}'
	return output


def load_csv(matrix_path):
	"""

	:param matrix_path: path of a csv file
	:return: spectra-matrix, error-vector
	"""
	import pandas as pd

	full_matrix = pd.read_csv(matrix_path, index_col=0).values
	return full_matrix[:, :-1].astype(np.int32), full_matrix[:, -1].reshape(-1,1)

######################################
# Barinel sub-functions
######################################

def no_obs_diag_pr(diag: tuple, m):
	"""
	calculate the probability of a diagnosis before any observation is observed.
	uses P global variable.

	:param diag: a diagnosis tuple
	:param m: number of components in system
	:return: probability
	"""

	diag_len = len(diag)

	return np.power(P, diag_len) * np.power(1-P, m-diag_len)


def ochiai_arguments(c, mat, vec):

	# extract component spectra (in which test it participated (1) and not (0))
	comp_spectra = mat[:, c].reshape(-1, 1)

	# compute Ochiai arguments
	n11 = np.sum(comp_spectra * vec)
	n10 = np.sum(comp_spectra * (1-vec))
	n01 = np.sum((1-comp_spectra) * vec)

	return n11, n10, n01


def ochiai_rank(c, mat, vec):
	# c is index of a component

	# compute Ochiai arguments
	n11, n10, n01 = ochiai_arguments(c, mat, vec)

	# compute Ochiai similarity
	if (n11 + n10) * (n11 + n01) != 0:
		return n11/np.sqrt((n11 + n10) * (n11 + n01))

	return 0


def multi_fault_maximize_probs(mat, vec, diags):
	"""
	Main Barinel maximization function.
	:param mat: mat: spectra matrix.
	:param vec: fuzzy error vector.
	:param diags: diagnosis candidates.
	:return: MLE maximized probabilities.
	"""
	def maximize(x):
		"""
		to be used in a scipy.optimize.minimize, as function.
		using outer-scope global parameters: error vector PARTIAL_VEC, and spectrum matrix PARTIAL_MAT
		:param x: initial guess = a list of a-priori probabilites of the diagnosis' components
		"""

		# multiply components occurrences with their health
		mat = np.multiply(PARTIAL_MAT, x)

		pr = []

		for test, obs in enumerate(PARTIAL_VEC):
			# extract participating components indices (from a given diagnosis) in a given test
			partic_ind = np.where(mat[test] != 0)[0]

			# multiply components' health coefficient
			prod_tmp = np.prod(mat[test][partic_ind])

			# fuzzy epsilon policy - Pr(obs|d)=e_i*(1-G(d,A_i)) + (1-e_i)*G(d,A_i)
			pr.append(obs * (1-prod_tmp) + (1-obs)*prod_tmp)

		# multiply all elements. minus is for maximizing instead of minimizing
		return -np.prod(pr)

	max_probs = []
	for diag in diags:
		# set a-priori probabilities to all diagnosis' components
		init_guess = np.ones(len(diag)) * P

		# create a partial matrix of only diagnosis' participant components
		PARTIAL_MAT = mat[:, diag]

		# extract indices of rows to leave (rows who contain at least one of diag's components)
		rows_to_include = np.where(np.sum(PARTIAL_MAT, axis=1) != 0)[0]

		# reduce partial matrix and partial vector according to rows to include
		PARTIAL_MAT = PARTIAL_MAT[rows_to_include, :]
		PARTIAL_VEC = vec[rows_to_include, :]

		# set tuples of ranges as boundaries for optimization
		boundaries = [(0, 1) for _ in range(len(diag))]

		# maximize probability
		sol = minimize(maximize, init_guess, method='L-BFGS-B', bounds=boundaries)
		max_probs.append(-sol.fun)

	return max_probs


######################################
# Barinel main algorithm
######################################

def barinel(mat, vec, diags, only_single=False):
	"""
	rate minimal-hitting-sets by their probability to be the correct diagnosis.
	STACCATO is deliberately outside this method (although it should be), for ease of use in out experiment.
	one can easily change the code so STACCATO will be called from here.
	:param mat: spectra-matrix
	:param vec: error vector (fuzzy or crisp)
	:param diags: minimal-hitting-sets, produced by STACCATO
	:param only_single: use this option if you can be sure there is only single-fault.
	:return: sorted list of diagnoses and their probability to be the correct ones
	"""

	# number of components
	num_comps = mat.shape[1]

	# not in use for the project
	if only_single:
		probs = [(diag, ochiai_rank(diag, mat, vec)) for diag in diags]

	else:
		probs = multi_fault_maximize_probs(mat, vec, diags)

	# multiply each diagnosis' post-priori probability with it's a-priori probability
	no_obs_probs = np.array([no_obs_diag_pr(diag, num_comps) for diag in diags])
	mxmul = probs * no_obs_probs

	# normalize probabilities
	normalized_probs = mxmul/np.sum(mxmul)

	output = sorted([(diag, prob) for diag, prob in zip(diags, normalized_probs)], key=lambda x: x[1], reverse=True)

	return output


def main():

	# load a matrix
	mat, vec = load_csv('spectra_error_example.csv')
	print(mat_pretty_print(mat, vec))

	# create hitting-sets using STACCATO
	stac_diags = staccato(mat, vec)

	# rank diagnoses using BARINEL
	ranked_diags = barinel(mat, vec, stac_diags, only_single=False)
	print(rank_pretty_print(ranked_diags))


if __name__ == '__main__':
	main()
