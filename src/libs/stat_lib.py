#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
@author: Flavio Lichtenstein
Updated on 2025/10/28, 2024/05/24
Created on 2010/06/10

@email:  flalix@gmail.com, flavio.lichtenstein@butantan.gov.br
@local:  Instituto Butantan / CENTD / Molecular Biology / Bioinformatics & Systems Biology
'''

from typing import Tuple, List # Optional, Iterable, Set, Any, List

import pandas as pd
import numpy as np

from scipy.special import factorial
from scipy import stats

import scikit_posthocs as sposthocs

from   statsmodels.formula.api import ols

from libs.Basic import *

def standardize_log2_x2_table(df, numfields) -> pd.DataFrame:
	df2 = df.copy()

	for col in numfields:
		df2["log2_"+col] = [np.log2(x) if x > 0 else -10 for x in df2[col]]
		df2["x2_"+col]   = [x*x for x in df2[col]]

		mu = df2[col].mean()
		sd = df2[col].std()
		df2["std_"+col] = (df2[col]- mu) / sd

	return df2

def standardize_table(df, numfields)-> pd.DataFrame:
	df2 = df.copy()

	for col in numfields:
		mu = df2[col].mean()
		sd = df2[col].std()

		df2["std_"+col] = (df2[col]- mu) / sd

	return df2

def calc_stat_chi2_mat(mat:List) -> Tuple[str, float, float, int, List]:

	try:
		res = stats.chi2_contingency(mat)

		statistic = float(res.statistic)
		pvalue = float(res.pvalue)
		dof = int(res.dof)
		expected = res.expected_freq
	except:
		s_error = f"Error: chi2_contingency {mat}"
		print(s_error)
		return s_error, -1, -1, 0, []
	

	s_stat = f"chi-square statistics = {statistic:5.3f}, p-value = {pvalue:.2e} dof = {dof}"

	return  s_stat, statistic, pvalue, dof, expected


def calc_stat_chi2(vals0:List, vals1:List) -> Tuple[str, float, float, int, List]:

	vals0 = list(vals0)
	vals1 = list(vals1)
	
	#--- at leat a table 2x2 -------------
	if len(vals0) < 2 or len(vals1) < 2:
		s_error = f"Error: len(val) < 2"
		print(s_error)
		return s_error, -1, -1, 0, []

	return calc_stat_chi2_mat([vals0, vals1])



def calc_fisher_exact_test(vals0:List, vals1:List) ->  Tuple[str, float, float, int, List]:

	vals0 = list(vals0)
	vals1 = list(vals1)
	
	if len(vals0) == 0 or len(vals0) != len(vals1):
		s_error = f"Error in vals0 {len(vals0)} x vals1 {len(vals1)}"
		print(s_error)
		return s_error, -1, -1, -1, []
	
	try:
		res = stats.fisher_exact([vals0, vals1])
		statistic, pvalue = float(float(res.statistic)), float(float(res.pvalue))
	except:
		statistic, pvalue = -1, -1
	
	s_stat = f"Fisher's exact test = {statistic:5.3f}, p-value = {pvalue:.2e}"
	
	return  s_stat, statistic, pvalue, 1, []

def chi2_or_fisher_exact_test(vals0:List, vals1:List) -> Tuple[str, float, float, int, List]:
	vals0 = list(vals0)
	vals1 = list(vals1)
	
	less5 = False
	for val in vals0:
		if val < 6:
			less5 = True
			break

	if not less5:
		for val in vals1:
			if val < 6:
				less5 = True
				break
	
	if less5 & len(vals0)==2 & len(vals1)==2:
		return calc_fisher_exact_test(vals0, vals1)
	
	return calc_stat_chi2(vals0, vals1)


def chisquare_gender01_female_male(df2: pd.DataFrame, measure:str, field:str='gender', alpha:float=0.05,
								   verbose:bool=False) -> Tuple[bool, str, pd.DataFrame, List, bool, int, float, float, str]:

	typetest = 'chi-square'

	valsf = df2[df2[field] == 0].groupby(measure).count().iloc[:,0:1]
	valsf['class'] = valsf.index
	valsf.columns = ['count', 'class']
	valsf = valsf[['class', 'count']].T

	valsm = df2[df2[field] == 1].groupby(measure).count().iloc[:,0:1]
	valsm['class'] = valsm.index
	valsm.columns = ['count', 'class']
	valsm = valsm[['class', 'count']].T

	cols = np.unique(list(valsf.columns) + list(valsm.columns))
	cols = np.sort(cols)
	dic = {}
	for col in cols:
		dic[col] = []

		try:
			dic[col].append(valsf[col][1])
		except:
			dic[col].append(0)

		try:
			dic[col].append(valsm[col][1])
		except:
			dic[col].append(0)

	dfgender = pd.DataFrame.from_dict(dic)
	# dfgender.index = ['female', 'male']
	dfgender.rename(index={0: 'female', 1: 'male'}, inplace=True)

	for col in dfgender.columns:
		dfgender[col] = dfgender[col].astype(int)
	nfemales = len([x for x in dfgender.iloc[0] if x != 0])
	nmales   = len([x for x in dfgender.iloc[1] if x != 0])

	sumf = np.sum( [True if x <= 5 else False for x in dfgender.iloc[0] if x != 0])
	summ = np.sum( [True if x <= 5 else False for x in dfgender.iloc[1] if x != 0])

	if sumf == 0 and summ == 0:
		add5 = False
	else:
		dfgender.loc[0, :] = [x+5 for x in dfgender.iloc[0] ]
		dfgender.loc[1, :] = [x+5 for x in dfgender.iloc[1] ]
		add5 = True

	#statistic, pvalue = stats.chisquare(np.array([dfgender.iloc[0].to_list(), dfgender.iloc[1].to_list()]), axis=None)
	try:
		statistic, pvalue, dof, _ = stats.chi2_contingency([dfgender.iloc[0].to_list(), dfgender.iloc[1].to_list()])
		statistic = float(statistic)
		pvalue = float(pvalue)
		dof = int(dof)
	except:
		s_error = f'Error chi2_contingency'
		print(s_error)
		return False, typetest, dfgender, [nfemales, nmales], add5, -1, -1, -1, s_error
			

	# statistic, pvalue = fisher_exact( [dfgender.iloc[0].to_list(), dfgender.iloc[1].to_list()] )  # statistic is oddsratio
	# typetest = 'fisher exact test'
	
	ret = True if pvalue < alpha else False
	stri = stat_asteristics(pvalue)
	stri_stat = f'p-value {pvalue:.2e}<br>{stri}'

	if verbose:
		print(dfgender)
		print(stri_stat)

	return ret, typetest, dfgender, [nfemales, nmales], add5, dof, statistic, pvalue, stri_stat


def chisquare_2series(df1:pd.DataFrame, df2:pd.DataFrame, measure, alpha = 0.05, verbose=False):
	try:
		vals1 = df1.groupby(measure).count().iloc[:,0:1]
	except:
		return [None]*9

	vals1['class'] = vals1.index
	vals1.columns = ['count', 'class']
	vals1 = vals1[['class', 'count']].T

	try:
		vals2 = df2.groupby(measure).count().iloc[:,0:1]
	except:
		return [None]*9

	vals2['class'] = vals2.index
	vals2.columns = ['count', 'class']
	vals2 = vals2[['class', 'count']].T

	cols = np.unique(list(vals1.columns) + list(vals2.columns))
	cols = np.sort(cols)
	dic = {}
	for col in cols:
		dic[col] = []

		try:
			dic[col].append(vals1[col][1])
		except:
			dic[col].append(0)

		try:
			dic[col].append(vals2[col][1])
		except:
			dic[col].append(0)

	dfgender = pd.DataFrame.from_dict(dic)
	dfgender.index = [0, 1]
	for col in dfgender.columns:
		dfgender[col] = dfgender[col].astype(int)
	nfemales = len([x for x in dfgender.iloc[0] if x != 0])
	nmales   = len([x for x in dfgender.iloc[1] if x != 0])

	sumf = np.sum( [True if x <= 5 else False for x in dfgender.iloc[0] if x != 0])
	summ = np.sum( [True if x <= 5 else False for x in dfgender.iloc[1] if x != 0])

	if sumf == 0 and summ == 0:
		add5 = False
	else:
		dfgender.loc[0, :] = [x+5 for x in dfgender.iloc[0] ]
		dfgender.loc[1, :] = [x+5 for x in dfgender.iloc[1] ]
		add5 = True

	#statistic, pvalue = stats.chisquare(np.array([dfgender.iloc[0].to_list(), dfgender.iloc[1].to_list()]), axis=None)
	statistic, pvalue, dof, _ = stats.chi2_contingency([dfgender.iloc[0].to_list(), dfgender.iloc[1].to_list()])
	typetest = 'chi-square'

	# statistic, pvalue = fisher_exact( [dfgender.iloc[0].to_list(), dfgender.iloc[1].to_list()] )  # statistic is oddsratio
	# typetest = 'fisher exact test'

	if pvalue < alpha:
		ret = True
	else:
		ret = False

	stri_stat = 'p-value %.2e<br>%s'%(pvalue, stat_asteristics(pvalue))

	if verbose:
		print(dfgender)
		print(stri_stat)

	return ret, typetest, dfgender, [nfemales, nmales], add5, dof, statistic, pvalue, stri_stat

def chisquare_2by2(df:pd.DataFrame, alpha:float=0.05, verbose:bool=False) -> Tuple[bool, int, float, float, str]:

	try:
		typetest = 'chi-square'
		# statistic, pvalue, dof, expected_freq
		statistic, pvalue, dof, _ = stats.chi2_contingency([df.iloc[0].to_list(), df.iloc[1].to_list()])
		ret = True if pvalue < alpha else False
	except:
		typetest = 'chi-square Error'
		statistic = -1; pvalue = 1; dof = 1
		ret = False

	stri_stat = f"{typetest} pval {pvalue:.2e} {stat_asteristics(pvalue)}"
	if verbose: print(stri_stat)

	return ret, dof, statistic, pvalue, stri_stat


def calc_mannwhitneyu(vals1: List, vals2:List, use_continuity:bool=True, alternative:str='two-sided', 
					  alpha:float=0.05, language:str='English') -> Tuple[str, float, float]:

	if len(vals1) < 2 or len(vals2) < 2:
		s_error = f'Insuficient data to calculate t-test {len(vals1)} x {len(vals2)}'
		print(s_error)
		return s_error, -1, -1

	statistic, pval = stats.mannwhitneyu(vals1, vals2, use_continuity=use_continuity, alternative=alternative)

	text_stat = "t-test statistic = %.2f, p-val = %.2e"%(statistic, pval)

	if alternative == 'two-sided':
		alpha /= 2

	if pval < alpha:
		if language == 'English':
			text_stat = 'Statiscally distinct distributions, H0 must be refuted ' + text_stat
		else:
			text_stat = 'Distribuições estatiscamente diferentes, H0 precisa ser refutado\n' + text_stat
	else:
		if language == 'English':
			text_stat = 'Statiscally similar distributions, H0 must be accepetd ' + text_stat
		else:
			text_stat = 'Distribuições estatiscamente similares, H0 aceita\n' + text_stat

	return text_stat, statistic, pval


def calc_ttest(vals1:List, vals2:List, alternative:str='two-sided', 
			   alpha:float=0.05, language:str='English') -> Tuple[str, float, float]:
	'''
		calculus of ttest
		input 2 val list
		removed: equal_var:bool=True
	'''

	vals1 = list(vals1)
	vals2 = list(vals2)

	if len(vals1) < 2 or len(vals2) < 2:
		return 'Insuficient data to calculate t-test', -1, -1

	try:
		res = stats.ttest_ind(vals1, vals2)
		statistic = float(res.statistic)
		pval = float(res.pvalue)
		df   = int(res.df)
		
	except:
		s_error = f"Error: ttest {vals1} {vals2}"
		print(s_error)
		return s_error, -1, -1

	text_stat = f"t-test statistic = {statistic:.3f}, pval = {pval:.2e}, and DF = {df}"

	if alternative == 'two-sided':
		alpha /= 2

	if pval < alpha:
		if language == 'English':
			text_stat = 'Statiscally distinct distributions, H0 must be refuted: ' + text_stat
		else:
			text_stat = 'Distribuições estatiscamente diferentes, H0 precisa ser refutado: ' + text_stat
	else:
		if language == 'English':
			text_stat = 'Statiscally similar distributions, H0 must be accepetd: ' + text_stat
		else:
			text_stat = 'Distribuições estatiscamente similares, H0 aceita: n' + text_stat

	return text_stat, statistic, pval

def ttest_2series(df1: pd.DataFrame, df2: pd.DataFrame, measure:str, alpha:float=0.05, 
				  verbose:bool=False) -> Tuple[bool, List, List, float, float, str]:

	vals1 = df1[measure].to_list()
	vals1 = [x for x in vals1 if not pd.isnull(x) ]
	n1 = len(vals1)

	vals2 = df2[measure].to_list()
	vals2 = [x for x in vals2 if not  pd.isnull(x) ]
	n2 = len(vals2)

	if n1 < 3 or n2 < 3:
		s_error = f'Error: insuficient number of vals: {n1} x {n2}'
		print(s_error)
		return False, [vals1, vals2], [n1, n2], -1, -1, s_error


	try:
		res = stats.ttest_ind(vals1, vals2)
		statistic = float(res.statistic)
		pvalue = float(res.pvalue)
		df   = int(res.df)
	except:
		s_error = f"Error: ttest {vals1} {vals2}"
		print(s_error)
		return False, [vals1, vals2], [n1, n2], -1, -1, s_error
	
	if pvalue < alpha:
		ret = True
	else:
		ret = False

	stri_stat = 'p-value %.2e<br>%s'%(pvalue, stat_asteristics(pvalue))

	if verbose:
		if pvalue < alpha:
			sig = stat_asteristics(pvalue)
			print("Null hypothesis is rejected, both distributions may be different.")
			print(f"p-value {pvalue:.2e} significance = {sig} dof = {df}")
		else:
			print("fail to reject the null hypothesis, both distributions may be similar.")
			print('p-value %.2e'%(pvalue))

	return ret, [vals1, vals2], [n1, n2], statistic, pvalue, stri_stat

def ttest_2vals(vals1, vals2, alpha = 0.05, verbose=False) -> Tuple[bool, List, List, float, float, str]:

	vals1 = [x for x in vals1 if not pd.isnull(x)]
	n1 = len(vals1)

	vals2 = [x for x in vals2 if not pd.isnull(x)]
	n2 = len(vals2)

	if n1 < 3 or n2 < 3:
		return None, [vals1, vals2], [n1, n2], None, None, None

	try:
		res = stats.ttest_ind(vals1, vals2)
		statistic = float(res.statistic)
		pvalue = float(res.pvalue)
		df   = int(res.df)
	except:
		s_error = f"Error: ttest {vals1} {vals2}"
		print(s_error)
		return False, [vals1, vals2], [n1, n2], -1, -1, s_error

	ret = True if pvalue < alpha else False

	stri_stat = 'p-value %.2e<br>%s'%(pvalue, stat_asteristics(pvalue))

	if verbose:
		if pvalue < alpha:
			print("Null hypothesis is rejected, both distributions may be different.")
			print('p-value %.2e significance = "%s"'%(pvalue, stat_asteristics(pvalue)))
		else:
			print("fail to reject the null hypothesis, both distributions may be similar.")
			print('p-value %.2e'%(pvalue))

	return ret, [vals1, vals2], [n1, n2], statistic, pvalue, stri_stat

def ttest_1val(vals:List, popmean:float=1.00, 
			   alternative:str='two-sided') -> (float, float, int, float, float, str):
	
	vals = list(vals)
	mu  = np.mean(vals)
	std = np.std(vals)

	try:
		res = stats.ttest_1samp(vals, popmean=popmean, alternative=alternative)
		statistic = float(res.statistic)
		pvalue = float(res.pvalue)
		df   = int(res.df)
	except:
		s_error = f"Error: ttest {vals}"
		print(s_error)
		return -1, -1, 0, mu, std, s_error

	msg = f"mean={mu:.3f} and std={std:.3f}, pval = {pval:.3e} and dof={df}"

	return statistic, pval, dof, mu, std, msg



def ttest_gender01_female_male(df2:pd.DataFrame, measure:str, field:str='gender', 
							   alpha:float=0.05, verbose:bool=False) -> Tuple[bool, List, List, float, float, str]:
	valsf = df2[df2[field] == 0][measure].to_list()
	valsf = [x for x in valsf if not try_float(x) is None ]
	nfemales = len(valsf)

	valsm = df2[df2[field] == 1][measure].to_list()
	valsm = [x for x in valsm if not try_float(x) is None]
	nmales = len(valsm)

	if nfemales < 3 or nmales < 3:
		return None, [valsf, valsm], [nfemales, nmales], None, None, None

	try:
		res = stats.ttest_ind(valsf, valsm)
		statistic = float(res.statistic)
		pvalue = float(res.pvalue)
		df   = int(res.df)
	except:
		s_error = f"Error: ttest {vals1} {vals2}"
		print(s_error)
		return False, [valsf, valsm], [nfemales, nmales], -1, -1, s_error


	ret = True if pvalue < alpha else False

	aster = stat_asteristics(pvalue)
	stri_stat = f'p-value {pvalue:2e}<br>{aster}'

	if verbose:
		print(stri_stat)

	return ret, [valsf, valsm], [nfemales, nmales], statistic, pvalue, stri_stat


def test_normality(vals:List, alpha:float=0.05, sample_limit:int=3) -> Tuple[bool, float, float, str]:
	'''
		at least 3 different points
	'''
	vals = [x for x in vals if not pd.isnull(x)]

	if len(vals) < sample_limit:
		return False, -1, -1, f'No sufficient data different values n={len(vals)}'
	if len(np.unique(vals)) == 1:
		return False, -1, -1, f'All data have the same values.'

	# normality test
	# print(">>>", len(vals), vals)
	try:
		statistic, pvalue = stats.shapiro(vals)
	except:
		s_error = f'Error to calculate shapiro {len(vals)}'
		print(s_error)
		return False, -1, -1, s_error

	normal = True if pvalue >= alpha else False

	aster = stat_asteristics(pvalue)
	stri_stat = f'p-value {pvalue:.2e}<br>{aster}'

	return normal, statistic, pvalue, stri_stat

def test_normality_params(vals:List, alpha:float=0.05):
	normal, statistic, pvalue, stri_stat = test_normality(vals=vals, alpha=alpha)

	if statistic is None:
		return None, None, None, None, normal, statistic, pvalue, stri_stat

	mu	   = np.mean(vals)
	median = np.median(vals)
	std	   = np.std(vals)
	n	   = len(vals)

	return mu, median, std, n, normal, statistic, pvalue, stri_stat

def test_normality_desc_eng(vals, alpha = 0.05, verbose=False):
	if len(vals) < 3:
		if verbose: print('Sample does not look Gaussian (reject H0), length < 3')
		stri_stat = 'p-value - impossible to evaluate, length < 3'
		return False, 999, 0.000001, stri_stat

	# normality test
	statistic, pvalue = stats.shapiro(vals)
	if verbose: print('Statistics=%.3f, p-value=%.3f' % (statistic, pvalue))
	# interpret
	if pvalue > alpha:
		if verbose: print('Sample looks Gaussian (fail to reject H0)')
		ret = True
	else:
		if verbose: print('Sample does not look Gaussian (reject H0)')
		ret = False

	aster = stat_asteristics(pvalue)
	stri_stat = f'p-value {pvalue:.2e}<br>{aster}'

	return ret, statistic, pvalue, stri_stat

def test_normality_desc(vals, alpha = 0.05, verbose=False):
	if len(vals) < 3:
		text = 'A distribuição não se assemelha a uma distriuição normal (rejeita-se a H0) - menos que 3 valores.'
		text_stat = 'p-value - impossible to evaluate, length < 3'
		if verbose: print(text, '\n', text_stat)
		return False, text, text_stat, 999, 0.000001

	# teste de normalidade de Shapiro-Wilkis
	statistic, pvalue = stats.shapiro(vals)

	if pvalue > alpha:
		text = 'A distribuição se assemelha a uma distriuição normal (não se rejeita a H0)'
		ret = True
	else:
		text = 'A distribuição não se assemelha a uma distriuição normal (rejeita-se a H0)'
		ret = False

	text_stat = 'p-value %.2e (%s)'%(pvalue, stat_asteristics(pvalue))
	if verbose: print(text, '\n', text_stat)

	return ret, text, text_stat, statistic, pvalue

def calc_vc(mu, std):
	if mu is None or mu == 0:
		return None

	VC = std/mu

	return VC


def calc_confidence_interval(vals:List, alpha:float=0.05, two_tailed:bool=True):
	vals = [x for x in vals if not pd.isnull(x)]
	n = len(vals)

	if n < 2:
		print("Error: at least 3 values.")
		return None, None, None, None, None, None, None, None

	mu  = np.mean(vals)
	std = np.std(vals)

	SEM = std/np.sqrt(n)
	df  = n-1

	s_two_tailed = 'two-tailed' if two_tailed else 'one-tailed'
	alpha2 = alpha/2 if two_tailed else alpha

	try:
		error = stats.t.ppf(1-alpha2,df)*SEM
		cinf = mu-error
		csup = mu+error
		stri = f"mu={mu:.2f} std={std:.2f} n={n} SEM={SEM:.2f} CI=[{cinf:.2f},{csup:.2f}] for alpha={alpha2:.3f} {s_two_tailed}"
	except:
		error, cinf, csup = None, None, None
		stri = f"mu={mu:.2f} std={std:.2f}  n={n} SEM={SEM:.2f} C CI=ERROR!!!"

	return mu, std, n, SEM, error, cinf, csup, stri

'''
alpha = .1
confidence = np.round(1-alpha*2, 2)
print(alpha, confidence)
'''
def calc_confidence_interval_param(mu:float, std:float, n:int, alpha:float=0.05, two_tailed:bool=True):
	SEM = std/np.sqrt(n)
	df  = n-1

	s_two_tailed = 'two-tailed' if two_tailed else 'one-tailed'
	alpha2 = alpha/2 if two_tailed else alpha

	try:
		error = stats.t.ppf(1-alpha2,df)*SEM
		cinf = mu-error
		csup = mu+error
		stri = f"mu={mu:.2f} std={std:.2f} n={n} SEM={SEM:.2f} CI=[{cinf:.2f},{csup:.2f}] for alpha={alpha2:.3f} {s_two_tailed}"
	except:
		error, cinf, csup = None, None, None
		stri = f"mu={mu:.2f} std={std:.2f}  n={n} SEM={SEM:.2f} C CI=ERROR!!!"

	''' cinf may go below zero '''
	return error, cinf, csup, SEM, stri


def test_one_way_ANOVA (ser1, ser2, alpha = 0.05):
	# teste de variancias de Fisher - one way ANOVA (analysis of variance)
	statistic, pvalue = stats.f_oneway(samp1, samp2)

	if pvalue > alpha:
		text = 'As distribuições têm variâncias similares (não se rejeita a H0)'
		ret = True
	else:
		text = 'As distribuições não têm variâncias similares (rejeita-se a H0)'
		ret = False

	text_stat = 'p-value %.2e (%s)'%(pvalue, stat_asteristics(pvalue))

	return ret, text, text_stat, statistic, pvalue

def calc_ols_regression(df, regression_model=None, col='day', depvar=None, withBIC=False):

	if regression_model is None:
		regression_model = f'%s ~ {col}'%(depvar)

	try:
		model = ols(regression_model, data=df).fit()
	except:
		print("Regression model should be wrong: '%s'"%(regression_model))
		return None, None, None, None, None, None, None, None, None, None

	# anova_table = sm.stats.anova_lm(model, typ=2)

	array_summaries = model.summary()

	dfcoeff1 = pd.read_html(array_summaries.tables[0].as_html(),index_col=0)[0]
	dfcoeff1 = dfcoeff1.reset_index()
	dfcoeff1.columns = ['param1', 'val1', 'param2', 'val2']

	dfcoeff2 = pd.read_html(array_summaries.tables[1].as_html(),index_col=0)[0]
	dfcoeff2 = dfcoeff2.reset_index()
	dfcoeff2.columns = ['param', 'coef', 'std_err', 't', 'p-value', 'ci_inf025', 'ci_sup975']
	dfcoeff2 = dfcoeff2.iloc[1:, :]

	dfcoeff3 = pd.read_html(array_summaries.tables[2].as_html(),index_col=0)[0]
	dfcoeff3 = dfcoeff3.reset_index()
	dfcoeff3.columns = ['param1', 'val1', 'param2', 'val2']

	r2,	r2val	= dfcoeff1.iloc[0, 2:4].to_list()
	adjr2, adjr2val = dfcoeff1.iloc[1, 2:4].to_list()
	fstat, fstatval = dfcoeff1.iloc[2, 2:4].to_list()
	spval, pval	 = dfcoeff1.iloc[3, 2:4].to_list()
	bic,   bicval   = dfcoeff1.iloc[6, 2:4].to_list()

	if adjr2val < 0:
		adjr2val = r2val

	bs = [float(x) for x in dfcoeff2.iloc[:,1]]
	r = np.sqrt(adjr2val)
	try:
		if bs[1] < 0: r *= -1
	except:
		pass

	if withBIC:
		stri = "adj-r = %.3f, adj-r2 = %.3f, F-stat = %.2f, p-value = %.2e, BIC = %.1f"%(r, adjr2val, fstatval, pval, bicval)
	else:
		stri = "adj-r = %.3f, adj-r2 = %.3f, F-stat = %.2f, p-value = %.2e"%(r, adjr2val, fstatval, pval)

	return stri, adjr2val, r, bs, fstatval, pval, bicval, dfcoeff1, dfcoeff2, dfcoeff3

def stat_list_asteristics(pvals, threshold=0.05, NS='NS'):
	lista_aster = []

	for pval in pvals:
		try:
			if pval >= threshold:
				lista_aster.append('NS')
			elif pval > threshold/5:
				lista_aster.append('*')
			elif pval > threshold/50:
				lista_aster.append('**')
			else:
				lista_aster.append('***')
		except:
			lista_aster.append('???')

	return lista_aster

def stat_asteristics(pval:float, threshold:float=0.05, NS:str='NS') -> str:
	if pval < 0: return '??'

	if threshold > 0.05:
		if pval >= threshold: return NS
		if pval > 0.05:  return '*'
		if pval > 0.01:  return '**'
		if pval > 0.001:  return '***'
		return '****'
	else:
		if pval >= threshold: return NS
		if pval > threshold/5:  return '*'
		if pval > threshold/50:  return '**'
		if pval > threshold/500:  return '***'

	return '****'

def ttest(vals0, vals1, alternative='two-sided', equal_var_thr=1.5, threshold=0.05, isBonferroni=False, nround=2):

	mu1  = np.round(np.mean(vals0), nround)
	std1 = np.round(np.std(vals0), nround)

	mu2  = np.round(np.mean(vals1), nround)
	std2 = np.round(np.std(vals1), nround)

	if std1 > std2:
		frac = std1/std2
	else:
		frac = std2/std1

	equal_var = True if frac < equal_var_thr else False

	statistic, pvalue = stats.ttest_ind(vals0, vals1, alternative=alternative, equal_var=equal_var )

	s_aster = stat_asteristics(pvalue, threshold)
	stri_stat = 'p-value %.2e<br>%s'%(pvalue, s_aster)

	if isBonferroni:
		stri_stat2 = 'equal_var = %s, Bonferroni threshold=%.2e'%(equal_var, threshold)
	else:
		stri_stat2 = 'equal_var = %s, threshold=%.2e'%(equal_var, threshold)

	return mu1, std1, mu2, std2, statistic, pvalue, stri_stat, s_aster, equal_var, stri_stat2

# https://machinelearningmastery.com/how-to-code-the-students-t-test-from-scratch-in-python/
def ttest_mu_std_calc(mu1, std1, n1, mu2, std2, n2, alternative='two-sided', alpha=0.05):

	sem1, sem2 = std1/np.sqrt(n1), std2/np.sqrt(n2)
	SEM = np.sqrt(sem1**2.0 + sem2**2.0)
	t_stat = (mu1 - mu2) / SEM

	df = n1+n2-2

	if alternative == 'two-sided':
		crit = stats.t.ppf(1.0 - alpha/2, df)
		pvalue = 1.0 - stats.t.cdf(abs(t_stat), df)
	else:
		crit = stats.t.ppf(1.0 - alpha, df)
		pvalue = (1.0 - stats.t.cdf(abs(t_stat), df)) * 2.0

	return t_stat, df, crit, pvalue

def ttest_mu_std(mu1, std1, n1, mu2, std2, n2, alternative='two-sided', threshold=0.05, isBonferroni=False):

	statistic, df, crit, pvalue = ttest_mu_std_calc(mu1, std1, n1, mu2, std2, n2, alternative=alternative, alpha=threshold)

	s_aster = stat_asteristics(pvalue, threshold)
	stri_stat = 'p-value %.2e<br>%s, df=%d'%(pvalue, s_aster, df)

	if isBonferroni:
		stri_stat2 = 'Bonferroni threshold=%.2e'%(threshold)
	else:
		stri_stat2 = 'threshold=%.2e'%(threshold)

	return statistic, pvalue, stri_stat, s_aster, stri_stat2


def how_many_comparisons(npoints):
	ncomparisons = 0

	while(True):
		npoints -= 1
		if npoints == 0: break
		ncomparisons += npoints

	return ncomparisons

def how_get_params(df):
	pd.concat([df.groupby('species').temperature.mean().round(1), df.groupby('species').temperature.std().round(1)], axis=1, sort=False)
	df.groupby('species').temperature.agg(['mean', 'std', 'count']).round(2)

def how_get_params_species(df):
	dfa = df.groupby(['species']).temperature.agg(['median', 'mean', 'std', 'count', 'min', 'max']).round(2)
	dfa['species'] = dfa.index
	dfa.index = np.arange(0, dfa.shape[0])
	dfa['SEM'] = dfa['std'] / np.sqrt(dfa['count'])

	normals=[]; stats=[]; pvals=[]
	for species in dfa.species:
		vals = df[df['species'] == species]['temperature']
		normal, statistic, pvalue, stri_stat = test_normality(vals)

		normals.append(normal)
		stats.append(statistic)
		pvals.append(pvalue)

	dfa['normal'] = normals
	dfa['stat'] = stats
	dfa['pvalue'] = pvals

	dfa = dfa[['species', 'count', 'min', 'max', 'median', 'mean', 'std', 'SEM', 'normal', 'stat', 'pvalue']]
	dfa.columns = ['species', 'n', 'min', 'max', 'median', 'mean', 'std', 'SEM', 'normal', 'stat-shapiro', 'pvalue']
	dfa

def how_get_params_pops(df):
	dfa = df.groupby(['ecopop']).temperature.agg(['median', 'mean', 'std', 'count', 'min', 'max']).round(2)
	dfa['ecopop']  = dfa.index
	dfa.index = np.arange(0, dfa.shape[0])
	dfa['SEM'] = dfa['std'] / np.sqrt(dfa['count'])


	normals=[]; stats=[]; pvals=[]
	for ecopop in dfa.ecopop:
		vals = df[df['ecopop'] == ecopop]['temperature']
		normal, statistic, pvalue, stri_stat = test_normality(vals)
		normals.append(normal)
		stats.append(statistic)
		pvals.append(pvalue)

	dfa['normal'] = normals
	dfa['stat'] = stats
	dfa['pvalue'] = pvals

	dfa = dfa[['ecopop', 'count', 'min', 'max', 'median', 'mean', 'std', 'SEM', 'normal', 'stat', 'pvalue']]
	dfa.columns = ['ecopop', 'n', 'min', 'max', 'median', 'mean', 'std', 'SEM', 'normal', 'stat-shapiro', 'pvalue']
	dfa

def how_get_params_species_pops(df):
	dfa = df.groupby(['species', 'ecopop']).temperature.agg(['median', 'mean', 'std', 'count', 'min', 'max']).round(2)
	dfa['species'] = [x[0] for x in dfa.index]
	dfa['ecopop']  = [x[1] for x in dfa.index]
	dfa.index = np.arange(0, dfa.shape[0])
	dfa['SEM'] = dfa['std'] / np.sqrt(dfa['count'])

	normals=[]; stats=[]; pvals=[]
	for i in range(dfa.shape[0]):
		ecopop  = dfa.iloc[i].ecopop
		species = dfa.iloc[i].species

		vals = df[ (df['species'] == species) & (df['ecopop'] == ecopop)]['temperature']
		normal, statistic, pvalue, stri_stat = test_normality(vals)
		normals.append(normal)
		stats.append(statistic)
		pvals.append(pvalue)

	dfa['normal'] = normals
	dfa['stat'] = stats
	dfa['pvalue'] = pvals

	dfa = dfa[['species', 'ecopop', 'count', 'min', 'max', 'median', 'mean', 'std', 'SEM', 'normal', 'stat', 'pvalue']]
	dfa.columns = ['species', 'ecopop', 'n', 'min', 'max', 'median', 'mean', 'std', 'SEM', 'normal', 'stat-shapiro', 'pvalue']
	dfa


def calc_params(vals, name=None, type=None, canNeg = False, ndig=2):
	try:
		vals = list(vals)
	except:
		print("Error: vals")
		return None, None

	meds=[]; mus=[]; stds=[]; ses=[]; ns=[]; minis=[]; maxis=[]
	normals=[]; pvalues=[]; normalities=[]
	q1s=[]; q2s=[]; q3s=[]; iqs=[]; qsups=[]; qinfs=[]
	corrections=['wo_correction', 'corrected']

	for i in range(2):
		med = np.median(vals).round(ndig)
		mu  = np.mean(vals).round(ndig)
		std = np.std(vals).round(ndig)
		n   = len(vals)
		SEM  = (std / np.sqrt(n)).round(ndig)
		mini = np.min(vals).round(ndig)
		maxi = np.max(vals).round(ndig)

		try:
			q1, q2, q3 = np.percentile(vals, [25,50,75]).round(ndig)
			iq = q3 - q1
			qsup = q3 + 1.5*iq
			qinf = q1 - 1.5*iq

			if not canNeg and qinf < 0: qinf = 0
		except:
			q1, q2, q3, iq = None, None, None, None

		q1s.append(q1)
		q2s.append(q2)
		q3s.append(q3)
		iqs.append(iq)
		qsups.append(qsup)
		qinfs.append(qinf)

		if len(vals) > 2 and std > 0 and type == 'continuous':
			normal, statistic, pvalue, stri_stat = test_normality(vals)
			stri_stat = stri_stat.replace('<br>', ' ')
		else:
			ret = None; statistic=None; pvalue=None; stri_stat = ''

		meds.append(med); mus.append(mu); stds.append(std)
		ses.append(SEM); ns.append(n); minis.append(mini); maxis.append(maxi)
		normals.append(ret); pvalues.append(pvalue); normalities.append(stri_stat)

		if i == 0 and iq is not None:
			#-- only remove outliers if quantiles are different
			if q1 != q2 and q2 != q3:
				#--- removing outliers
				vals = [x for x in vals if x <= qsup]
				vals = [x for x in vals if x >= qinf]


	# print(vals, meds, mus, q2s, qsups, normals, pvalues)
	# print(meds, mus, stds, ses, ns, minis, maxis, normals, pvalues, normalities)
	dfqq =  pd.DataFrame(np.array([corrections, meds, mus, stds, ses, ns, minis, maxis,
								   q1s, q2s, q3s, iqs, qsups, qinfs,
								   normals, pvalues, normalities]).T,
						columns=['correction', 'median', 'mean', 'std', 'SEM', 'n', 'mini', 'maxi',
								 'q1', 'q2', 'q3', 'iq', 'qsup', 'qinf',
								 'normal', 'pvalue', 'normal_obs'])

	if not name is None:
		meds=[]; mus=[]; stds=[]; ses=[]; ns=[]; minis=[]; maxis=[]
		normals=[]; pvalues=[]; normalities=[]
		q1s=[]; q2s=[]; q3s=[]; iqs=[]; qsups=[]; qinfs=[]
		corrections=['wo_correction', 'corrected']

		for i in range(2):
			med = np.median(vals).round(ndig)
			mu  = np.mean(vals).round(ndig)
			std = np.std(vals).round(ndig)

			dfqq['case'] = name
			cols = list(dfqq.columns[:-1])
			dfqq = dfqq[ ['case'] + cols]

	return vals, dfqq

def calc_params_outliers(vals:List, remove_outliers:bool=True, ndig:int=4):

	vals0 = [x for x in vals if not pd.isnull(x)]

	vals = vals0
	if len(vals) < 2:
		return vals, None, None, None, None, None, None, None, None, None, None, None, None, None

	try:
		q1, q2, q3 = np.percentile(vals, [25,50,75]).round(ndig)
		iq = q3 - q1
		qsup = q3 + 1.5*iq
		qinf = q1 - 1.5*iq

		#-- only remove outliers if quantiles are different
		if remove_outliers:
			if q1 != q2 and q2 != q3:
				#----- removing outliers --------------------
				vals = [x for x in vals if x <= qsup]
				vals = [x for x in vals if x >= qinf]

				#--- recalc ---------------------------------
				q1, q2, q3 = np.percentile(vals, [25,50,75]).round(ndig)
				iq = q3 - q1
				qsup = q3 + 1.5*iq
				qinf = q1 - 1.5*iq
	except:
		q1, q2, q3, iq = None, None, None, None

	if len(vals) < 2:
		return vals0, None, None, None, None, None, None, None, None, None, None, None, None, None

	med = np.median(vals).round(ndig)
	mu  = np.mean(vals).round(ndig)
	std = np.std(vals).round(ndig)
	n   = len(vals)
	SEM  = (std / np.sqrt(n)).round(ndig)
	mini = np.min(vals).round(ndig)
	maxi = np.max(vals).round(ndig)

	return vals, med, mu, std, n, SEM, mini, maxi, q1, q2, q3, iq, qsup, qinf


def calc_params_outliers_new(vals:List, remove_outliers:bool=True, ndig:int=4):
	vals, med, mu, std, n, SEM, mini, maxi, q1, q2, q3, iq, qsup, qinf = \
	calc_params_outliers(vals=vals, remove_outliers=remove_outliers, ndig=ndig)

	if mu is None or mu == 0:
		VC = None
	else:
		VC = std / mu

	return vals, med, mu, std, VC, n, SEM, mini, maxi, q1, q2, q3, iq, qsup, qinf


def calc_params_all_values(vals:List, ndig:int=4):

	vals = [x for x in vals if not pd.isnull(x)]

	if len(vals) < 2:
		return vals, None, None, None, None, None, None, None, None, None, None, None, None, None

	try:
		q1, q2, q3 = np.percentile(vals, [25,50,75]).round(ndig)
		iq = q3 - q1
		qsup = q3 + 1.5*iq
		qinf = q1 - 1.5*iq
	except:
		q1, q2, q3, iq = None, None, None, None

	med = np.median(vals).round(ndig)
	mu  = np.mean(vals).round(ndig)
	std = np.std(vals).round(ndig)
	n   = len(vals)
	SEM = (std / np.sqrt(n)).round(ndig)
	mini = np.min(vals).round(ndig)
	maxi = np.max(vals).round(ndig)

	if mu == 0:
		VC = None
	else:
		VC = std / mu

	return vals, med, mu, std, VC, n, SEM, mini, maxi, q1, q2, q3, iq, qsup, qinf


def clean_df_outliers(dfq, test):

	vals, med, mu, std, n, SEM, mini, maxi, q1, q2, q3, iq, qsup, qinf = calc_params_outliers(dfq[test])

	return dfq[ [True if x in vals else False for x in dfq[test]] ]

def filter_outlier(df, exam, group_id):
	df2 = df[df.group_id == group_id]
	vals = [x for x in df2[exam] if not pd.isnull(x)]

	vals, med, mu, std, n, SEM, mini, maxi, q1, q2, q3, iq, qsup, qinf = calc_params_outliers(vals)
	df2 = df[(df.group_id == group_id) & (pd.notna(df[exam])) & (df[exam] >= qinf) & (df[exam] <= qsup)].copy()
	# print(len(df2))
	return df2

def filter_outlier_groups(df, exam, groups):
	dfa = None
	for grp in range(len(groups)):
		df2 = filter_outlier(df, exam, group_id=grp)
		if dfa is None:
			dfa = df2
		else:
			dfa = dfa.append(df2)

	return dfa

'''
conda install -y -c "r/label/archive" r

conda install -y -c r r
conda install -y -c "r/label/archive" r

conda uninstall r

conda activate main
conda search r-base

conda install -c conda-forge r-base=4.2.0

R
install.packages('DescTools')

dunnett_test.R

library(DescTools)

file = './tmp/dunnett.tsv'
df = read.table(file, header = TRUE, sep = "\t")

control = df[1,'group']

dn <- DunnettTest(df$exam, df$group)

# print(dn[control])
dn <- dn[control]

file = './tmp/dunnet_result.tsv'
write.table(dn, file, sep='\t')

print("Ok dunnet")
#-------------------

dfn = get_dunnett_test(dfgroup, 'val')
'''
def get_dunnett_test_old(df, exam):

	dfa = df[[exam, 'group']].copy()
	dfa = dfa[dfa[exam].notna()]
	dfa.columns = ['exam', 'group']

	try:
		os.mkdir('./tmp/')
	except:
		pass

	pdwritecsv(dfa, "./tmp/dunnett.tsv")

	# install.packages('DescTools')
	if os.path.exists('./tmp/dunnett_result.tsv'):
		os.remove('./tmp/dunnett_result.tsv')

	os.system("Rscript ../src/dunnett_test.R ./tmp/dunnett.tsv 0.8")

	if os.path.exists('./tmp/dunnet_result.tsv'):
		dfdn = pdreadcsv('./tmp/dunnet_result.tsv')
	else:
		dfdn = None

	return dfdn


def get_dunnett_test(df:pd.DataFrame, exam:str, confidence:float=0.95):

	dfa = df[[exam, 'group']].copy()
	dfa = dfa[dfa[exam].notna()]
	dfa.columns = ['exam', 'group']

	try:
		os.mkdir('./tmp/')
	except:
		pass

	pdwritecsv(dfa, "./tmp/dunnett.tsv")

	# install.packages('DescTools')
	if os.path.exists('./tmp/dunnett_result.tsv'):
		os.remove('./tmp/dunnett_result.tsv')

	# new command 2024/08/26
	cmd = f"Rscript ../src/dunnett_test.R ./tmp/dunnett.tsv {confidence}"
	os.system(cmd)

	if os.path.exists('./tmp/dunnett_result.tsv'):
		dfdn = pdreadcsv('./tmp/dunnett_result.tsv')
	else:
		dfdn = None

	return dfdn


def scipy_dunnet(samples:List, control:List, alternative:str='two-sided'):
	'''
		dunnet test:
			assume normality
			samples against control
		
		samples: array of lists
		control: list
	'''
	return stats.dunnett(*samples, control=control, alternative=alternative)

def scipy_dunn(samples:List, p_adjust:str='holm') -> pd.DataFrame:
	df_pval = sposthocs.posthoc_dunn(samples, p_adjust=p_adjust)
	return df_pval

def scipy_dunn_ref(samples:List, ref_col:int=0, p_adjust:str='holm') -> list:
	df_pval = sposthocs.posthoc_dunn(samples, p_adjust=p_adjust)
	return df_pval.iloc[0].to_list()


def scipy_kruskal(*samples):
	return stats.kruskal(*samples)

def scipy_kruskal_mat(mat):
	return stats.kruskal(*mat)


def fdr_df(df, col):
	# df --> col is pval
	df2 = df.copy()
	if 'id' in df2.columns:
		df2.id = np.arange(0, df2.shape[0])
	else:
		df2['id'] = df2.index

	df2 = df2[['id', col]]

	df2.columns = ['id', 'pval']
	df2 = df2.sort_values('pval', ascending=True)

	pvals = df2[ [(not x is None) for x in df2.pval] ].pval
	nones = df2[ [(	x is None) for x in df2.pval] ].pval

	ranked_p_values = stats.rankdata(pvals)
	fdr = pvals * len(pvals) / ranked_p_values
	fdr[fdr > 1] = 1
	pvalsadj = list(fdr) + list(nones)

	dfr = pd.DataFrame( data={'fdr': pvalsadj, 'id': df2.id} )

	return dfr.sort_values('id').fdr

def fdr(p_vals):
	#--- cannot send p_vals with Nones
	#--- fixing it - to review
	p_vals = [1 if pd.isnull(x) else x for x in p_vals]
	p_vals = np.array(p_vals)

	ranked_p_values = stats.rankdata(p_vals)
	fdr = p_vals * len(p_vals) / ranked_p_values
	fdr[fdr > 1] = 1

	return fdr


# statsmodel result --> to DataFrame
# https://www.statsmodels.org/stable/generated/statsmodels.stats.multicomp.pairwise_tukeyhsd.html
# https://www.statsmodels.org/stable/generated/statsmodels.sandbox.stats.multicomp.TukeyHSDResults.html#statsmodels.sandbox.stats.multicomp.TukeyHSDResults
def results_summary_to_dataframe(results):
	'''take the result of an statsmodel results table and transforms it into a dataframe'''
	mat = []
	groups = results.groupsunique
	for i in range(len(groups)-1):
		for j in range(i+1, len(groups)):
			mat.append((groups[i], groups[j]))

	mat = np.array(mat)
	group1 = np.array(mat)[:, 0]
	group2 = np.array(mat)[:, 1]

	pvals = results.pvalues
	meandiffs = results.meandiffs
	conf_lower = np.array(results.confint)[:, 0]
	conf_higher = np.array(results.confint)[:, 1]
	reject = results.reject

	results_df = pd.DataFrame({"group1": group1,
							   "group2": group2,
							   "pvalue_adj":pvals,
							   "meandiffs":meandiffs,
							   "conf_lower":conf_lower,
							   "conf_higher":conf_higher,
							   "reject": reject })

	#Reordering...
	results_df = results_df[["group1", "group2", "meandiffs","pvalue_adj","conf_lower","conf_higher"]]
	results_df['abs(meandiffs)'] = [np.abs(x) for x in results_df.meandiffs]
	results_df = results_df.sort_values(["pvalue_adj", 'abs(meandiffs)'], ascending=[True, False])
	return results_df

def poisson(k, lamb):
	return (lamb**k/factorial(k)) * np.exp(-lamb)


def propag_errors_mult_div_ref(result:float, mu_val:float, std_val:float, n_val:int, 
											 mu_ref:float, std_ref:float, n_ref:int, alpha2:float=0.025):

	# result = result * np.sqrt( std_val*std_val/mu_val + std_ref*std_ref/mu_ref)

	std = result * np.sqrt( (std_val/mu_val)**2 + (std_ref/mu_ref)**2 )
	
	SEM = std / np.sqrt(n_val)
	# supposing normal distributed
	error = stats.t.ppf(1-alpha2, n_val-1)*SEM
	
	return SEM, error, std

def calc_normal_dist_errors(vals:List, std_list:List, is_ref_fold_change:bool=False, alpha2:float=0.025, 
							sample_limit:int=2, strict_normality:bool=False):

	vals = [x for x in vals if not pd.isnull(x)]
	n = len(vals)

	if n < sample_limit:
		return None, None, None, n, None, None, None, None, None, None, None

	mu  = np.mean(vals)
	med = np.median(vals)

	std_list = [x for x in std_list if not pd.isnull(x)]

	# sum of variances
	total_var = 0
	for i in range(n):
		total_var += std_list[i]**2

	# sqrt of <total variance>
	std = np.sqrt(total_var/n)
	SEM = std / np.sqrt(n)

	normal, stat_normal, pval_normal, stri_stat_normal = test_normality(vals)

	if is_ref_fold_change: normal = True


	if not strict_normality:
		try:
			error = stats.t.ppf(1-alpha2, n-1)*SEM
		except:
			error = None
	else:
		if normal == True:
			# if normal distributed
			try:
				error = stats.t.ppf(1-alpha2, n-1)*SEM
			except:
				error = None
		else:
			error = None
	
	VC = std/mu if mu != 0 else None


	return mu, med, std, n, SEM, error, VC, normal, stat_normal, pval_normal, stri_stat_normal



def propag_errors_mult_div_list(mu_list:List, std_list:List, alpha2:float=0.025):

	mu_list  = list(mu_list)
	std_list = list(std_list)

	mu  = np.mean(mu_list)
	med = np.median(mu_list)
	n   = len(mu_list)

	total = 0
	for i in range(n):
		total += (std_list[i]/mu_list[i])**2

	std = mu * np.sqrt(total)
	SEM = std / np.sqrt(n)
	error = stats.t.ppf(1-alpha2, n-1)*SEM
	
	return mu, med, std, n, SEM, error

def stat_dunnett_test(df:pd.DataFrame, col_list:List, col_control:str):

	vals0 = df[col_control]

	if isinstance(col_list, str):
		col_list = [col_list]

	if len(col_list) == 1:
		res = stats.dunnett(df[col_list[0]], control=vals0)

	elif len(col_list) == 2:
		res = stats.dunnett(df[col_list[0]], df[col_list[1]], control=vals0)

	elif len(col_list) == 3:
		res = stats.dunnett(df[col_list[0]], df[col_list[1]], df[col_list[2]], control=vals0)

	else:
		print("Open stat_lib and review stat_dunnett_test()")
		res = None


	return res

def stat_dunn_test(data:List, cols:List, p_adjust:str='bonferroni') -> pd.DataFrame:
	res = sposthocs.posthoc_dunn(data, p_adjust=p_adjust)

	dfd = pd.DataFrame(res)
	dfd.index = cols
	dfd.columns = cols

	return dfd

