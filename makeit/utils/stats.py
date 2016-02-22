from __future__ import print_function
import numpy as np

def mse(x, y):
	'''Mean-squared error of two 1D numpy arrays or Python lists'''
	x = np.array(x)
	y = np.array(y)
	return np.mean((x - y) ** 2)

def mae(x, y):
	'''Mean absolute error of two 1D numpy arrays or Python lists'''
	x = np.array(x)
	y = np.array(y)
	return np.mean(np.abs(x - y))

def q(y_true, y_pred):
	'''q value as described in Tropsha, Gramatica, Gombar:
	The Importance of Being Earnest'''
	y_true = np.array(y_true)
	y_pred = np.array(y_pred)
	y_mean = np.mean(y_true)
	return 1 - np.sum((y_true - y_pred) ** 2) / np.sum((y_true - y_mean) ** 2)

def linreg(x, y):
	'''Computes a linear regression through the origin using OLS'''
	x = np.array(x)
	y = np.array(y)
	x = x[:, np.newaxis]
	a, _, _, _ = np.linalg.lstsq(x, y)
	r2 = q(y, (x * a)[:, 0])
	return (r2, a)