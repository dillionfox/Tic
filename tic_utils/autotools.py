import numpy as np

def auto_filter(dG,opts):
	#---Compute autocorrelation for energy function and determine if result is 'periodic enough'
	ac = corr(dG,dG)
	if any(ac) == False or min(ac) > opts['ac_cutoff']:
		return False
	return True

def corr(x,y):
	#---Simple autocorrelation function
	try:
		result = np.correlate(x,y,mode='full')
	except ValueError:
		return [False]
	return result[int(result.size/2):]/max(result)
