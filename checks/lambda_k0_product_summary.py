"""LAMBDA_K0_PRODUCT_SUMMARY  Summary of checks/lambda_k0_product_recovery.m (RESULTS.md section 89).
Recovered = within +-20% of the true value; true lambda*k0 = 60 * lambda.
Usage: py -3 checks/lambda_k0_product_summary.py"""
import pandas as pd, numpy as np
T = pd.read_csv('analysis_outputs/lambda_k0_product/recovery.csv')
T['true'] = 60 * T.lambdaTrue
ok = lambda est, tru: (np.abs(est / tru - 1) <= 0.2)
T['lamOK'] = ok(T.lambdaWS, T.lambdaTrue)
T['prodOK'] = ok(T.productWS, T['true'])
T['meanOK'] = ok(T.meanBurst, T['true'])
T['k0OK'] = ok(T.k0WS, 60)
print('Recovered within 20% (share of fits), by SNR:')
print(T.groupby('snr')[['lamOK', 'k0OK', 'prodOK', 'meanOK']].mean().round(2).to_string())
print()
print('By SNR and true lambda (lambda / product / mean burst):')
print(T.groupby(['snr', 'lambdaTrue'])[['lamOK', 'prodOK', 'meanOK']].mean().round(2).unstack(0).to_string())
print()
g = T.groupby(['snr', 'lambdaTrue'])
print('Mean estimate / true:')
print(pd.DataFrame({'lambda': g.lambdaWS.mean() / g.lambdaTrue.mean(), 'product': g.productWS.mean() / g['true'].mean(), 'meanBurst': g.meanBurst.mean() / g['true'].mean()}).round(2).unstack(0).to_string())
