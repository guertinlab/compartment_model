# used `Nilson2017_Quantify_Kinit_TRP.ai` to roughly quantify the fold changes 
## the files to load is `kinit_Nilson_2017.txt`
import pandas as pd
import numpy as np
import scipy.optimize
import matplotlib.pyplot as plt

def monoExp(x, a, b, c):
        return a * np.exp(-b * x) + c
def genLog(x, base, k, m, b):
    return m * (base**(x*k)) + b

initial_guesses = [1.0, 0.001, 0.0]

kinit_factor_df = pd.read_csv("kinit_Nilson_2017.txt", sep='\t')
paramsSat, cv = scipy.optimize.curve_fit(monoExp,
                                         kinit_factor_df['time_sec'],
                                         kinit_factor_df['kinit_AU'],  [1.0, 0.001, 0.0])
a_fit, b_fit, c_fit = paramsSat

#paramsSat, cv = scipy.optimize.curve_fit(genLog, kinit_factor_df['time_sec'],
#                                         kinit_factor_df['kinit_AU'], [2.0, -1,1,0] )

## again curve fitting failing, so used Gemini
#y = 0.848⋅e ^ (−0.005x)+0.157
#dy/dx = - 0.0045869 . e ^ (−0.005x)
a = 0.848
b = 0.005
c = 0.157
y_fit = monoExp(kinit_factor_df['time_sec'], a_fit, b_fit, c_fit)
derivative_fit = monoExp(kinit_factor_df['time_sec'], -1 * a_fit * b_fit, 0.005, 0)
kinitFC_at_12andhalfmins = monoExp(750, a_fit, b_fit, c_fit)
slop_fit__at_12andhalfmins = monoExp(750, -1 * a_fit * b_fit, b_fit, 0)
kinit_factor_df['fittedPoints'] = y_fit
kinit_factor_df['derivative']= derivative_fit
kinit_factor_df.to_csv("Nilson2017_Kinit_fit_and_deriv.txt", sep="\t", index=False)

plt.clf()
plt.figure(figsize=(10, 6))
plt.scatter( kinit_factor_df['time_sec'],  kinit_factor_df['kinit_AU'], label='Original Data', color='blue', marker='o')
plt.plot(kinit_factor_df['time_sec'], y_fit, label=f'Fitted Curve: $y = {a_fit:.3f}e^{{-{b_fit:.3f}x}} + {c_fit:.3f}$', color='red')
plt.title('Non-linear Curve Fit (Exponential Decay)')
plt.xlabel('time_sec')
plt.ylabel('kinit (AU)')
plt.legend()
plt.grid(True)
plt.show()

## the ggplot2 version for plotting `Nilson2017_Kinit_fit_and_deriv.txt` is in `plotFitsForNilson2017.R`

