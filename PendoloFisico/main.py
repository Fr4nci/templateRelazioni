import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
buco = {"buco_1": [15.67, 15.80, 15.84, 15.98, 16.00, 16.03, 16.12], 
"buco_2": [15.59, 15.73, 15.57, 15.47, 15.54, 15.74, 15.71],
"buco_3": [15.85, 15.93, 15.96, 16.03, 15.93, 15.77, 15.93],
"buco_4": [18.65, 18.65, 18.58, 18.51, 18.73, 18.47, 18.52],
"buco_5": [38.20, 39.40, 42.49, 38.69, 38.01, 37.44, 39.27],
"buco_6": [22.77, 22.80, 22.67, 22.77, 22.50, 22.54, 22.34],
"buco_7": [16.70, 16.93, 16.65, 16.67, 16.70, 16.69, 16.53],
"buco_8": [15.56, 15.55, 15.63, 15.53, 15.57, 15.44, 15.89], 
"buco_9": [15.61, 15.72, 15.83, 15.89, 15.67, 15.81, 15.73],
"buco_10": [16.49, 16.42, 16.43, 16.45, 16.46, 16.37, 16.43]
}
dist = np.array([10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0], dtype=float)
dist = 105.6 - dist
l_cm = 54.55
d = abs(l_cm - dist)
err_dist = np.full(shape=len(dist), fill_value= np.sqrt( (0.1) ** 2 + (0.1) ** 2 ))
err_dist = err_dist * 10 ** (-2)
T = np.array([], dtype=float)
err_T = np.array([], dtype=float)
sigma = 0
media = 0
for el in buco:
    media = np.mean(buco[el])
    T = np.append(T, media)
    for i in buco[el]:
        sigma = sigma + (i-media)**2
    sigma = np.sqrt((1./(len(el) * (len(el)-1)))*sigma)
    err_T = np.append(err_T, sigma)
    media = 0
    sigma = 0
T = T/10
err_T = err_T/10
d = d * (10 ** (-2))
g = 9.80665

def period_model(d, l):
    return 2.0 * np.pi * np.sqrt((l**2.0 / 12.0 + d**2.0) / (g * d))

def derivative(d, l):
    return abs((np.pi * ( (1./12.) * l**2 - d**2))/(d**2 * g * np.sqrt((d**2 + (1./12.) * l**2)/(d * g))))
err_eff = np.array(err_T)
for i in range(0, 100):
    print(f"{i} iterazione")
    print(err_eff)
    popt, pcov = curve_fit(period_model, d, T, sigma=err_eff)
    for el in range(len(err_eff)):
        err_eff[el] = np.sqrt(err_T[el] ** 2 + (derivative(d[el], popt[0]) * err_dist[el]) ** 2 )
# ...and calculate the residuals with respect to the best-fit model.
res = T - period_model(d, *popt)
 # Create the main figure...
fig = plt.figure("Grafico di best-fit e residui")
print(f"l: {popt[0]} +- {np.sqrt(np.diag(pcov))}")
# ...and make space for the two plots. Note that ‘gridspec_kw‘ and ‘hspace‘
# control the arrangements of the two sub-panels within the figure, see
# https://matplotlib.org/stable/api/_as_gen/matplotlib.gridspec.GridSpec.html
ax1, ax2 = fig.subplots(2, 1, sharex=True, gridspec_kw=dict(height_ratios=[2, 1], hspace=0.05))
# Main plot: the scatter plot of x vs. y, on the top panel.
ax1.errorbar(d, T, err_T, fmt="o", label="Dati")
# Plot the best-fit model on a dense grid.
xgrid = np.linspace(0.0, 0.5, 100)
ax1.plot(xgrid, period_model(xgrid, *popt), label="Modello di best-fit")
# Setup the axes, grids and legend.
ax1.set_ylabel("T [s]")
ax1.grid(color="lightgray", ls="dashed")
ax1.legend()
# And now the residual plot, on the bottom panel.
ax2.errorbar(d, res, err_T, fmt="o")
# This will draw a horizontal line at y=0, which is the equivalent of the best-fit
# model in the residual representation.
ax2.plot(xgrid, np.full(xgrid.shape, 0.0))
# Setup the axes, grids and legend.
ax2.set_xlabel("d [m]")
ax2.set_ylabel("Residui [s]")
ax2.grid(color="lightgray", ls="dashed")
# The final touch to main canvas :-)
plt.xlim(0.0, 0.5)
fig.align_ylabels((ax1, ax2))
plt.savefig("Fit_e_residui.pdf")
plt.show()
chi2 = 0
for el in zip(T, err_T, d):
    chi2 = chi2 + ((el[0] - period_model(el[2], popt[0]))**2)/(el[1] ** 2)
print(chi2)
