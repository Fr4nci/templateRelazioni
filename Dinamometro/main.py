import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

dati_T = {
    "m_1": np.array([5.65, 5.56, 5.77, 5.73, 5.85, 5.86, 5.81, 5.73, 5.77, 5.77], dtype=np.float64),
    "m_2": np.array([7.59, 7.50, 7.53, 7.53, 7.57, 7.48, 7.46, 7.55, 7.45, 7.46], dtype=np.float64),
    "m_3": np.array([6.21, 6.25, 6.27, 6.22, 6.29, 6.29, 6.25, 6.27, 6.20, 6.25], dtype=np.float64),
    "m_4": np.array([4.71, 4.77, 4.74, 4.73, 4.86, 4.77, 4.52, 4.91, 4.66, 4.75], dtype=np.float64),
    "m_5": np.array([8.83, 8.87, 8.82, 8.98, 8.83, 8.69, 8.82, 8.98, 8.88, 8.76], dtype=np.float64)
}
m = np.array([10.652, 28.945, 16.541, 4.694, 34.956], dtype=np.float64)
m = m * (10 ** (-3))
m_p = 7.874 * (10 ** (-3))
m_m = 8.093 * (10 ** (-3))
m = m + m_p + (m_m/3)
sigma_m = np.full(len(m), dtype=np.float64, fill_value=np.sqrt((0.001 * 10** (-3))**2 + (0.001 * 10 ** (-3))**2 + ((0.001/3) * 10 ** (-3))**2))
T = []
sigma_T = []
for el in dati_T:
    T.append((np.mean(dati_T[el]))/9.25)
    sigma_T.append(np.std(dati_T[el])/9.25)
T = np.array(T, dtype=np.float64)
sigma_T = np.array(sigma_T, dtype=np.float64)

for i in range(len(T)):
    print(f"Per la massa m_{i+1} il periodo è: {T[i]} +- {sigma_T[i]}")
def modello_periodo(x, k):
    return(2 * np.pi * np.sqrt(x/k))

popt, pcov = curve_fit(modello_periodo, m, T, sigma=sigma_T)
# sigma_eff = alluminio_sigma_v
# for i in range(10):
#    sigma_eff = np.sqrt(np.power(alluminio_sigma_v, 2) + np.power((1./popt[0]) * alluminio_sigma_m, 2))
#    popt, pcov = curve_fit(mod_d, alluminio_m, alluminio_v, [1./2.710], sigma=sigma_eff)
#    chisq = np.sum(((alluminio_v - mod_d(alluminio_m, *popt))/sigma_eff) ** 2)
#    print(f"Step {i}")
#    print(popt, np.sqrt(pcov.diagonal()))
#    print(f"Chisquare = {chisq:.2f}")


res = T - modello_periodo(m, *popt)

fig = plt.figure("Un grafico dei residui")

ax1, ax2 = fig.subplots(2, 1, sharex=True, gridspec_kw=dict(height_ratios=[2, 1], hspace=0.05))

ax1.errorbar(m, T, sigma_T, fmt='o', label='Dati')
xgrid = np.linspace(0.0, 0.2, 100)
ax1.plot(xgrid, modello_periodo(xgrid, *popt), label='Modello di best-fit')
ax1.set_ylabel("T [s]")
ax1.grid(color="lightgray", ls="dashed")
ax1.legend()


ax2.errorbar(m, res, sigma_T, fmt='o')
ax2.plot(xgrid, np.full(xgrid.shape, 0.0))
ax2.set_xlabel("m [g]")
ax2.set_ylabel("Residui [s]")
plt.xlim(0.0, 0.2)
fig.align_ylabels((ax1, ax2))
plt.savefig("Modello_best_fit_k_residui.pdf")
chisq = np.sum((res/(sigma_T))**2)
print(f"{chisq:.2f}")
plt.show()

print(f"{popt[0]} +- {np.sqrt(pcov.diagonal()[0])}")
k = popt[0]
sigma_k = np.sqrt(pcov.diagonal())
# Parte due dell'esperienza
m = np.array([10.652, 28.945, 16.541, 4.694, 34.956], dtype=np.float64)
m = m * (10 ** (-3))
m_p = 7.874 * (10 ** (-3))
m = m + m_p
sigma_m = np.full(len(m), dtype=np.float64, fill_value=np.sqrt((0.001 * 10** (-3))**2 + (0.001 * 10 ** (-3))**2))
def modello(x, c):
    return(c * np.ones_like(x))

l0 = 11.3 * (10 ** (-2))
sigma_l0 = 1 * (10 ** (-2))
lf = np.array([19.3, 27.10, 21.5, 16.6, 29.3], dtype=np.float64)
lf = lf * (10 ** (-2))
delta_l = lf - l0
sigma_l_f = np.full(len(lf), fill_value=0.1 * (10 ** (-2)))
sigma_delta_l = np.sqrt(sigma_l0**2 + sigma_l_f ** 2)
for i in range(len(delta_l)):
    print(f"Massa {i+1}: {delta_l[i]} +- {sigma_delta_l[i]}")
g = (k/m) * (lf - l0)
sigma_g = g * np.sqrt((sigma_k/k) ** 2 + (sigma_m/m)**2 + (sigma_delta_l/(delta_l))**2)
for i in range(len(m)):
    print(f"{sigma_g[i]} +- {(k/m[i]) * delta_l[i] * sigma_m[i]}")
popt, pcov = curve_fit(modello, m, g, sigma=sigma_g,p0=9.81)
res = g - modello(m, *popt)

fig = plt.figure("Un grafico dei residui")

ax1, ax2 = fig.subplots(2, 1, sharex=True, gridspec_kw=dict(height_ratios=[2, 1], hspace=0.05))

ax1.errorbar(m, g, fmt='o', label='Dati', yerr=sigma_g)
xgrid = np.linspace(0.0, 0.2, 100)
ax1.plot(xgrid, modello(xgrid, *popt), label='Modello di best-fit')
ax1.set_ylabel("g [m/s^2]")
ax1.grid(color="lightgray", ls="dashed")
ax1.legend()


ax2.errorbar(m, res, yerr=sigma_g, fmt='o')
ax2.plot(xgrid, np.full(xgrid.shape, 0.0))
ax2.set_xlabel("m [g]")
ax2.set_ylabel("Residui [m/s^2]")
plt.xlim(0.0, 0.2)
fig.align_ylabels((ax1, ax2))
plt.savefig("Modello_best_fit_residui_g.pdf")
plt.show()
print(f"{popt[0]} +- {np.sqrt(pcov.diagonal())}")

