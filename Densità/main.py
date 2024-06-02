import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
plt.rcParams.update({
    'font.size': 8,
    'text.usetex': True,
    'text.latex.preamble': r'\usepackage{amsfonts} \\ \usepackage{siunitx}'
})
def mod_d(x, m, q):
    return(m*x + q)


print(np.pi)
class Parallelepipedo():
    def __init__(self, materiale, m, sigma_m, base, larghezza, altezza, sigma_b, sigma_l, sigma_h):
        self.materiale = materiale
        self.m = m
        self.sigma_m = sigma_m
        self.b = base
        self.l = larghezza
        self.h = altezza
        self.sigma_b = sigma_b
        self.sigma_l = sigma_l
        self.sigma_h = sigma_h        
        self.calc_volume()
    def get_mass(self):
        print(f"La massa: {self.m} +- {self.sigma_m}")
    def calc_volume(self):
        self.volume = self.b * self.l * self.h
        self.err = np.sqrt(np.power(self.sigma_b/self.b, 2) + np.power(self.sigma_l/self.l, 2) + np.power(self.sigma_h/self.h, 2)) * self.volume
    def get_volume(self):
        print(f"{self.volume} +- {self.err}")

class Sfera():
    def __init__(self, materiale, m, sigma_m, d, sigma_d):
        self.materiale = materiale
        self.m = m
        self.sigma_m = sigma_m
        self.d = d
        self.sigma_d = sigma_d
    def get_mass(self):
        print(f"{self.m} +- {self.sigma_m}")
    def calc_volume(self):
        self.volume = (1./6.) * np.pi * np.power(self.d, 3)
        self.err = 3 * (self.sigma_d/self.d) * self.volume
    def get_volume(self):
        print(f"Il volume: {self.volume} +- {self.err}")

class Cilindro():
    def __init__(self, materiale, m, sigma_m, d, sigma_d, h, sigma_h):
        self.materiale = materiale
        self.m = m
        self.sigma_m = sigma_m
        self.d = d
        self.sigma_d = sigma_d
        self.h = h
        self.sigma_h = sigma_h
    def get_mass(self):
        print(f"{self.m} +- {self.sigma_m}")
    def calc_volume(self):
        self.volume = np.pi * np.power(self.d, 2) * (self.h / 4)
        self.err = np.sqrt(4 * np.power(self.sigma_d/self.d, 2) + np.power(self.sigma_h/self.h)) * self.volume
    def get_volume(self):
        print(f"Il volume: {self.volume} +- {self.err}")

p1 = Parallelepipedo("alluminio", 0.250, 0.001, 0.20, 0.10, 0.30, 0.01, 0.01, 0.01)

print("Informazioni sul primo parallelepipedo")
p1.get_mass()
p1.get_volume()

# Dati alluminio
alluminio_v = [p1.volume, p2.volume, p3.volume]
alluminio_sigma_v = [p1.err, p2.err, p3.err]
alluminio_m = [p1.m, p2.m, p3.m]
alluminio_sigma_m = [p1.sigma_m, p2.sigma_m, p3.sigma_m]

# Grafico e modello di best-fit alluminio

popt, pcov = curve_fit(model, alluminio_m, alluminio_v, [1./2710., 0], sigma=alluminio_sigma_v)
res = alluminio_v - mod_d(alluminio_m, *popt)

fig = plt.figure("Un grafico dei residui")

ax1, ax2 = fig.subplots(2, 1, sharex=True, gridspec_kw=dict(height_ratios=[2, 1], hspace=0.05))

ax1.errorbar(alluminio_m, alluminio_v, alluminio_sigma_v, fmt='o', label='Dati')
xgrid = np.linspace(0.0, 10.0, 100)
ax1.plot(xgrid, mod_d(xgrid, *popt), label='Modello di best-fit')
ax1.set_ylabel(r"$V$ $[\si{\meter^3}]$")
ax1.grid(color="lightgray", ls="dashed")
ax1.legend()


ax2.errorbar(alluminio_m, res, alluminio_sigma_v, fmt='o')
ax2.plot(xgrid, np.full(xgrid.shape, 0.0))
ax2.set_xlabel("$m$ $[\si{\kilo\gram}]$")
ax2.set_ylabel("Residui $[\si{\meter^3}]$")
plt.xlim(0.0, 10.0)
fig.align_ylabels((ax1, ax2))

plt.show()


