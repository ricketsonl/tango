import numpy as np
import XGCfluxModel as XfM
import matplotlib.pyplot as plt

print('Running flux model')

dx = 0.01
fluxModel = XfM.dummyXGCFluxModel(dx)
x = np.linspace(0.,1.,num=100)
n = 1. + 0.5*np.sin(2.*np.pi*x)
profiles = {'n': n}
flux = fluxModel.get_flux(profiles)

print('Flux I got is:')
print(flux)

plt.figure(1)
plt.plot(x,profiles['n'],linewidth=2)
plt.plot(x,flux['n'],linewidth=2)
plt.savefig('myfig.png')
print('Figure saved')
plt.show()
