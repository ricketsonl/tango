import numpy as np
import XGCfluxModel as XfM
import matplotlib.pyplot as plt

print('Running flux model')
num = 20
L = 1.
dx = L/num

fluxModel = XfM.XGCFluxModel(dx,imode='text',omode='adios')
x = np.linspace(0.,1.,num=num)
n = 1. + 0.5*np.sin(2.*np.pi*x)
profiles = {'n': n}
#flux = fluxModel.get_flux(profiles)
fluxModel.writeProfiles(profiles)

dumXGC = XfM.dummyXGC(dx,imode='text',omode='adios')
dumXGC.readProfiles()

readProf = dumXGC.profiles

print(readProf['n'])
print(profiles['n'])
print(readProf['n'] - profiles['n'])
'''print('Flux I got is:')
print(flux)

plt.figure(1)
plt.plot(x,profiles['n'],linewidth=2)
plt.plot(x,flux['n'],linewidth=2)
plt.savefig('myfig.png')
print('Figure saved')
plt.show()'''
