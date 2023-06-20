import numpy as np
import XGCfluxModel as XfM

print('Running dummy XGC')

dx = 0.01
XGC = XfM.dummyXGC(dx)
XGC.writeFlux()
