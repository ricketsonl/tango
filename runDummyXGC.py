import numpy as np
import sys
sys.path.append("XGCcouplingTools")
import XGCfluxModel as XfM

dx = 1./499.

dummyXGC = XfM.dummyXGC(dx,imode='text',omode='adios')
dummyXGC.simulateXGClooping()
