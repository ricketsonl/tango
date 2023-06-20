import numpy as np
import adios2

N = 20
x = np.linspace(0.,1.,num=N)
print(x)

hi = "Hello there"

with adios2.open('dumbadiosout.bp', 'w') as fh:
    fh.write('array',hi)


with adios2.open('dumbadiosout.bp','r') as fh:
    y = fh.read('array')
    print(y)
