import numpy as np
from scipy import interpolate

class InteriorBCHandler(object):
    def __init__(self, x, fields, profiles, params):
        self.x = x
        self.fields = fields
        self.profiles = profiles
        self.interior_BC_locs = params['interior_BC_locs']
        self.interior_BC_desired_vals = params['interior_BC_vals']
        self.current_actual_BCs = params['actual_BC_vals']
        self.relaxation_params = params['relaxation']

        self.current_interior_vals = {}
        self.getCurrentInteriorVals
        

    def getCurrentInteriorVals(self):
        for field in self.fields:
           interpolated_profile = interpolate.interp1d(self.x, self.profiles[field.label])
           self.current_interior_vals[field.label] = interpolated_profile(self.interior_BC_locs[field.label])


    def updateActualBCs(self):
        self.getCurrentInteriorVals()
        updatedBCs = {}
        for field in self.fields:
            interior_error = self.interior_BC_desired_vals[field.label] - self.current_interior_vals[field.label]
            self.current_actual_BCs[field.label] += self.relaxation_params[field.label]*interior_error

    def checkConvergence(self,tol):
        max_error = 0.
        for field in self.fields:
            interior_error = self.interior_BC_desired_vals[field.label] - self.current_interior_vals[field.label]
            max_error = max(max_error, np.fabs(interior_error))

        if max_error < tol:
            return True
        else:
            return False
        
