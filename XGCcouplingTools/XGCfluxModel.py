"""See https://github.com/LLNL/tango for copyright and license information"""

from __future__ import division
import numpy as np
import os
import adios2
import os.path
import time
import re

"""Shestakov Test Module"""
class XGCFluxModel(object):
    """Class-based interface to the flux model in the Shestakov analytic
    test problem.
    
    Slight modification of boundary conditions so that n(1) = nL rather than n(1)=0,
    in order to make the computation of the flux at the boundary much easier by
    avoiding a divide by zero condition.
    """
    def __init__(self, dx, XGCpath=None, writePath=None, imode='npy', omode='npy'):
        self.dx = dx
        self.XGCpath = XGCpath
        self.writePath = writePath
        self.imode = imode # input mode
        self.omode = omode # output mode

        ## File extenstion for input (profiles)
        if imode == 'adios':
            input_ext = '.bp'
        elif imode == 'npy':
            input_ext = '.npy'
        else:
            input_ext = '.txt'

        ## File extension for output (fluxes)
        if omode == 'adios':
            output_ext = '.bp'
        elif omode == 'npy':
            output_ext = '.npy'
        else:
            output_ext = '.txt'

        if (writePath is None):
            self.profile_marker_file_name = 'profilesWritten.done'
            self.flux_marker_file_name = 'fluxWritten.done'
            self.profile_data_name = 'profiles' + input_ext
            self.flux_data_name = 'flux' + output_ext
        else:
            self.profile_marker_file_name = writePath + 'profilesWritten.done'
            self.flux_marker_file_name = writePath + 'fluxWritten.done'
            self.profile_data_name = writePath + 'profiles' + input_ext
            self.flux_data_name = writePath + 'flux' + output_ext

    def writeProfiles(self, profiles):
        nx = profiles['n'].shape[0]
        x_grid = np.linspace(0.,self.dx*(nx-1),num=nx)
        if self.imode == 'adios':
            ## Write adios profile files
            with adios2.open(self.profile_data_name, 'w') as fh:
                shape = []; start = []; count = [profiles['n'].shape[0]]
                fh.write('n',profiles['n'], shape, start, count)
        elif self.imode == 'npy':
            np.save(self.profile_data_name,profiles['n'])
        elif self.imode == 'text':
            with open(self.profile_data_name,'w') as f:
                f.write(str(nx) + '\n')
                for i in range(profiles['n'].shape[0]):
                    f. write(str(x_grid[i]) + '       ' + str(profiles['n'][i]))
                    f.write('\n')
                f.write('-1')
        else:
            print('MAYDAY: This profile writing mode not implemented yet')

        #print('Wrote profile')

        ## Write file indicating successful write of profiles
        f = open(self.profile_marker_file_name,'w')
        f.close()
    
    def get_flux(self, profiles,verbose=False):
        self.writeProfiles(profiles)
        ## Wait to see a marker file indicating flux has been written
        while not os.path.exists(self.flux_marker_file_name):
            #print('Waiting for flux to exist')
            time.sleep(0.1)
        if verbose:
            print('Found flux file')

        flux = {}
        if self.omode == 'adios':
            with adios2.open(self.flux_data_name,'r') as fh:
                flux['n'] = fh.read('flux_n')
        elif self.omode == 'npy':
            flux['n'] = np.load(self.flux_data_name)
        else:
            print('MAYDAY: This flux reading mode is not implemented yet')
        
        if os.path.exists(self.flux_marker_file_name):
            os.remove(self.flux_marker_file_name) # Having fetched flux, delete flux marker file

        return flux


class dummyXGC(object):
    def __init__(self, dx, TangoPath=None, writePath=None, imode='npy', omode='npy'):
        self.TangoPath=TangoPath
        self.writePath=None
        self.imode = imode # input mode
        self.omode = omode # output mode

        self.dx = dx

        ## File extenstion for input (profiles)
        if imode == 'adios':
            input_ext = '.bp'
        elif imode == 'npy':
            input_ext = '.npy'
        else:
            input_ext = '.txt'

        ## File extension for output (fluxes)
        if omode == 'adios':
            output_ext = '.bp'
        elif omode == 'npy':
            output_ext = '.npy'
        else:
            output_ext = '.txt'

        if (writePath is None):
            self.profile_marker_file_name = 'profilesWritten.done'
            self.flux_marker_file_name = 'fluxWritten.done'
            self.profile_data_name = 'profiles' + input_ext
            self.flux_data_name = 'flux' + output_ext
        else:
            self.profile_marker_file_name = writePath + 'profilesWritten.done'
            self.flux_marker_file_name = writePath + 'fluxWritten.done'
            self.profile_data_name = writePath + 'profiles' + input_ext
            self.flux_data_name = writePath + 'flux' + output_ext
        self.profiles = {}

    def readProfiles(self,termination_marker='converged.done'):
        ## Wait for Tango to notify you that profiles have been written
        while not os.path.exists(self.profile_marker_file_name) and not os.path.exists(termination_marker):
            time.sleep(1)

        if self.imode == 'adios':
            with adios2.open(self.profile_data_name,'r') as fh:
                self.profiles['n'] = fh.read('n')
        elif self.imode == 'npy':
            self.profiles['n'] = np.load(self.profile_data_name)
        elif self.imode == 'text':
            with open(self.profile_data_name,'r') as f:
                num = int(f.readline().split('\n')[0])
                self.profiles['n'] = np.zeros(num)
                for i in range(num):
                    myline = f.readline()
                    myval = re.split(' ', myline)[-1]
                    myval = myval.split('\n')[0]
                    self.profiles['n'][i] = float(myval)
                
        else:
            print('MAYDAY: Selected read mode not implemented')

    def writeFlux(self,termination_marker='converged.done',verbose=False):
        self.readProfiles(termination_marker=termination_marker)
        if verbose:
            print('Successfully read profiles')
        if os.path.exists(self.profile_marker_file_name):
            os.remove(self.profile_marker_file_name) # Having read profiles, delete the profile marker file
        n = self.profiles['n']
        flux = {}
        flux['n'] = get_flux(n, self.dx)

        if self.omode == 'adios':
            with adios2.open(self.flux_data_name,'w') as fh:
                shape = []; start = []; count = [flux['n'].shape[0]]
                fh.write('flux_n', flux['n'],shape,start,count)
        elif self.omode == 'npy':
            np.save(self.flux_data_name, flux['n'])
        else:
            print('MAYDAY: Selected write mode not implemented')

        ## Write file that indicates fluxes have been written
        f = open(self.flux_marker_file_name,'w')
        f.close()

    def simulateXGClooping(self,maxits=100,termination_marker_file='converged.done'):
        it = 0
        while it < maxits and not os.path.exists(termination_marker_file):
            self.writeFlux()
        



#class shestakov_analytic_fluxmodel(object):
#    """Class-based interface to the flux model in the Shestakov analytic
#    test problem.
#    
#    Slight modification of boundary conditions so that n(1) = nL rather than n(1)=0,
#    in order to make the computation of the flux at the boundary much easier by
#    avoiding a divide by zero condition.
#    
#    Alias to AnalyticFluxModel which conforms to style guidleines.  This class is kept around for backwards compatibility reasons.
#    """
#    def __init__(self, dx):
#        self.dx = dx
#    
#    def get_flux(self, n):
#        return get_flux(n, self.dx)
    
    
    
#==============================================================================
#   Functions specifying the Shestakov analytic test problem start here
#==============================================================================
def H7contrib_Source(x):
    S = GetSource(x)
    H7 = S
    return H7


def get_flux(n, dx):
    """Test problem from Shestakov et al. (2003)
    Return the flux Gamma, which depends on the density profile n as follows:
       Gamma[n] = -(dn/dx)^3 / n^2
    """
    Gamma = np.zeros_like(n)
    
    # Return flux Gamma on the same grid as n
    dndx = _dxCenteredDifference(n, dx)
    Gamma = - dndx**3 / n**2    
    return Gamma
    
def GetSource(x, S0=1, delta=0.1):
    """Test problem from Shestakov et al. (2003).
    Return the source S."""
    S = np.zeros_like(x)
    S[x < delta] = S0
    return S
    
def _dxCenteredDifference(u, dx):
    """Compute du/dx.
      du/dx is computed using centered differences on the same grid as u.  For the edge points, one-point differences are used.
    
    Inputs:
      u         profile (array)
      dx        grid spacing (scalar)
    
    Outputs:
      dudx      (array, same length as u)
    """
    dudx = np.zeros_like(u)
    dudx[0] = (u[1] - u[0]) / dx
    dudx[1:-1] = (u[2:] - u[:-2]) / (2*dx)
    dudx[-1] = (u[-1] - u[-2]) / dx
    return dudx
    
def steady_state_solution(x, nL, S0=1, delta=0.1):
    """Return the exast steady state solution for the Shestakov test problem
    
    Inputs:
      x             Spatial coordinate grid (array)
      nL            boundary condition n(L) (scalar)
      S0            parameter in source term --- amplitude (scalar)
      delta         parameter in source term --- location where it turns off (scalar)
    Outputs:
    """
    nright = ( nL**(1/3) + 1/3 * (S0 * delta)**(1/3) *(1-x) )**3
    # nleft = (L - delta + 0.75*(delta - x**(4/3) / delta**(1/3)))**3-
    nleft = ( nL**(1/3) + 1/3 * (S0 * delta)**(1/3) * (1 - delta + (3/4) * (delta - x**(4/3) / delta**(1/3))))**3
    nss = nright
    nss[x < delta] = nleft[x < delta]
    return nss
    
def GetSteadyStateSolution(x, nL, S0=1, delta=0.1):
    """alias to steady_state_solution which conforms to style guidelines, but this function is kept around for backwards compatibility."""
    return steady_state_solution(x, nL, S0=S0, delta=delta)
