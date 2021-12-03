"""QuantumTransmissionCoefficientCalculator

:platform: Python 2.7

This python code is to calculate quantum tunneling transmission coefficient
using piece wise constant method.

To use, first specify a barrier structure in myStructure.py
myStructure.py will be imported as below
Then, create an object, call compute method with numpy array specifying
desired energy range of interest. It will return transmission coefficient.

See main function in this module for usage
"""

from __future__ import division
import numpy as np
from scipy import constants as sc
from matplotlib import pyplot as plt
import importlib

# custom modules
import nu
import myutil
# A module containing a structure of interest
# module name can be chosen arbitrary but set it to 'as myStr'
desiredModule = input("Which module would you like to import?")
myStr = importlib.import_module(desiredModule)


class QuantumTransmissionCoefficientCalculator():
    """This class provides quantum transmission coefficient.

    Creates an object and call a function to get started.
    """
    def __init__(self, N=1):
        """Set up potential structured

        Set up barriers based on descriptions in myStructure.py.
        What it does is to create self.strcuture and self.biased_structure

        It has an ability to further subdivide the given structure.
        It is useful when applying a bias as it smoothes out the shape of the
        potential structure.

        Args:
            N: Number of subdivision to perform
        """
        self._set_up_structure(N)

    def _set_up_structure(self, N):
        U = np.array([])
        x = np.array([])
        m_e = np.array([])
        # set proper units
        myStr_t = np.array(myStr.thickness)*nu.nm
        myStr_m_e = np.array(myStr.effective_mass)*sc.m_e
        myStr_U = np.array(myStr.potential)*nu.eV
        myStr.position = np.array(myStr.position)*nu.nm
        # perform subdivision
        # thickness t -> t/N of N sub-regions
        x = [np.linspace(myStr_t[0:_].sum(), myStr_t[0:_+1].sum(),
             num=N, endpoint=False) for _ in range(myStr_t.size)]
        x = np.ravel(x)
        x = np.append(x, x[-N]+myStr_t[-1])
        x = np.delete(x, 0)
        m_e = np.append(m_e, [np.ones(N)*myStr_m_e[_]
                        for _ in range(myStr_m_e.size)])
        U = np.append(U, [np.ones(N)*myStr_U[_]
                      for _ in range(myStr_U.size)])
        # keep data in structured array
        self.structure = np.zeros(len(U), dtype=[('U', 'float'),
                                                 ('x', 'float'),
                                                 ('m_e', 'float')])
        self.structure['U'] = U
        self.structure['x'] = x
        self.structure['m_e'] = m_e
        # make a copy of the original in use in applying bias
        self.biased_structure = np.copy(self.structure)
        return

    def apply_bias(self, bias):
        """Apply bias to the given potential strcuture.

        Args:
            bias: Specify a value in (V)
        """
        x = self.structure['x']
        x = np.insert(x, 0, 0)
        # E field V/m
        E_field = bias / (x[-1]-x[0])
        self.biased_structure['U'] = (self.structure['U'] /
                                      nu.eV-E_field*x[1:])*nu.eV
        return

    def plot_structure(self):
        """Plot potential structure
        """
        U = self.biased_structure['U']
        x = self.biased_structure['x']
        # double up for plotting purpose
        x = np.ravel(np.dstack((x, x)))
        x = np.insert(x, 0, 0)
        x = np.delete(x, -1)
        U = np.ravel(np.dstack((U, U)))
        # set up max and min value for plotting purpose
        Vmax = np.max(U)
        Vmax = 1.05*Vmax/sc.e
        Vmin = np.min(U)
        Vmin = 1.05*Vmin/sc.e
        xmin = x[0]/nu.nm
        xmax = x[-1]/nu.nm
        # plot
        plt.plot(x/nu.nm, U/nu.eV)
        plt.grid()
        plt.xlabel('position (nm)')
        plt.ylabel('potential (eV)')
        plt.xlim(xmin, xmax)
        plt.ylim(Vmin, Vmax)
        plt.show()

    def compute(self, E):
        """Compute transmission coefficient in the energy range E

        To compute transmission coefficient, piece wise constant method is
        performed.

        Args:
            E: specify desired energy range in (J)
            e.g. E = np.linspace(0, 1.0)*1.60e-19

        Return:
            TC: transmission coefficient (no unit)
        """
        U = self.biased_structure['U']
        x = self.biased_structure['x']
        x = np.insert(x, 0, 0)
        m_e = self.biased_structure['m_e']
        k = np.array([np.sqrt(2*m_e*(E[_]-U+0j)) / sc.hbar
                     for _ in range(len(E))])
        cns_m_e = (m_e[1:]/m_e[:-1])
        cns_k = (k[:, :-1]/k[:, 1:])
        cns_m_e_k = cns_m_e*cns_k
        M11 = (0.5*(1+cns_m_e_k)*np.exp(-1j*(k[:, 1:]-k[:, :-1])*x[1:-1]))
        M12 = (0.5*(1-cns_m_e_k)*np.exp(-1j*(k[:, 1:]+k[:, :-1])*x[1:-1]))
        M21 = (0.5*(1-cns_m_e_k)*np.exp(1j*(k[:, 1:]+k[:, :-1])*x[1:-1]))
        M22 = (0.5*(1+cns_m_e_k)*np.exp(1j*(k[:, 1:]-k[:, :-1])*x[1:-1]))

        m11, m12, m21, m22 = M11[:, -1], M12[:, -1], M21[:, -1], M22[:, -1]
        for __ in range(len(U)-2):
            func = lambda m1, m2, m3, m4: m1*m2+m3*m4
            a = func(m11, M11[:, -__-2], m12, M21[:, -__-2])
            b = func(m11, M12[:, -__-2], m12, M22[:, -__-2])
            c = func(m21, M11[:, -__-2], m22, M21[:, -__-2])
            d = func(m21, M12[:, -__-2], m22, M22[:, -__-2])
            m11, m12, m21, m22 = a, b, c, d

        MT22 = m22
        ret = ((m_e[-1]/m_e[0])*(k[:, 0]/k[:, -1]) *
               (MT22*np.conjugate(MT22))**-1)
        TC = np.where(np.isnan(ret), 0, ret.real)
        return TC


def main():
    # set up output plot font size
    myutil.setup_fonts()
    # set up energy range of interest
    # Note that 0 eV gives a warning -> division by zero
    E = np.linspace(0.01, 1.0, 2000)*nu.eV
    # create an object RTD with N=10 subdivision
    RTD = QuantumTransmissionCoefficientCalculator(20)
    # check configured potential structure described in myStructure.py
    # which is imported at the beginnning of this module
    RTD.plot_structure()
    # compute method gives a transmission coefficient TC
    TC1 = RTD.compute(E)
    # apply bias .2 (V)
    RTD.apply_bias(0.2)
    RTD.plot_structure()
    # recalculate transmission coefficient
    TC2 = RTD.compute(E)
    # Results
    plt.plot(TC1, E/nu.eV, label='0 V')
    plt.plot(TC2, E/nu.eV, label='0.2 V')
    plt.xlabel('Transmission coefficient')
    plt.ylabel('Energy (eV)')
    plt.legend(loc='best')
    plt.grid()
    plt.show()
    return


if __name__ == '__main__':
    main()
