# from tensor.convensions import Haeberlen
import numpy as np

class Frame:

    __slots__ = (
        '_cartesian',
        '_spherical',
        'tensor'
    ) 

    def __init__(self, tensor):
        self.tensor = tensor

    def _print_tensor(self, frame='cartesian'):
        self._cartesian = {
            'xx': 1,
            'xy': 1,
            'xz': 1,
            'yx': 1,
            'yy': 1,
            'yz': 1,
            'zx': 1,
            'zy': 1,
            'zz': 1
        }

        self._spherical = {
            'R00': 1,
            'R1+1': 1,
            'R10': 1,
            'R1-1': 1,
            'R2-2': 1,
            'R2-1': 1,
            'R20': 1,
            'R21': 1,
            'R22': 1,
         }

    def to(self, frame='cartesian'):
        """
        Return a dictionary with tensor elements in the given frame.

        The valid frames are 'cartesian' and 'spherical'.
        """
        if frame == 'cartesian':
            return self._cartesian
        if frame == 'spherical':
            return self._spherical
        raise KeyError(
            ( f"{frame} is not a valid frame. The valid frames "
               "are 'cartesian' and 'spherical'")
        )
    

class MRSite:

    __slots__ = (
        '_spin',
        '_isotope',
        '_nuclear_shielding_tensor',
        '_electric_quadrupole_tensor',
        '_R0_n',
        '_R1_n',
        '_R2_n',
        '_R0_e',
        '_R1_e',
        '_R2_e'
    )
    def __init__(
            self,
            isotope='1H',
            nuclear_shielding_tensor=[],
            electric_quadrupole_tensor=[]
            ):
        

        self._spin = self.get_spin_quantum_number(isotope)
        self._isotope = isotope
        self._nuclear_shielding_tensor = nuclear_shielding_tensor
        self._electric_quadrupole_tensor = electric_quadrupole_tensor

    @property
    def isotope(self):
        """Return the isotope name of the nuclear spin."""
        return self._isotope

    @property
    def spin(self):
        """Return the spin quantum number of the nuclear spin."""
        return self._spin
    
    @property
    def nuclear_shielding_tensor(self):
        """Return the nuclear shielding tensor."""
        return self._nuclear_shielding_tensor

    @nuclear_shielding_tensor.setter
    def nuclear_shielding_tensor(self, value):
        self._nuclear_shielding_tensor = value

    
    @property
    def electric_quadrupole_tensor(self):
        """Return the electric quadrupole tensor."""
        return self._electric_quadrupole_tensor

    @electric_quadrupole_tensor.setter
    def electric_quadrupole_tensor(self, value):
        """Return the electric quadrupole tensor."""
        self._electric_quadrupole_tensor = value

    
    def _update(self):
        pass

    def _get_PAS_values_for_nuclear_shielding(self, iso, aniso, eta):
        self._R0_n = np.asarray([iso])
        # not required.
        self._R1_n = np.asarray([])
        temp = -(aniso*eta)/np.sqrt(6)
        self._R2_n = np.asarray([temp, 0.0, aniso, 0.0, temp])

    def _get_PAS_values_for_quandrupole_coupling(self, Cq, eta):
        two_i = 2.0*self._spin
        wq = 6*np.pi*Cq/(two_i*(two_i-1))

        self._R0_e = np.asarray([0])
        # not required.
        self._R1_e = np.asarray([])
        temp = -eta/6
        self._R2_e = np.asarray([temp, 0.0, 1.0/np.sqrt(6), 0.0, temp])*wq

    def set_nuclear_shielding_tensor(self):
        pass
