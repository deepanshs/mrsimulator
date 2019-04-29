import numpy as np
from collections import namedtuple
from dataclasses import dataclass
from .unit import string_to_quantity
from astropy.units import Quantityv

@dataclass
class Hamiltonian:
    _first_order_nuclear_shielding_: int = 0
    _first_order_electric_quadrupole_coupling_: int =0
    _second_order_electric_quadrupole_coupling_: int=0

    def __add__(self, other):
        i_ = self._first_order_nuclear_shielding_
        i_+= other._first_order_nuclear_shielding_
        j_ = self._first_order_electric_quadrupole_coupling_
        j_+= other._first_order_electric_quadrupole_coupling_
        k_ = self._second_order_electric_quadrupole_coupling_
        k_+= other._second_order_electric_quadrupole_coupling_
        return Hamiltonian(i_, j_, k_)

class FirstOrderNuclearShielding:
    def __new__(self):
        return Hamiltonian(_first_order_nuclear_shielding_ = 1)

class FirstOrderElectricQuadrupoleCoupling:
    def __new__(self):
        return Hamiltonian(_first_order_electric_quadrupole_coupling_ = 1)

class SecondOrderElectricQuadrupoleCoupling:
    def __new__(self):
        return Hamiltonian(_second_order_electric_quadrupole_coupling_ = 1)


def _get_haeberlen_values_from_PAS(xx, yy, zz):
    unit = xx.unit
    lst = np.asarray([xx.value, yy.value, zz.value])
    iso = lst.mean()
    lst-=iso
    i = np.argmax(np.abs(lst))
    zeta = lst[i]
    np.delete(lst, i)
    eta = np.abs((lst[0]-lst[1])/zeta)
    return HaeberlenConvension(iso*unit, zeta*unit, eta)


def _get_PAS_from_haeberlen_values(iso, zeta, eta):
    zz = zeta + iso
    temp1 = 0.5*eta*zeta
    temp2 = 1.5*iso - 0.5*zz
    yy = temp2 + temp1
    xx = temp2 - temp1
    return PrincipalAxisSystem(xx, yy, zz)


@dataclass(frozen=True)
class HaeberlenConvension:
    iso: str = 0.0*u.Hz
    zeta: str = 0.0*u.Hz
    eta: str = 0.0

    def __post_init__(self):
        iso = self.iso
        if not isinstance(iso, Quantity):
            iso = string_to_quantity(str(iso))
        iso = iso.to('Hz')
        
        zeta = self.zeta
        if not isinstance(zeta, Quantity):
            zeta = string_to_quantity(str(zeta))
        zeta = zeta.to('Hz')
        
        eta = float(self.eta)
        
        super(HaeberlenConvension, self).__setattr__('eta', eta)
        super(HaeberlenConvension, self).__setattr__('zeta', zeta)
        super(HaeberlenConvension, self).__setattr__('iso', iso)
        
    @property
    def principal_values(self):
        return _get_PAS_from_haeberlen_values(*astuple(self))


@dataclass(frozen=True)
class PrincipalAxisSystem:
    xx: str = 0.0*u.Hz
    yy: str = 0.0*u.Hz
    zz: str = 0.0*u.Hz
        
    def __post_init__(self):
        xx = self.xx
        if not isinstance(xx, Quantity):
            xx = string_to_quantity(str(xx))
        xx = xx.to('Hz')
        
        yy = self.yy
        if not isinstance(yy, Quantity):
            yy = string_to_quantity(str(yy))
        yy = yy.to('Hz')
        
        zz = self.zz
        if not isinstance(zz, Quantity):
            zz = string_to_quantity(str(zz))
        zz = zz.to('Hz')
        
        super(PrincipalAxisSystem, self).__setattr__('xx', xx)
        super(PrincipalAxisSystem, self).__setattr__('yy', yy)
        super(PrincipalAxisSystem, self).__setattr__('zz', zz)
    
    @property
    def haeberlen_values(self):
        return _get_haeberlen_values_from_PAS(*astuple(self))


@dataclass(frozen=True)
class EulerAngle:
    alpha: str = 0.0*u.rad
    beta: str = 0.0*u.rad
    gamma: str = 0.0*u.rad
        
    def __post_init__(self):
        alpha = self.alpha
        if not isinstance(alpha, Quantity):
            alpha = string_to_quantity(str(alpha))
        alpha = alpha.to('rad')
        
        beta = self.beta
        if not isinstance(beta, Quantity):
            beta = string_to_quantity(str(beta))
        beta = beta.to('rad')
        
        gamma = self.gamma
        if not isinstance(gamma, Quantity):
            gamma = string_to_quantity(str(gamma))
        gamma = gamma.to('rad')
        
        super(EulerAngle, self).__setattr__('alpha', alpha)
        super(EulerAngle, self).__setattr__('beta', beta)
        super(EulerAngle, self).__setattr__('gamma', gamma)
        


class ParameterizedTensor:
    @dataclass
    class QuadupoleTensor:
        Cq: str = '0.0 Hz'
        eta: float = 0.0
        PAS: PrincipalAxisSystem = PrincipalAxisSystem('0 Hz', '0 Hz', 0) 
        euler_angle: EulerAngle = EulerAngle('0.0 rad','0.0 rad','0.0 rad')
        
        def __post_init__(self):
            self.eta = float(self.eta)
            if not isinstance(self.eta, float):
                raise Exception(f'Expecting a float instance for eta, {type(self.eta)} found.')
            if not isinstance(self.Cq, str) and not isinstance(self.Cq, Quantity):
                raise Exception(f'Expecting a str or Quantity instance for Cq, {type(self.Cq)} found.')
            else:
                self.Cq = string_to_quantity(self.Cq)
            if not isinstance(self.euler_angle, EulerAngle):
                self.euler_angle = EulerAngle(*self.euler_angle)

    class NuclearShieldingTensor:
        PAS: PrincipalAxisSystem = PrincipalAxisSystem('0 Hz', '0 Hz', 0)
        euler_angle: EulerAngle = EulerAngle('0.0 rad','0.0 rad','0.0 rad')
        pass

class spin(Tensor):
    __slots__=(
        '_nuclear_shielding_tensor',
        '_electric_quadrupole_tensor'
    )

    def __init__(
            self, 
            isotope='1H',
            _nuclear_shielding_tensor = None,
            _electric_quadrupole_tensor = None,
            *args,
            **kwargs):


        self._isotope = isotope

        # if _nuclear_shielding_tensor 
        # self._nuclear_shielding_tensor






    def __set__tensor(self, value):
        if isinstance(vlaue, )