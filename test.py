
# import NMRLineshape
from nmr.lib import (
    EulerAngles,
    EulerAnglesInRadians,
    wigner_D_matrix,
    wigner_d_element,
    wigner_D_element,
)

from nmr.methods import (
    CSA_spinning_sideband,
    CSA_static_lineshape
)

import numpy as np
import matplotlib.pyplot as plt

a = EulerAngles('2.3 rad','2 rad','1 rad')
b = EulerAngles('2.53 rad','2 rad','1 rad')
print(a)
print(a.alpha)
print(a.beta)
print(a.gamma)
print (a+b)

a = EulerAnglesInRadians(2.3, 2 , 1)
b = EulerAnglesInRadians(2.53, 2, 1)
print(a)
print(a.alpha)
print(a.beta)
print(a.gamma)
print (a+b)

array = wigner_D_matrix(2, omega=EulerAnglesInRadians(2.,4.,0.))
print (array)
# print(array.ravel())

print(wigner_d_element(2, -2, -1, 0.5))
print(wigner_D_element(2, -2, -1, b))


# a = PowderScheme('lebedev', 38)
# print (a.weights)

# freq, data = CSA_spinning_sideband(
#                     start_frequency = -64.0,
#                     frequency_bandwidth = 128.0,
#                     number_of_points = 128,

#                     iso = 10,
#                     aniso = 10,
#                     eta = 0.0,

#                     ph_step = 30,
#                     spin_angular_freq = 1*2*np.pi,
#                     rotor_angle = (34.7356/180.)*np.pi,

#                     omega_PM = [0,0,0],

#                     averageing_scheme=1,
#                     averagin_size=6)

# plt.plot(freq, data)
# plt.show()
fig, ax = plt.subplots(3,1, sharex=True)
freq, data, time = CSA_spinning_sideband(
                    start_frequency = -512.0,
                    frequency_bandwidth = 1024.0,
                    number_of_points = 1024,

                    iso = 0.,
                    aniso = 120.,
                    eta = 0.7,

                    ph_step = 256,
                    spin_frequency = 10.,
                    rotor_angle = np.arccos(1/np.sqrt(3)),# (54.7356/180.)*np.pi,

                    omega_PM = np.asarray([0.,0.,0.]),

                    averaging_scheme=2,
                    averaging_size=95)

print(f'time {time} s')
ax[0].plot(freq, data)

freq, data, time = CSA_spinning_sideband(
                    start_frequency = -512.0,
                    frequency_bandwidth = 1024.0,
                    number_of_points = 1024,

                    iso = 0.,
                    aniso = 160.,
                    eta = 0.3,

                    ph_step = 1024*2,
                    spin_frequency =5.,
                    rotor_angle = np.arccos(1/np.sqrt(3)), #= (54.7356/180.)*np.pi,

                    omega_PM = np.asarray([0.,0.,0.]),

                    averaging_scheme=2,
                    averaging_size=95)

print(f'time {time} s')
ax[1].scatter(freq, data, s=0.2)


freq, data = CSA_static_lineshape(
                    (0*2*np.pi,
                     0.0,
                     160,
                     0.3),
                    start_frequency = -512.0,
                    frequency_bandwidth = 1024.0,
                    number_of_points = 1024,
)


ax[2].plot(freq, data)
plt.show()

# freq, data = CSA_spinning_sideband(
#                     start_frequency = -64.0,
#                     frequency_bandwidth = 128.0,
#                     number_of_points = 128,

#                     iso = 0,
#                     aniso = 10,
#                     eta = 0.0,

#                     ph_step = 30,
#                     spin_angular_freq = 1*2*np.pi,
#                     rotor_angle = (54.7356/180.)*np.pi,

#                     omega_PM = [0,0,0],

#                     averageing_scheme=1,
#                     averagin_size=1730)


# plt.plot(freq, data)
# plt.show()