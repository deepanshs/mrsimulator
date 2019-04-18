
# coding: utf-8

# In[1]:


import NMRmethods as nmr
import numpy as np
import matplotlib.pyplot as plt
from pymatgen.analysis.nmr import ChemicalShielding


# In[2]:


print(nmr.CSA_static_lineshape.__doc__)


# Passing a symmetric tensor to the ChemicalShielding class

# In[3]:


Si29_tensor = [[443.2366, -10.9178, 0.0133  ],
               [9.7144  , 435.4115, 0.0330  ],
               [0.1192  , 0.6163  , 433.3699]]

Si29_tensor = np.asarray(Si29_tensor)
Si29_tensor = (Si29_tensor+Si29_tensor.T)/2

tensor = ChemicalShielding(Si29_tensor)


# Simulate the NMR shielding/chemical shift line shape spectrum using

# In[4]:


freq, amp = nmr.CSA_static_lineshape(tensor.haeberlen_values,
                        number_of_points=512,
                        start_frequency=420.,
                        frequency_bandwidth=50)


# In[5]:


plt.plot(freq,amp)
plt.xlabel('$^{29}$Si nuclear shielding frequency / ppm')
plt.ylabel('intensity')
plt.show()


# In[ ]:




