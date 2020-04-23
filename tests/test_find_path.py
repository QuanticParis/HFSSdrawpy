# -*- coding: utf-8 -*-
"""
Created on Thu Apr 12 23:30:42 2018

@author: checkhov
"""

import matplotlib.pyplot as plt
import numpy as np
from scripts.designer import ConnectElt, Vector, Circuit
plt.close('all')

def val(self, nb):
    return nb
ConnectElt.val = val
ConnectElt.to_bond = []


A_pos = Vector([0,0])
A_ori = Vector([0, 1])
B_pos = Vector([1, 1])
B_ori = Vector([0, 1])

c = Circuit()
c.key_elt('port_A', A_pos, A_ori)
c.key_elt('port_B', B_pos, B_ori)
c.port_A.create_port()
c.port_B.create_port()

cable = c.connect_elt('cable', 'port_A', 'port_B')


found_cable = cable.find_path(fillet=0.05, is_meander=True, to_meander=[0,2,0,0,0,0,0], meander_length=0.4, meander_offset=0.1)
to_bond_points = ConnectElt.to_bond

fig, ax = plt.subplots(figsize = (6,12))

start = [A_pos+ A_ori.orth()*0.05, A_pos- A_ori.orth()*0.05]
end = [B_pos+ B_ori.orth()*0.05, B_pos- B_ori.orth()*0.05]


ax.plot(np.array(start).T[0],np.array(start).T[1], color='red')
ax.plot(np.array(end).T[0],np.array(end).T[1], color='red')
ax.plot(np.array(found_cable).T[0],np.array(found_cable).T[1], 'o-')

for ii, points in enumerate(to_bond_points):
    ax.plot(np.array(points).T[0],np.array(points).T[1], color='r')

ax.set_xlim((-1,1))
ax.set_ylim((-0.5,1.5))
ax.axis('equal')
