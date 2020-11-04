import numpy as np
from itertools import product
import math


def rotation_matrix(axis, angle):
	a = math.cos(angle/2.0)
	s = math.sin(angle/2.0)
	b = axis[0]*s
	c = axis[1]*s
	d = axis[2]*s

	return np.array([[a*a + b*b - c*c - d*d, 2*(b*c - a*d), 2*(b*d+a*c)],
					 [2*(b*c+a*d), a*a+c*c-b*b-d*d, 2*(c*d-a*b)],
					 [2*(b*d-a*c), 2*(c*d+a*b), a*a+d*d-b*b-c*c]])



lattice_parameter = 4.046

basis = np.array([[1.0, 0.0, 0.0],
                  [0.0, 1.0, 0.0],
                  [0.0, 0.0, 1.0]])*lattice_parameter

base_atoms = np.array([[0.0, 0.0, 0.0],
                       [0.5, 0.5, 0.0],
                       [0.5, 0.0, 0.5],
                       [0.0, 0.5, 0.5]])*lattice_parameter

axis = [1, 0, 0]
angle = np.pi/6
rot = rotation_matrix(axis,angle)

basis = np.matmul(rot.T, basis)
base_atoms = np.array([np.matmul(rot,atom) for atom in base_atoms])


system_size = 20
sustembasis = np.array([[1.0, 0.0, 0.0],
                        [0.0, 1.0, 0.0],
                        [0.0, 0.0, 1.0]])*system_size

invbasis = np.linalg.inv(basis)
invsystem = np.linalg.inv(systembasis)
corners = np.array(list(product([0.0, system_size], repeat=3)))
extentes = np.array([np.matmul(invbasis.T, corner) for corner in corners])

mins = [int(np.min([m,0])) for m in np.floor(np.min(extents, axis=0))]
maxs = [int(np.min([m,0]))+1 for m in np.floor(np.max(extents, axis=0))]

positions = []
for i in range(mins[0], maxs[0]):
    for j in range(mins[1], maxs[1]):
        for k in range(mins[2], maxs[2]):
            base_position = np.array([i,j,k])
            cart_position = np.inner(basis.T, base_position)
            for atom in base_atoms:
                atompos = cart_position + atom
                systempos = np.matmul(invsystem.T, atompos)
                if all(systempos>=0.0) and all(systempos<1.0):
                    positions.append(atompos)

with open('rotatedcrystalline.data', 'w') as fdata:
    fdata.write('Random atoms - written for Asker \n\n')

    fdata.write('{} atoms\n'.format(len(positions)))
    fdata.write('{} atom types\n'.format(1))

    fdata.write('{} {} xlo xhi\n'.format(0.0, system_size))
    fdata.write('{} {} ylo yhi\n'.format(0.0, system_size))
    fdata.write('{} {} zlo zhi\n'.format(0.0, system_size))
    fdata.write('\n')

    fdata.write('Atoms\n\n')

    for i,pos in enumerate(positions):
        fdata.write('{} 1 {} {} {}\n'.format(i+1,*pos))
