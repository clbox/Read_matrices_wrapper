import numpy as np
import numpy.f2py as f2py
import friction_read_matrices_H_S_p1
 

n_basis = 2094
n_spin = 1 
n_k_points = 13
n_cells_in_hamiltonian = 1
friction_n_active_atoms = 2
real_eigenvectors = True
friction_index_list = np.array([37,38])


S = np.empty((3, friction_n_active_atoms, n_k_points, n_basis*n_basis, n_spin))
H = np.empty((3, friction_n_active_atoms, n_k_points, n_basis*n_basis, 1, n_spin))
S,H,S_c,H_c = friction_read_matrices_H_S_p1.friction_read_matrices_h_s_p1(n_spin,n_basis,n_k_points,n_cells_in_hamiltonian,friction_n_active_atoms,real_eigenvectors,friction_index_list)
