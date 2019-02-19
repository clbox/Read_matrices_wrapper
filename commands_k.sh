f2py -m friction_read_matrices_H_S_p1_k -h friction_read_matrices_H_S_p1_k.pyf friction_read_matrices_H_S_p1_k.f90
f2py -c friction_read_matrices_H_S_p1_k.pyf friction_read_matrices_H_S_p1_k.f90 --fcompiler=intel --f90exec=mpiifort
