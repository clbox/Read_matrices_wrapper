f2py -m friction_read_matrices_H_S_p1 -h friction_read_matrices_H_S_p1.pyf friction_read_matrices_H_S_p1.f90
f2py -c friction_read_matrices_H_S_p1.pyf friction_read_matrices_H_S_p1.f90 --fcompiler=intel --f90exec=mpiifort
