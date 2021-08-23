# This module has some Jacobian-related methods
import numpy as np

def jacobian_update(J0, x0, x1, ovec0, ovec1):
    "updates the Jacobian using Broyden's method"
    deltao = ovec1-ovec0
    deltax = x1-x0
    deltaJ = np.outer((deltao - np.dot(J0,deltax)), deltax.reshape(1,4))/np.dot(deltax,deltax)
    J1 = J0 + deltaJ
    return J1

def export_jac(J):
    "Exports the jacobian to a text file"
    N = J.shape[0]; M = J.shape[1]
    f = open('Jacobian', 'w')
    for i in range(N):
        for j in range(M):
            f.write('%13.6e ' % J[i,j])
        f.write('\n')
    f.close()
    return