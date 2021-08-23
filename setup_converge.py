# This module has methods that distribute the user input options and initialise the algorithm
import numpy as np
from scipy.stats import linregress
from argparse import ArgumentParser
from moognmodels import run_moog, run_model
from parsemoog import parse_moog_out, parse_moog_fe, diff_ovec

def initial_options():
    "Returns info from command-line arguments"
    parser = ArgumentParser()
    parser.add_argument('p0', help='First guess of atmospheric parameters', nargs='+', type=float)
    parser.add_argument('-p', '--plot', help='Plot iteration history.', action='store_true')
    parser.add_argument('-a', '--alpha', action='store_true', help='Correct for alpha abundances using the Salaris formula. Disabled for differential analysis.')
    parser.add_argument('-d', '--differential', action='store_true', help='Differential analysis. Reference star input file must be named reference.moog.')
    args = parser.parse_args()
    
    makeplot  = args.plot
    argsalpha = args.alpha
    return args.p0, makeplot, argsalpha, args.differential

def cmatrix():
    "Initial guess of the jacobian"
    
    c = np.array([[-2.41357143e-04, -2.53571429e-02, -1.19523810e-02, 5.75714286e-02],
                          [ 1.56907143e-03, -3.54476190e-01, -1.75142857e-01, 6.59285714e-02],
                          [ 2.29500000e-04,  1.87166667e-01, -8.90809524e-01,-4.04142857e-01],
                          [ 3.14428571e-04, -7.40476190e-03, -5.82380952e-02,-9.94571429e-01]])
    
    return c

def setup_atmpar(p0):
    "Configure iterations 0 and 1 of the atmospheric parameters"
    if p0 != None:
        if len(p0) == 4:
            p_new = np.array(p0)
        else:
            print('Wrong number of parameters!')
            raise SystemExit()
    else:
        p_new = np.array([4300., 1.60, -0.55, 1.72])
    
    p1 = p_new + np.array([1.,0.001,0.001,0.001])
    return p_new, p1

def std_anal_init(p0, argsalpha, ovec_ref=np.zeros(4)):
    "Initialisation for a non-differential (standard) analysis"
    C = cmatrix()
    Cinv = np.linalg.inv(C)
    
    p_new, p1 = setup_atmpar(p0)
    
    Jn, ovec = jacobian_initial_update(C, p_new, p1, argsalpha, ovec_ref)
    pn = p1
    quadr = np.linalg.norm(ovec)
    return p1, Jn, ovec, pn, quadr

def jacobian_initial_update(J0, x0, x1, args_alpha, ovec_ref):
    "Creates the first iteration of the Jacobian (standard analysis)"
    run_model(x0)
    run_moog()
    ovec0 = np.asarray(parse_moog_out(x0[2], alpha=args_alpha, native=True)) - ovec_ref
    run_model(x1)
    run_moog()
    ovec1 = np.asarray(parse_moog_out(x1[2], alpha=args_alpha, native=True)) - ovec_ref
    deltao = ovec1-ovec0
    deltax = x1-x0
    deltaJ = np.outer((deltao - np.dot(J0,deltax)), deltax.reshape(1,4))/np.dot(deltax,deltax)
    J1 = J0 + deltaJ
    return J1, ovec1

def dif_anal_init(p0, fe1ref, fe2ref, refmh):
    "Initialisation for a differential analysis"
    C = cmatrix()
    Cinv = np.linalg.inv(C)
    
    p_new, p1 = setup_atmpar(p0)
    
    Jn, ovec = jacobian_initial_differential_update(C, p_new, p1, fe1ref, fe2ref, refmh)
    pn = p1
    quadr = np.linalg.norm(ovec)
    return p1, Jn, ovec, pn, quadr

def jacobian_initial_differential_update(J0, x0, x1, fe1ref, fe2ref, refmh, stdout='a.out'):
    "Creates the first iteration of the Jacobian (differential analysis)"
    run_model(x0)
    run_moog()
    ovec0 = np.asarray(diff_ovec(stdout, fe1ref, fe2ref, x0[2], refmh))
    run_model(x1)
    run_moog()
    ovec1 = np.asarray(diff_ovec(stdout, fe1ref, fe2ref, x1[2], refmh))
    deltao = ovec1-ovec0
    deltax = x1-x0
    deltaJ = np.outer((deltao - np.dot(J0,deltax)), deltax.reshape(1,4))/np.dot(deltax,deltax)
    J1 = J0 + deltaJ
    return J1, ovec1
