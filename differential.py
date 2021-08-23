# This module contains methods related to differential analysis
import os
import numpy as np
import matplotlib.pyplot as plt
from moognmodels import run_moog, run_model
from parsemoog import parse_moog_fe, diff_ovec
from setup_converge import dif_anal_init
from jacobian import jacobian_update
from converge_plots import make_dif_feplot

def reference_star(finref='ref.par', finrefatm='refatm', stdout='ref.out'):
    "Returns the observables from the reference star"
    "This method needs existing finref and finrefatm files to work"
    "finref: MOOG input file for reference star"
    "finrefatm: ASCII file containing the reference atmospheric parameters. Must be written as a single-line text file in the format 'Teff logg metal micro'"
    "stdout: MOOG summary_out for reference star"
    
    # backing-up the MOOG input file
    os.system('cp batch.par backupforrefrun.par')
    os.system('cp %s batch.par' % finref)
    # runs reference star
    refatm = np.genfromtxt(finrefatm)
    run_model(refatm) # teff logg [M/H] ([Fe/H] here) micro
    run_moog()
    # restores par file
    os.system('cp batch.par %s' % finref)
    os.system('mv backupforrefrun.par batch.par')
    # loads Fe data
    fe1data, fe2data = parse_moog_fe(stdout)
    return fe1data[:,6], fe2data[:,6], refatm[2]

def diff_main(par0, makeplot, stdout='a.out'):
    "The main routine for differential analysis"
    "par0 (list): Initial guess of the atmospheric parameters given as command-line argument"
    "makeplot (boolean): command-line option. If true, iteration sequence is plotted"
    
    # loads reference star info
    fe1ref, fe2ref, refmh = reference_star()
    
    # initialisation of Broyden's method
    p1, Jn, ovec, pn, quadr = dif_anal_init(par0, fe1ref, fe2ref, refmh)
    
    # === using Broyden's method to converge the atmospheric parameters ===
    old_quadr = 6e23
    counter = 0
    # creating matplotlib object if the option of plotting the iterations is requested
    if makeplot == True:
        fig, ax = plt.subplots(2,2)
    while quadr > 0.001 and abs(old_quadr - quadr) > 1e-7 and counter < 20:
        # updates the atmospheric parameters vector
        Jinv = np.linalg.inv(Jn)
        delta_p = -np.dot(Jinv, ovec)
        pn = pn + delta_p
        old_ovec = ovec
        
        # external stuff: creates an input model for MOOG and runs MOOG
        run_model(pn)
        run_moog()
        
        # loading results from MOOG
        ovec = np.asarray(diff_ovec(stdout, fe1ref, fe2ref, pn[2], refmh))
        
        # the minimization parameter
        old_quadr = quadr
        quadr = np.linalg.norm(ovec)
        
        # this will be useful if convergence is not yet achieved
        counter = counter + 1
        if p1.shape[0] == 4:
            pn_1 = pn-delta_p
            Jn = jacobian_update(Jn, pn_1, pn, old_ovec, ovec)
        
        # print info about the current iteration
        print('%7.1e %2i' % (quadr, counter))
        print(delta_p)
        
        # populate the plots if requested
        if makeplot == True:
            ax[0,0].scatter(counter, ovec[0])
            ax[0,1].scatter(counter, ovec[1])
            ax[1,0].scatter(counter, ovec[2])
            ax[1,1].scatter(counter, ovec[3])
            ax[0,0].hlines(0.0, 0.0, counter, linestyles='dotted')
            ax[0,1].hlines(0.0, 0.0, counter, linestyles='dotted')
            ax[1,0].hlines(0.0, 0.0, counter, linestyles='dotted')
            ax[1,1].hlines(0.0, 0.0, counter, linestyles='dotted')
    
    # Creates (differential) Boltzmann plots    
    make_dif_feplot(fe1ref, fe2ref, pn[2], refmh)
    return pn, quadr, Jn
