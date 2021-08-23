import numpy as np
import matplotlib.pyplot as plt
from parsemoog import parse_moog_out, parse_moog_fe
from moognmodels import run_moog, run_model
from jacobian import jacobian_update, export_jac

def converge_main(p1, Jn, ovec, pn, quadr, makeplot, argsalpha):
    "Uses Broyden's method to converge the atmospheric parameters"
    
    old_quadr = 6e23
    counter = 0
    # creating matplotlib object if requested
    if makeplot == True:
        fig, ax = plt.subplots(2,2)
    while quadr > 0.001 and abs(old_quadr - quadr) > 1e-7 and counter < 20:
        # updates the atmospheric parameters vector
        Jinv = np.linalg.inv(Jn)
        delta_p = -np.dot(Jinv, ovec)
        old_ovec = ovec
        pn = pn + delta_p
        
        # external stuff: creates an input model for MOOG and runs MOOG
        run_model(pn)
        run_moog()
        
        # loading results from MOOG
        ovec = parse_moog_out(pn[2], alpha=argsalpha)
        ovec = np.asarray(ovec)
        old_quadr = quadr
        
        # the minimization parameter
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
    return pn, quadr, Jn

def final_remarks(pn, quadr, Jn, makeplot, argsalpha, analtype, datm):
    "The results..."
    
    final_ovec = parse_moog_out(pn[2], alpha=argsalpha, final=True)
    
    # print print print...
    print('\nFINAL RESULTS:')
    if analtype == True:
        print('Method: differential')
    else:
        print('Method: standard')
    print('-- Output parameters quadrature: %7.1e' % quadr)
    print('-- Output parameters: %6.3f %6.3f %6.3f %6.3f' % (final_ovec[0],final_ovec[1],final_ovec[2],final_ovec[3]))
    print('-- Atmospheric parameters: %6.1f %6.3f %6.3f %7.4f' % (pn[0],pn[1],pn[2],pn[3]))
    print('-- Atmospheric uncertainties: %6.1f %6.3f %6.3f %7.4f' % (datm[0],datm[1],datm[2],datm[3]))
    print('-- alpha: %s' % str(argsalpha))
    if argsalpha == True:
        print('-- [alpha/Fe] = %.3f\t[Fe/H] = %.3f' % (final_ovec[4], final_ovec[5]))
    
    # saving the atmospheric parameters to a text file
    f = open('atmparam', 'w')
    f.write('%4.1f %.3f %.3f %.4f ' % (pn[0],pn[1],pn[2],pn[3]))
    if argsalpha == True:
        f.write('%.3f %.3f\n' % (final_ovec[4], final_ovec[5]))
    else:
        f.write('\n')
    f.close()
    
    # saving the final Jacobian to a text file
    export_jac(Jn)
    
    # saving the uncertainties of the atmospheric parameters
    f = open('unc_atmparam', 'w')
    f.write('%4.1f %.3f %.3f %.4f ' % (datm[0],datm[1],datm[2],datm[3]))
    f.write('\n')
    f.close()
    
    # if plotting is requested
    if makeplot == True:
        plt.show()