# This module has a method for calculating the uncertainties
import numpy as np
from scipy.stats import linregress
from parsemoog import parse_moog_fe
from differential import reference_star
from setup_converge import cmatrix

def make_errors(C, fin='a.out', diff_analysis=False, usecmatrix=True):
    "Calculates the internal uncertainties of the atmospheric parameters"
    Cinv = np.linalg.inv(C)
    fe1, fe2 = parse_moog_fe(fin)
    if diff_analysis == False:
        a1 = fe1[:,6]
        a2 = fe2[:,6]
    else:
        fe1_ref, fe2_ref, ref_metal = reference_star()
        a1 = fe1[:,6] - fe1_ref
        a2 = fe2[:,6] - fe2_ref
    ep = fe1[:,2]
    rw = fe1[:,5]
    deltaion = []
    mm = []
    for i in range(1000):
        stoc = np.random.normal(np.mean(a1), np.std(a1, ddof=1), a1.shape[0])
        sto2 = np.random.normal(np.mean(a2), np.std(a2, ddof=1), a2.shape[0])
        deltaion.append((np.mean(a1+stoc)-np.mean(a2+sto2)))
        mm.append((np.mean(stoc-a1) + np.mean(sto2-a2))*0.5)
    ovec = -np.array([linregress(ep, a1)[-1], np.std(deltaion, ddof=1), np.std(mm, ddof=1), linregress(rw, a1)[-1]])
    dp = np.dot(Cinv, ovec)
    if np.any(dp<=0) or usecmatrix==True:
        jac = cmatrix()
        Cinv = np.linalg.inv(jac)
        dp = np.dot(Cinv, ovec)
    else:
        pass
    return dp