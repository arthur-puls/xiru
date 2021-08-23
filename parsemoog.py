import numpy as np
from scipy.stats import linregress

def mini_parse(skip_h, skip_f, spid, spstr, fin):
    "Auxiliary method for parse_moog_out"
    for i in range(len(spid)):
        if spid[i] == spstr:
            skip_idx = i
            break
    if spstr == 'Fe I':
        sp = np.genfromtxt(fin,
                           usecols=(2,5,6),
                           skip_header=skip_h[skip_idx],
                           skip_footer=skip_f[skip_idx])
    else:
        try:
            sp = np.genfromtxt(fin,
                           usecols=6,
                           skip_header=skip_h[skip_idx],
                           skip_footer=skip_f[skip_idx])
        except:
            sp = np.zeros(0)
    return sp

def parse_moog_out(minput, feI=True, feII=True, alpha=True, native=True, final=False, fesolar=7.50, mgsolar=7.60, sisolar=7.51, casolar=6.34, fin='a.out'):
    "This method retrieves information from MOOG's summary_out file"
    "Returns the convergence parameters EP slope, RW slope, difference in FeI/II abundances and difference between model and estimated metallicity"
    "If alpha-enhancement is considered it returns alpha and Fe abundances separately"
    "Solar abundances from Asplund et al. 2009ARA&A..47..481A"
    "The 'fin' variable MUST coincide with the summary_out in MOOG input file"
    asol = np.array([mgsolar, sisolar, casolar])
    f = open(fin, 'r')
    t = f.readlines()
    f.close()
    if native == False:
        N = len(t)
        spnames = []
        avabd = []
        idxl = []
        for i in range(N):
            if 'Abundance Results for Species' in t[i]:
                spnames.append(t[i])
            if 'average abundance' in t[i]:
                avabd.append(t[i])
                idxl.append(i)
        if feI == True:
            for i in range(len(spnames)):
                ionid = spnames[i].split()[4]+' '+spnames[i].split()[5]
                if ionid == 'Fe I':
                    fe1 = float(avabd[i].split()[3])
                    epslope = float(t[idxl[i]+1].split()[4])
                    rwslope = float(t[idxl[i]+2].split()[4])
                    break
            feh = fe1 - fesolar
        else:
            return None, None, None, None
        if feII == True:
            for i in range(len(spnames)):
                ionid = spnames[i].split()[4]+' '+spnames[i].split()[5]
                if ionid == 'Fe II':
                    fe2 = float(avabd[i].split()[3])
                    break
            gparam = fe1 - fe2
        else:
            gparam = None
        if alpha == True:
            mg = None; si = None; ca = None
            for i in range(len(spnames)):
                ionid = spnames[i].split()[4]+' '+spnames[i].split()[5]
                if ionid == 'Mg I':
                    mg = float(avabd[i].split()[3])
                if ionid == 'Si I':
                    si = float(avabd[i].split()[3])
                if ionid == 'Ca I':
                    ca = float(avabd[i].split()[3])
                    break # needs alpha elements to be ordered by atomic number
            alphal = [mg, si, ca]
            try:
                ahl = alphal - asol
            except:
                for i in range(len(alphal)-1, -1, -1):
                    if alphal[i] == None:
                        del alphal[i]
                        asol = np.delete(asol, i)
                ahl = alphal - asol
            afel = ahl - feh
            afe = np.mean(afel)
            deltam = np.log10(0.638*10.**afe + 0.362)
            mparam = np.around(feh + deltam - minput, decimals=3)
        else:
            mparam = np.around(feh - minput, decimals=3)
    else:
        skip_h = []
        skip_f = []
        spid = []
        hcounter = 0
        fcounter = 0
        for i in range(len(t)):
            if 'Abundance Results for Species' in t[i]:
                hcounter = hcounter + 1
                skip_h.append(i+hcounter+1)
                ion = t[i].split()[4]+' '+t[i].split()[5]
                spid.append(ion)
            if 'average abundance' in t[i]:
                skip_f.append(len(t)-i)
        offset_h = np.arange(len(spid))
        offset_f = np.arange(len(spid)-1, -1, -1)
        skip_h = np.asarray(skip_h) - offset_h
        skip_f = np.asarray(skip_f) - offset_f
        if feI == True:
            fe1 = mini_parse(skip_h, skip_f, spid, 'Fe I', fin)
        if feII == True:
            fe2 = mini_parse(skip_h, skip_f, spid, 'Fe II', fin)
        if alpha == True:
            mg = mini_parse(skip_h, skip_f, spid, 'Mg I', fin)
            si = mini_parse(skip_h, skip_f, spid, 'Si I', fin)
            if si.shape != () and si.shape[0] == 0:
                si = mini_parse(skip_h, skip_f, spid, 'Si II', fin)
            ca = mini_parse(skip_h, skip_f, spid, 'Ca I', fin)
        epslope = np.around(linregress(fe1[:,0], fe1[:,2])[0], decimals=4)
        rwslope = np.around(linregress(fe1[:,1], fe1[:,2])[0], decimals=4)
        fe1_mean = np.mean(fe1[:,2])
        feh = fe1_mean - fesolar
        if feII == True:
            gparam = np.around(fe1_mean - np.mean(fe2), decimals=4)
        else:
            gparam = None
        if alpha == True:
            prealpha = [mg, si, ca]; alphal = []; alphadelarg = []
            for i in range(len(prealpha)):
                if prealpha[i].shape != () and prealpha[i].shape[0] == 0:
                    alphadelarg.append(i)
                else:
                    alphal.append(prealpha[i])
            if len(alphadelarg) > 0:
                asol = np.delete(asol, alphadelarg)
            for i in range(len(alphal)):
                try:
                    alphal[i] = np.mean(alphal[i])
                except:
                    pass
            try:
                ahl = alphal - asol
            except:
                for i in range(len(alphal)-1, -1, -1):
                    if alphal[i] == None:
                        del alphal[i]
                        asol = np.delete(asol, i)
                ahl = alphal - asol
            afel = ahl - feh
            afe = np.mean(afel)
            deltam = np.log10(0.684*10.**afe + (1.-0.684))
            mparam = np.around(feh + deltam - minput, decimals=4)
        else:
            mparam = np.around(feh - minput, decimals=4)
    if alpha == False or final == False:
        return epslope, gparam, mparam, rwslope
    else:
        return epslope, gparam, mparam, rwslope, afe, feh

def parse_moog_fe(fin):
    "Gets MOOG output for Fe lines only"
    f = open(fin, 'r')
    t = f.readlines()
    f.close()
    skip_h = []
    skip_f = []
    spid = []
    hcounter = 0
    for i in range(len(t)):
        if 'Abundance Results for Species' in t[i]:
            hcounter = hcounter + 1
            skip_h.append(i+hcounter+1)
            ion = t[i].split()[4]+' '+t[i].split()[5]
            spid.append(ion)
        if 'average abundance' in t[i]:
            skip_f.append(len(t)-i)
    offset_h = np.arange(len(spid))
    offset_f = np.arange(len(spid)-1, -1, -1)
    skip_h = np.asarray(skip_h) - offset_h
    skip_f = np.asarray(skip_f) - offset_f
    for i in range(len(spid)):
        if spid[i] == 'Fe I':
            fe1idx = i
        if spid[i] == 'Fe II':
            fe2idx = i
    fe1data = np.genfromtxt(fin, skip_header=skip_h[fe1idx], skip_footer=skip_f[fe1idx])
    fe2data = np.genfromtxt(fin, skip_header=skip_h[fe2idx], skip_footer=skip_f[fe2idx])
    return fe1data, fe2data

def diff_ovec(stdout, fe1ref, fe2ref, mh, refmh, mkplot=False):
    "Generates differential ovec, i.e., the vector of the minization parameters, but with respect to the reference star"
    fe1data, fe2data = parse_moog_fe(stdout)
    dfe1 = fe1data[:,6] - fe1ref
    dfe2 = fe2data[:,6] - fe2ref
    dfe1m = np.mean(dfe1)
    dfe2m = np.mean(dfe2)
    epsl = linregress(fe1data[:,2], dfe1)[0]
    rwsl = linregress(fe1data[:,5], dfe1)[0]
    gpar = dfe1m - dfe2m
    mpar = np.mean(fe1data[:,6]) - np.mean(fe1ref) - mh + refmh
    if mkplot == False:
        return epsl, gpar, mpar, rwsl
    else:
        epint = linregress(fe1data[:,2], dfe1)[1]
        rwint = linregress(fe1data[:,5], dfe1)[1]
        return epsl, gpar, mpar, rwsl, dfe1, dfe2, fe1data[:,2], fe1data[:,5], fe2data[:,5], epint, rwint
