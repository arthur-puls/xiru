# Downloads tables of temperature-pressure relations from castelli website
# Default values for alpha, vturb and metallicity can be changed in lines 117-119
from urllib.request import urlopen
from pathlib import Path
from os import mkdir
import numpy as np

def get_params(metal, alpha=0.0, vturb=2.0):
    "metal: string following the syntax from kurucz/castelli websites, e.g., p05 for +0.5, m15 for -1.5"
    "alpha: alpha-enhancement, 0.0 or 0.4 dex"
    "vturb: 0.0 or 2.0"
    if alpha == 0.0:
        astr = 'k'
    elif alpha == 0.4:
        astr = 'ak'
    else:
        print('Extractor: Invalid alpha!')
        raise SystemExit()
    if vturb == 0.0:
        vstr = '0'
    elif vturb == 2.0:
        vstr = '2'
    else:
        print('Extractor (get_params): Invalid vturb!')
    mpath = 'https://wwwuser.oats.inaf.it/castelli/grids/'
    grid_set = 'grid%s%s%sodfnew' % (metal, astr, vstr)
    model_set  = 'a%s%s%sodfnew.dat' % (metal, astr, vstr)
    return '%s%s/%s' % (mpath, grid_set, model_set)

def cut_extras(model, fout):
    "Cuts header/footer of each model"
    "model: string variable"
    lines = model.split('\n')
    # counts number of teff/logg models
    N = len(lines)
    ilines = []; tlines = []; ntau = []; teff = []; logg = []
    for i in range(N):
        try:
            if lines[i].split()[0] == 'READ':
                ilines.append(i+1)
                ntau.append(lines[i].split()[2])
            elif lines[i].split()[0] == 'PRADK':
                tlines.append(i-1)
            elif (lines[i].split()[0] == 'TEFF') and (lines[i].split()[2] == 'GRAVITY'):
                teff.append(lines[i].split()[1])
                logg.append(lines[i].split()[3])
            else:
                pass
        except:
            pass
    # checks if number of headers/tails is the same
    if len(ilines) == len(tlines):
        M = len(ilines)
    else:
        print('Extractor (cut_extras): header/footer count error!')
        raise SystemExit()
    # proceeds to select only the lines/columns relevant to MOOG and print them to an output file
    with open(fout, 'w') as f:
        for i in range(N):
            for j in range(M):
                if (i == ilines[j]):
                    print('MODEL %s' % ntau[j], file=f)
                    print('%8s %8s' % (teff[j], logg[j]), file=f)
                if (i >= ilines[j]) and (i <= tlines[j]):
                    line = lines[i].split()[:7]
                    for k in range(len(line)):
                        print('%s ' % line[k], end='', file=f)
                    print('', file=f)
    return

def get_vmp():
    teff = np.arange(3500, 6250, 250)
    logg = np.array(['00','05','10','15','20','25','30','35','40','45','50'])
    mpath = 'https://wwwuser.oats.inaf.it/castelli/grids/gridm40k2odfnew/'
    model = ''
    for i in range(teff.shape[0]):
        for j in range(logg.shape[0]):
            u = '%sam40t%4ig%2sk2odfnew.dat' % (mpath, teff[i], logg[j])
            with urlopen(u) as f:
                output = f.read().decode('utf-8')
            model = model+output
    return model

def check_dir():
    "Checks if output directory exists"
    "Creates one if false"
    dir_exists = Path('models').is_dir()
    if dir_exists:
        pass
    else:
        print('Creating output directory...', end='')
        mkdir('models')
        print('Done!')
    return

def extractor(metal, vturb, alpha):
    "Loads files in string, selects relevant information, converts to interpolator-friendly format, and saves to disc"
    for i in range(len(metal)):
        u = get_params(metal[i], alpha=alpha, vturb=vturb)

        try:
            if (metal[i] == 'm40') and (alpha == 0.0) and (vturb == 2.0):
                output = get_vmp()
                fout = './models/am40'
                cut_extras(output, fout)
            else:
                with urlopen(u) as f:
                    output = f.read().decode('utf-8')
                fout = './models/a%s' % metal[i]
                cut_extras(output, fout)
        except:
            pass # grid not found

if __name__ == "__main__":
    # list with the required metallicities
    # follows the castelli/kurucz naming convention
    metal = ['m40', 'm25', 'm20', 'm15', 'm10', 'm05', 'p00', 'p02', 'p05']
    vturb = 2.0
    alpha = 0.0
    
    check_dir()
    
    extractor(metal, vturb, alpha)
