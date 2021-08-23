import os
#from kurtype import run_kurtype_interpolator
from interpolator import run_interpolator
from custom_conv import run_custom_interpolator

def run_moog(moogpath='/path/to/MOOG'):
    "Basically...runs MOOG"
    "The 'moogpath' string must be the path to MOOG in your system"
    os.system('%s' % moogpath)

def run_model(atmin, modeltype='kurucz'):
    "Runs the interpolator that creates the atmospheric model for MOOG using some grid"
    
    # if you are going to use castelli/kurucz models, a built-in interpolator is available for you :)
    # if your favourite model grid is different (MARCS, etc.), edit the modeltype string in the method argument
    # as well as the run_custom_interpolator method in custom_conv.py
    if modeltype == 'kurucz':
        run_interpolator(atmin)
    elif modeltype == 'custom':
        run_custom_interpolator(atmin)
    else:
        print('Invalid option for model interpolator!')
        raise SystemExit()