# this module contains a custom method to be written by the user
# the main purpose here is to implement atmospheric model interpolators other than the native one (interpolator.py) provided with this code
# what the user needs to do:
# edit the run_custom_interpolator method as you need
# the atmparam variable MUST be a numpy array
# if you already have a grid interpolator (written in any language),
# you can, for instance, instruct the run_custom_interpolator to tell the system to run it
# in any case, a MOOG 'model_in' file must be created in the system

def run_custom_interpolator(atmparam):
    # example:
    # if the user has a code called 'interpolator.a' that
    # interpolates the atmospheric model from some grid on some 'model_in' file defined in the MOOG input file
    # and uses the atmospheric parameters teff logg metal micro as command-line arguments
    # write:
    # import os
    # os.system('interpolator.a %.1f %.3f %.3f %.4f' % (atmparam[0], atmparam[1], atmparam[2], atmparam[3]))
    return void