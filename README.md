# Xiru - a Python wrapper for MOOG

**Xiru** (pronounced *"she-ROO"*)[\*] is a wrapper for [MOOG](https://www.as.utexas.edu/~chris/moog.html) driver *abfind*, to be used in MOOG silent mode, written in Python3. Its main purpose is to find the spectroscopic parameters of a star (Teff, log g, [M/H], and microturbulent velocity) from excitation/ionisation balance of Fe I/II spectral lines. It uses an implementation of [Broyden's method](https://www.ams.org/journals/mcom/1965-19-092/S0025-5718-1965-0198670-6/home.html) to achieve convergence as quick as possible. Both classic and differential analysis are supported.

It currently runs natively with model atmospheres in KURUCZ and KURTYPE formats (see MOOG documentation for detailed info). Native support for MARCS models is planned for future versions, nevertheless the user can add a custom plugin for them in the current version. It is assumed that the user already has some previous experience with MOOG.

[\*] *Xiru* is a word in Gaúcho dialect of Brazilian Portuguese, meaning "friend" or "mate" (from Guarani *che irũ*). It is also used to denote a wise and experienced person.

## Installation and configuration

(1) Clone the files to any folder in your system.

If you intend to use [Castelli/Kurucz models](https://wwwuser.oats.inaf.it/castelli/grids.html) go to step (2), otherwise jump to step (4):

---

(2) Run the routine kurucz_extractor.py

$ python kurucz_extractor.py

It downloads the models with vturb=2.0 and alpha=0.0 by default and saves them to disc in a "models" subfolder in a format friendly to Xiru's interpolator. The default values for vturb and alpha can be changed according to user preferences in lines 117-119. 

(3) In interpolator.py, edit the path to the model grid (line 9). Jump to step (6).

---

(4) If you intend to use another grid of model atmospheres:

Customise the custom_conv.py module as instructed in its commented lines. Barely any knowledge of Python is required if the user has some external interpolator routine (see example in the module).

(5) In moognmodels.py module, *run_model* method, change the *modeltype* argument from 'kurucz' to 'custom'.

---

(6) In moognmodels.py module, *run_moog* method, edit the *moogpath* argument to the path of your silent mode version of MOOG.

### A few notes on filenames:

- Xiru assumes that the input file read by silent MOOG is called **batch.par**. If you use a different filename in your MOOG installation, please edit lines 19, 20, 26 and 27 in differential.py accordingly, replacing "batch.par" with the alternative filename.

- The native interpolator assumes that the *model_in* file in MOOG input is called **MODEL**. (interpolator.py, line 370)

- The assumed *summary_out* filename is **a.out**. If you intend to use a different filename for summary_out, use "$ grep -n a.out \*.py" in Xiru root folder to locate the lines that must be edited.

## Usage

First, prepare your *lines_in* file and batch.par as in a conventional usage of MOOG. It is suggested to create a separate folder for each star that the user is going to study. Xiru must be executed from this star-specific folder.

### Standard (non-differential) spectroscopic method:

You need to input a first guess for the atmospheric parameters.

$ python xiru.py [Teff] [logg] [metal] [micro]

Units are K, dex (cm/s), dex and km/s, respectively.

After few iterations the code will print the results and the respective internal uncertainties.

Equivalent width data from the spectrum of Arcturus ([Hinkle et al, 2000](http://ast.noao.edu/data/other)), measured by me, is included with the code as an example for the standard analysis. Simply run Xiru from its folder to see how to code works.

### Differential analysis:

Two files must be prepared for the reference star: (i) ref.par is the MOOG input file for reference star (e.g., the Sun), with its respective *lines_in* file containing equivalent widths measured by the user in the reference star spectrum, while (ii) refatm is an ASCII file contaning the atmospheric parameters of the reference star in a single line, e.g.:

5777 4.438 0.00 1.0

These files must be placed in the same folder as batch.par (where the user is executing Xiru) for each star in the analysis (i.e., if the user has a set of solar twins, each one must possess a folder with its respective batch.par, ref.par, and refatm). A third file related to the reference star is ref.out, which corresponds to the *summary_out* listed in ref.par for the reference star, and it does not need to be created by the user.

If the user decides to adopt different filenames for the reference star, line 11 from differential.py must be edited accordingly.

As positional argument you need to input a first guess for the atmospheric parameters. This time the command-line option -d is required to activate the differential analysis.

$ python xiru.py [Teff] [logg] [metal] [micro] -d

A differential Boltzmann diagram will be saved in a file called "feplot.png".

### Command-line options:

-d: activates differential analysis

-a: uses Salaris 1993ApJ...414..580S formula to correct metallicity for alpha-ehnancement (disabled for differential analysis).

-p: plots iteration history of EP Slope, Delta_Fe, Delta_M_H, and RW slope

## Citation:

Puls et al. in prep.

## Acknowledgement

Thanks to Dr. Lorenzo Spina for kindly providing tables with equivalent widths data from solar twins for differential analysis testing.

## Contact

arthur.alencastropuls at anu.edu.au
