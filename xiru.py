from setup_converge import initial_options, std_anal_init
from main_routines_conv import converge_main, final_remarks
from differential import diff_main
from converge_errors import make_errors

# interprets command-line arguments
par0, makeplot, argsalpha, differential = initial_options()

if differential == True:
    argsalpha = False # forces solar-scaled atmosphere for differential analysis
    # runs differential analysis:
    pn, qdr, final_jac = diff_main(par0, makeplot)
else:
    # creates x0 and x1 vectors; creates Jacobian
    p1, Jn, ovec, pn_init, quadr = std_anal_init(par0, argsalpha)
    # main algorithm for standard analysis
    pn, qdr, final_jac = converge_main(p1, Jn, ovec, pn_init, quadr, makeplot, argsalpha)

# calculates internal uncertainties
datm = make_errors(final_jac, diff_analysis=differential)

# prints information; saves the results; plots iterations if requested
final_remarks(pn, qdr, final_jac, makeplot, argsalpha, differential, datm)
