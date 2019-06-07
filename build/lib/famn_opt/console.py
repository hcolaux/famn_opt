# -*- coding: utf-8 -*-
"""
Console mode execution of the FAM-N optimisation program.
Fil-in the prameters execute this script (console)
or press the run button/F5 (Spyder/pyCharm/VSstudio …)
@author: Henri Colaux
"""

# ----------------------------------------------------------------------------
# Create a new class instance
from opt import FAMN_main
famn = FAMN_main()
famn.param['simpson_basename'] = 'famn'

# ----------------------------------------------------------------------------
# (Optionnal) Path where all output files are exported.
famn.param['export_filepath'] = 'R:/output/'

# ----------------------------------------------------------------------------
# Input parameters
famn.param['nucleus'] = '23Na' # Nucleus, as written in Simpson
famn.param['proton_larmor_frequency'] = 600e6 # Proton Larmor frequency (Hz) 
famn.param['MAS'] = 12500 # Magic angle spinning rate (Hz) 
famn.param['RF'] = 100000 # Radio-frequency field ν1 (Hz)  
famn.param['CQ'] = 1e6 # Quadrupolar coupling constant CQ (Hz) 
famn.param['ηQ'] = 0 # Asymmetry parameter ηQ 

# ----------------------------------------------------------------------------
# Transfer parameters

# To fill for a standard FAM-N optimisation
famn.param['transition'] = 3 # Transition of the FAM-N method.

# To fill in for complex FAM-N optimisation
famn.param['complex_optimisation'] = '' # '5>3;3>1'
famn.param['single_pulse_nutation_list'] = ''  # '120,75' # If complex optimisation, string list of the single pulse nutation angles for each step

# ---------------------------------------------
# Crystal file
famn.param['no_cryst'] = 66 # Number of crystallites
famn.param['no_γ'] = 4 # Number of γ-angles for the SIMPSON optimisation

# ---------------------------------------------
# Advanced user parameters.
famn.param['min_improvement'] = 1.1 # Minimum improvement between two subsequent pulses
                                    # for the simulation to be allowed to add additionnal pulses
famn.param['max_nopulses'] = 30 # Maximum number of pulses that the program is allowed to output.
famn.param['nop_360'] = 180 # Number of points par 360° nutation.
famn.param['min_time_increment'] = 1e-8 # Minimum time increment (seconds)
famn.param['maxdt_ratio'] = 20 # Ratio between MAS period and maxdt in SIMPSON.
famn.param['num_cores'] = 0 # Number of cores on which the optimisation is programmed

# ---------------------------------------------
# Debug
famn.param['no_SIMPSON_execution'] = True
famn.param['debug_figure'] = True
famn.param['do_every_cryst'] = True
famn.param['keep_simpson_files'] = True

# ----------------------------------------------------------------------------
# Execution of the code
famn.initialise()
famn.FAMN_loop()
famn.finalise()
famn.screen_report()


## Export Bruker
#famn.export_topspin_MQMAS('bruker_MQMAS')
#famn.export_topspin_pulselist('bruker_pulselist')
#famn.export_topspin_MQF('bruker_MQF')
#famn.export_topspin_MQCPMGMAS('bruker_MQCPMGMAS')
#
## Export Jeol
#famn.export_delta_MQMAS('jeol_MQMAS')
#famn.export_delta_MQF('jeol_MQF')
#famn.export_delta_pulselist('jeol_pulselist')

# ----------------------------------------------------------------------------