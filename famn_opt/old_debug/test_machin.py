# -*- coding: utf-8 -*-
"""
Console mode execution of the FAM-N optimisation program.
Fill in the prameters execute this script (console)
or press the run button/F5 (Spyder/pyCharm/VSstudio…)
"""

# Initialise the script
from FAMN_opt import FAMN_main
test_machin = FAMN_main()
test_machin.parameters['simpson_basename'] = test_machin
# ----------------------------------------------------------------------------
# Input parameters

 # Spin system input parameters
test_machin.parameters['nucleus'] = '23Na' # Nucleus, as written in Simpson
test_machin.parameters['proton_larmor_frequency'] = 600e6 # Proton Larmor frequency (Hz) 
test_machin.parameters['MAS'] = 12500 # Magic angle spinning rate (Hz) 
test_machin.parameters['RF'] = 100000 # Radio-frequency field ν1 (Hz)  
test_machin.parameters['CQ'] = 1.0e6 # Quadrupolar coupling constant CQ (Hz) 
test_machin.parameters['ηQ'] = 0 # Asymmetry parameter EtaQ 

# ----------------------------------------------------------------------------
# Transfer parameters

# To fill for a standard FAM-N optimisation
test_machin.parameters['transition'] = 3 # Transition of the FAM-N method.

# To fill in for complex FAM-N optimisation
test_machin.parameters['complex_optimisation'] = ''
test_machin.parameters['single_pulse_nutation_list'] = '', # If complex optimisation, string list of the single pulse nutation angles for each step

# ---------------------------------------------
# Advanced user parameters.

test_machin.parameters['min_improvement'] = 1.1 # Minimum improvement between two subsequent pulses
                                    # for the simulation to be allowed to add additionnal pulses
test_machin.parameters['max_nopulses'] = int(30) # Maximum number of pulses that the program is allowed to output.
test_machin.parameters['nop_360'] = 180 # Number of points par 360° nutation.
test_machin.parameters['min_time_increment'] = 1e-8 # Minimum time increment (seconds)
test_machin.parameters['maxdt_ratio'] = int(20) # Ratio between MAS period and maxdt in SIMPSON.
test_machin.parameters['crystfile'] = 'rep66' # Crystal file used for the SIMPSON optimisation
test_machin.parameters['no_γ'] = 4 # Number of γ-angles for the SIMPSON optimisation
test_machin.parameters['num_cores'] = int(0) # Number of cores on which the optimisation is programmed

# ----------------------------------------------------------------------------
# Execution of the code
test_machin.initialise()
test_machin.FAMN_loop()
test_machin.finalise()

# ----------------------------------------------------------------------------
