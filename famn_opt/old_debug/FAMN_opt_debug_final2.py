# -*- coding: utf-8 -*-
"""
Console mode execution of the FAM-N optimisation program.
Fil-in the prameters and press the run button/F5
@author: Henri
"""

import numpy
import subprocess
import FAMN_opt as fnf
import matplotlib.pyplot as plt

import copy
self = copy.deepcopy(famn.output['FAMNsteps'][1])
self.FAMN_main.parameters['no_SIMPSON_execution'] = False

#self = fnf.FAMN_main()
#self.loadstate('test_conc')


# Initialise some arrays
self.saved_density = []
self.saved_detect = []

self.temp['cur_name'] = self.finalise_filename(1)

print('Final simulation for pulse n°1')
# Do the simulations for the single pulse
            


# -----------------------------------------------------------------

# Write SIMPSON input file
self.write_simpson(self.FAMN_main.output['FAMNsteps']\
                   [self.opt['order_index']-1].data['final_density'], # Initial density matrix
                   self.opt['nop_CW'], # Number of points
                   self.FAMN_main.output['dim_density'], # Dimesion of the density matrix
                   ) # self.temp['pulse_phase'] # Pulse phase

# SIMPSON Execution
if not self.FAMN_main.parameters['no_SIMPSON_execution']:
    subprocess.call('simpson ' + './' + self.FAMN_main.parameters['simpson_tempfolder'] +\
          '/' + self.temp['cur_name'] + '.in')

# Read in the file
self.temp['simpson_out_full'] = fnf.readin_simpson_out(('./' + self.FAMN_main.parameters['simpson_tempfolder'] +\
      '/' + self.temp['cur_name'] + '.out'),
    self.opt['nop_CW'], self.FAMN_main.output['dim_density'],
    self.FAMN_main.output['size_density'],
    )

# Save the current density evolution
self.saved_density.append(self.temp['simpson_out_full'])

# Test figure
plt.figure(figsize = (12,6))
self.saved_detect.append(self.get_detect_element(self.saved_density[0],self.opt['nop_CW']))
plt.plot(1e6 * self.opt['time_increment'] * numpy.arange(0,self.opt['nop_CW']),self.saved_detect[0])

# -----------------------------------------------------------------
# Loop over the number of pulses
for truc in range(1,self.data['tnopules']):
    
    print('Final simulation for pulse n°' + str(truc+1))
    # self.temp['pulse_phase'] = not self.temp['pulse_phase'] # Inverse phase # Not used
    
    # Change filename
    self.temp['cur_name'] = self.finalise_filename(truc+1)
    
    # Write SIMPSON input file
    self.write_simpson(self.saved_density[truc-1][self.data['max_pos'][truc-1]+self.data['inv_pos'][truc-1]],
                       # Initial density matrix
                       self.data['pulse_dur'][truc], # Number of points
                       self.FAMN_main.output['dim_density'], # Dimesion of the density matrix
                       ) # self.temp['pulse_phase'] # Pulse phase
    
    # SIMPSON Execution
    if not self.FAMN_main.parameters['no_SIMPSON_execution']:
        subprocess.call('simpson ' + './' + self.FAMN_main.parameters['simpson_tempfolder'] +\
              '/' + self.temp['cur_name'] + '.in')
    
    # Read in the file
    self.temp['simpson_out_full'] = fnf.readin_simpson_out(('./' + self.FAMN_main.parameters['simpson_tempfolder'] +\
          '/' + self.temp['cur_name'] + '.out'),
        self.data['pulse_dur'][truc], self.FAMN_main.output['dim_density'],
        self.FAMN_main.output['size_density'],
        )
    
    # Save the current density evolution
    self.saved_density.append(self.temp['simpson_out_full'])
    self.saved_detect.append(self.get_detect_element(self.saved_density[truc],self.data['pulse_dur'][truc]))



    plt.plot(1e6 * self.opt['time_increment'] * numpy.arange(self.data['inv_pos_abs'][truc-1],
                          self.data['inv_pos_abs'][truc-1] + self.data['pulse_dur'][truc]),
                            self.saved_detect[truc])
    
    
    
    
    
    
    
    