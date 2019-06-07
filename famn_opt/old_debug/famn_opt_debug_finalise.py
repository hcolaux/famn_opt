# -*- coding: utf-8 -*-
"""
Console mode execution of the FAM-N optimisation program.
Fil-in the prameters and press the run button/F5
@author: Henri
"""

import subprocess
from famn_core import *
import matplotlib.pyplot as plt
import numpy
import copy

# self = copy.deepcopy(famn.output['FAMNsteps'][0])
self = copy.deepcopy(famn)

# ----------------------------------------------------------
# Final detect coherence
self.final_coherence = self.output['list_detect_density'][  self.output['no_FAMN_steps']-1  ]

# ----------------------------------------------------------
# Number of pulses
self.total_nopules = numpy.zeros(self.output['no_FAMN_steps'], dtype = 'int32')
for FAMNstepcount in range(self.output['no_FAMN_steps']): 
    self.total_nopules[FAMNstepcount] = int(self.output['FAMNsteps'][FAMNstepcount].data['tnopules'])
    if self.output['is_succesion'][FAMNstepcount]:
        self.total_nopules[FAMNstepcount-1] -= 1

# ----------------------------------------------------------
# Select the final coherence evolution
self.select_coherences = []
self.pulse_durations = []

for FAMNstepcount in range(self.output['no_FAMN_steps']):
    # Swipe all separated FAM-N of the optimisation
    if self.output['do_final_sim'][FAMNstepcount]:
        # 
        for FAMNpulses in range(self.output['FAMNsteps'][FAMNstepcount].data['tnopules']):
            self.select_coherences.append(get_all_detect_element(
                self.output['FAMNsteps'][FAMNstepcount].data['final_full_density'][FAMNpulses],
                self.final_coherence))
            self.pulse_durations.append(len(self.select_coherences[-1]))
    else:
            
        for FAMNpulses in range(self.output['FAMNsteps'][FAMNstepcount].data['tnopules']):
            self.select_coherences.append(get_all_detect_element(
                self.output['FAMNsteps'][FAMNstepcount].data['selected_full_density'][FAMNpulses],
                self.final_coherence))
            self.pulse_durations.append(len(self.select_coherences[-1]))


# -----------------------------------------------------
# Final signal improvements    
self.final_signal_improvement = numpy.real(self.select_coherences[-1][ self.output['FAMNsteps'][-1].data['max_pos'][-1] ])/\
                                numpy.real(self.select_coherences[0][ self.output['FAMNsteps'][0].data['max_pos'][0] ])

# -----------------------------------------------------
# Total number of points
self.total_nopoints = numpy.full((self.output['no_FAMN_steps']),int(0))  
for FAMNstepcount in range(self.output['no_FAMN_steps']):
    self.total_nopoints[FAMNstepcount] = numpy.sum(self.output['FAMNsteps'][FAMNstepcount].data['pulse_durations_nop'])

# Total duration: * self.output['FAMNsteps'][FAMNstepcount].opt['time_increment']

# -----------------------------------------------------
# Text fractions of the final pulse for output in Topspin and Delta
self.pulse_frac = []
for FAMNstepcount in range(self.output['no_FAMN_steps']): 
    # Pre-allocate
    temp = numpy.full( self.total_nopules[FAMNstepcount] , None )
    
    # Store in a temporary array
    for FAMNpulses in range(self.total_nopules[FAMNstepcount]):
        temp[FAMNpulses] = (str(self.output['FAMNsteps'][FAMNstepcount].data['pulse_durations_nop'][FAMNpulses]) +
                '/' + str(self.total_nopoints[FAMNstepcount]))
    # -----------------
    self.pulse_frac.append(temp)
    
# -----------------------------------------------------
# First absolute point of each pulse
self.initial_time = numpy.zeros(self.output['no_FAMN_steps'])
for FAMNstepcount in range(1,self.output['no_FAMN_steps']): 
    self.initial_time[FAMNstepcount] =\
    1.0e6 * self.output['FAMNsteps'][FAMNstepcount].opt['time_increment'] *\
    self.initial_time[FAMNstepcount-1] +\
    self.output['FAMNsteps'][FAMNstepcount].data['pulse_durations_nop'][-1]
    
    
# -----------------------------------------------------
# Pre-allocate figure
self.final_figure, self.final_axes = plt.subplots(1,1, figsize = (10,5))
self.final_axes.set_title('Final FAM-N pulse')
self.final_axes.set_xlabel('Time /μs')

counter = int(0)
for FAMNstepcount in range(self.output['no_FAMN_steps']):
    # Plot the first pulse
    self.final_axes.plot(self.initial_time +\
            1.0e6 *self.output['FAMNsteps'][FAMNstepcount].opt['time_increment'] *\
            numpy.arange(self.pulse_durations[counter]),
        self.select_coherences[counter])
    # Plot subsequent pulses
    counter =+ 1
    for FAMNpulses in range(self.output['FAMNsteps'][FAMNstepcount].data['tnopules']-1):
        # Plot curve
        self.final_axes.plot(self.initial_time + 1.0e6 * \
                self.output['FAMNsteps'][FAMNstepcount].opt['time_increment'] *numpy.arange(\
               self.output['FAMNsteps'][FAMNstepcount].data['inv_pos_abs'][FAMNpulses],
               self.output['FAMNsteps'][FAMNstepcount].data['inv_pos_abs'][FAMNpulses] +\
               self.pulse_durations[counter]),
               self.select_coherences[counter])
        # Trace point of phase inversion
        self.final_axes.axvline(x = self.initial_time[FAMNstepcount] +\
                        1.0e6 *self.output['FAMNsteps'][FAMNstepcount].opt['time_increment'] *\
                        self.output['FAMNsteps'][FAMNstepcount].data['inv_pos_abs'][FAMNpulses]        
                        , color = 'black')
        # -----
        counter =+ 1
# ---------------------------------------------------------------------
# Plot position of first maximum
self.final_axes.axvline(1.0e6 *self.output['FAMNsteps'][FAMNstepcount].opt['time_increment'] *\
                        self.output['FAMNsteps'][0].data['max_pos_abs'][0],
        color = 'red', linestyle = '--')
# Plot position of last maximum
self.final_axes.axvline(self.initial_time[-1] +\
                        1.0e6 *self.output['FAMNsteps'][FAMNstepcount].opt['time_increment'] *\
                        self.output['FAMNsteps'][-1].data['max_pos_abs'][-1],
        color = 'red', linestyle = '--')






