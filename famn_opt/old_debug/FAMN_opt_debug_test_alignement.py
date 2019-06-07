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
self = copy.deepcopy(famn.output['FAMNsteps'][0])
self.FAMN_main.parameters['no_SIMPSON_execution'] = False

#self = fnf.FAMN_main()
#self.loadstate('test_conc')


# Initialise some arrays
self.saved_density = []
self.saved_detect = []

self.temp['cur_name'] = 'debug_test1'



# -----------------------------------------------------------------

# Write SIMPSON input file
self.write_simpson(self.FAMN_main.output['list_initial_density'][self.opt['order_index']], # Initial density matrix
                   self.opt['nop_CW'], # Number of points
                   self.FAMN_main.output['dim_density'], # Dimesion of the density matrix
                   ) # self.temp['pulse_phase'] # Pulse phase

# SIMPSON Execution
if not self.FAMN_main.parameters['no_SIMPSON_execution']:
    subprocess.call('simpson ' + './' + self.FAMN_main.parameters['simpson_tempfolder'] +\
          '/' + self.temp['cur_name'] + '.in')

# Read in the file
self.simpson_out_full1 = fnf.readin_simpson_out(('./' + self.FAMN_main.parameters['simpson_tempfolder'] +\
      '/' + self.temp['cur_name'] + '.out'),
    self.opt['nop_CW'], self.FAMN_main.output['dim_density'],
    self.FAMN_main.output['size_density'],
    )


self.simpson_out_full1 = numpy.concatenate((
        [numpy.conj(self.FAMN_main.output['list_initial_density'][self.opt['order_index']])],
         self.simpson_out_full1))

# Get only the detect element from this numpy array
self.detect1 = numpy.real(self.get_detect_element(\
         self.simpson_out_full1,self.opt['nop_CW']+1))



plt.plot(self.detect1)



# -----------------------------------------------------------------

self.temp['cur_name'] = 'debug_test2'

nop = 5

# Write SIMPSON input file
self.write_simpson(numpy.conj(self.simpson_out_full1[nop]), # Initial density matrix
                   self.opt['nop_CW'], # Number of points
                   self.FAMN_main.output['dim_density'], # Dimesion of the density matrix
                   ) # self.temp['pulse_phase'] # Pulse phase

# SIMPSON Execution
if not self.FAMN_main.parameters['no_SIMPSON_execution']:
    subprocess.call('simpson ' + './' + self.FAMN_main.parameters['simpson_tempfolder'] +\
          '/' + self.temp['cur_name'] + '.in')




# Read in the file
self.simpson_out_full2 = fnf.readin_simpson_out(('./' + self.FAMN_main.parameters['simpson_tempfolder'] +\
      '/' + self.temp['cur_name'] + '.out'),
    self.opt['nop_CW'], self.FAMN_main.output['dim_density'],
    self.FAMN_main.output['size_density'],
    )


self.simpson_out_full2 = numpy.concatenate((
        [numpy.conj(self.simpson_out_full1[nop])],
         self.simpson_out_full2 ))

self.detect2 = self.get_detect_element(self.simpson_out_full2,self.opt['nop_CW']+1)



plt.plot(numpy.arange(nop,nop+self.opt['nop_CW']+1),self.detect2)
    
    
    