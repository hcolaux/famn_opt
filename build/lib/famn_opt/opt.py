# -*- coding: utf-8 -*-
'''
Program containing classes and methods for FAM-N optimisation. 
See "QuickStart.txt" for instructions.
# -----------------------------------------------------------
Henri Colaux
University of St Andrews
Fife, Scotland, United Kingdom
e-mail: henri.colaux@normalesup.org
'''

# Import packages
import sys, os
import subprocess
import numpy 
import lxml.etree as xmlp # To read xml files
import matplotlib.pyplot as plt # To plot curves
import pickle # To save classes

# GUI package
from PyQt5 import QtCore, QtGui, QtWidgets

# Function file
from famn_opt import functions as fmc

# -----------------------------------------------------------------------------
# Main class 

class FAMN_main:
    ''' Definition of the class that contains all informations, elements and methods
    for FAM-N optimisation.'''
    
    def __init__(self): 
        
        self.param = { # ---------------------------------- 
                # Instance parameters
                'FAMopt_version' : '0.3a', # Version of this program
                'simpson_tempfolder' : 'famn_temp', # Name of the SIMPSON folder
                'simpson_basename' : '', # Name of this instance
                
                # ---------------------------------------------
                # Main user parameters.
                'nucleus' : '', # Nucleus of the optimisation
                'proton_larmor_frequency' : float(0) , # Proton Larmor frequency (Hz) 
                'MAS' : float(0), # Magic angle spinning rate (Hz) 
                'RF': float(0), # Radio-frequency field ν1 (Hz) 
                'CQ' : float(0), # Quadrupolar coupling constant CQ (Hz) 
                'ηQ' : float(0), # Asymmetry parameter EtaQ 
                
                # --------------
                # Optimisation option
                'spin' : '', # Nucleus spin, filled automatically
                'transition' : int(0), # Nucleus spin, filled automatically
                'complex_optimisation' : '', # Advanced parameter, to make looped optimisation.
                'single_pulse_nutation_list' : '', # If complex optimisation, string list of the single pulse nutation angles for each step
                
                # ---------------------------------------------
                # Advanced user parameters.
                'nop_360' : int(180), # Number of points par 360° nutation.
                'maxdt_ratio' : int(20), # Ratio between MAS period and maxdt in SIMPSON.
                'max_nopulses' : int(26), # Maximum number of pulses that the program is allowed to output.
                'min_time_increment' : float(1e-8), # Minimum time increment for FAM-N optimisation
                'num_cores' : int(0), # Number of cores on which the optimisation is programmed
                'set_offset_to_qis' : True, # Set the offset to the quadrupolar induced shift (QIS)
                'export_filepath' : '', # Path where all output files are exported
                # Tolerance counters. Set the number of time increments
                #the programme carries on calculation despite a decreasing amount of single quantum.
                'tol_count+' : 5, # Tolerance for increasing pulse length
                'tol_count−' : 3, # Tolerance for decreasing pulse length
                
                # Crystal file data
                'do_every_cryst' : True, # 1 <=> Do the crystalite-per-crystalite optimisation
                'cryst_name' : 'famn_cryst', # Name of the crystal file if applicable
                'crystfile' : '', # Name of the crystal file written in SIMPSON simulations
                'only_beta' : False, # 1 <=> Use only beta in the simulation
                'cryst_gen_method' : 'golden', # Method used to genenerate the crystal file. Can be 'golden' or 'zcw'(NB: add repulsion)
                'no_gamma' : int(20), # Number of gamma-angles for the SIMPSON optimisation
                'no_cryst' : int(34), # Number of crystallites used to generate the file (input)
                
                # Debug
                'console_warning' : True, # Display warnings from the console
                'keep_simpson_files' : True, # Keep the intermediate simpson file
                'no_SIMPSON_execution' : False, # Disable simpson execution
                'debug_figure' : False, #Show debug figures
                'klupolnghsietd' : 'vdrffuadbpt', # dltrsaiqeqldostsrei
                # Others
                'tab' : '     ', # Tabulation before numbers in the output
                'opt_method' : 'cheby2', # Calculation method used internally in SIMPSON
                }

        # Main data arrays 
        self.output = {
                'Larmor_frequency' : float(0), # Larmor frequency of the nucleus used if the optimisations
                'no_FAMN_steps' : int(0), # Number of FAM-N steps
                'transfer_type' : int(0), # Returns the transfer type of the current optimisation.
                                          # It can be 'conversion', 'excitation' or 'unknown'
                'simple_opt' : True, # True <=> a single optimisation was carried out.
                'crystal_list' : [] , # List of crystals file
                'do_final_sim' : [], # If complex optimisation: list of pulses for which to redo the simulation
                                    # with the final density of the previous step as the initial density
                'is_succesion' : [], # If complex optimisation, turns to True if
                                    # the current step is the succession of the previous
                'completed' : False, # Turns to "True" if the optimisation completed successfully
                
                # Other parameters
                'offset' : float(0), # Offset during the optimisations
                'no_cryst' : int(0), # Number of crystallite used in the simulation
                'dim_density' : int(0), # Dimension of density matrix
                'size_density' : int(0), # Size of the density matrix
                'list_initial_density' : [], # List of the initial density matrices
                'list_detect_density' : None, # Detect element of the current FAM-N step
                'list_single_pulse_nutation' : [[int(0)]], # List of the total nutation used for the first pulse.
                'total_execution' # Timer
                # NB: The convention for the detection of this program is inverse to that of SIMPSON
                
                # Warning bools
                'using_min_timeinc' : False, # Turn to True minimum time increment

                # Error booleans
                'SIMPSON_installed' : True, # 0 <=> SIMPSON is not installed
                'non_existing_nucleus' : False, # 1 <=> Non-existing nucleus
                'impossible_transition' : False, # =1 <=> The transition is not possible for a given spin
                'crystallites_not_integer' : False, # =1 <=> Turns to True if the number of crystallite is not integer 
                'not_enough_crystallites' : False, # =1 <=> Not enough crystallite were given
                'partial_infos_complex' : False, # Partial informations were provided for complex FAM-N
                'complex_parse_error' : False, # Indicate that there was an error in parsing the input for complex FAM-N
                'nonexisting_transition' : False, # =1 <=> The transition used as input is not valid
                'non_quad_spin' : False, # =1 <=> The nucleus is not half-integer quadruporal spin 
                'Non_consistent_number_of_famn_steps' : False, # The number of FAM-N steps if different
                # between 'complex_optimisation' and 'single_pulse_nutation_list'
                'loading_error' : False, # 1 <=> Attempt to load a non-valid .pkl file
                'existing_allexport_filepath' : True, # Turns to False is the filepath is not correct
                
                # Array containing the very final parameters
                'final_coherence' : None, # Matrix with the final,final detected elements
                'total_nopules' : int(0), # Total number of pulses in all steps
                'select_coherences' : [], # List of all select coherence (unsorted)
                'pulse_durations' : [], # List of all pulse durations in number of points (unsorted)
                'final_TIR' : float(0), # Signal improvement for all steps
                'total_nop' : None, # Total number of points of all steps (sortes)
                'pulse_frac' : None, # Pulse fraction of for all steps for output in Topspin and Delta
                'max_amplitudes' : [], # Max amplitudes calculated from the final detect element
                'initial_time' : [], # Stores the first initial time for each step
                }
        
        self.figures = { # Contains figure properties 
                'allsteps_fig' : None, # Figure containing all steps
                'allsteps_ax' : None, # corresponding axis
                }
        
        self.temp = { 'noprev_pulses' : int(0), # Number of pulses before the previous one
                     'time_increment_min' : float(0), # Number of points determined from the minimum time increment
                     'time_increment_calc' : float(0), # Number of points determined from the nutation
                     }

        # Number one most important:
        self.FAMNsteps = [] # List stocking every step of FAM-N optimisation
        
    # -------------------------------------------------------------------------
    # Define FAMN_step class
    class FAMN_step():
        ''' Contains the results for one optimisation for one transfer'''
        
        def __init__(self,FAMN_main):
            
            # Get the data from the main class into the step class
            self.FAMN_main = FAMN_main
            
            self.opt = { # 
                    # Commun elements
                    'order_index' : int(0), # The order of this step in the totpal optimisation
                    'detect_density' : None, # Detected density matrix of the FAM-N step
                    'maxdt' : float(0) , # Same as the maxdtparameter as defined by SIMPSON
                    'time_increment' : float(0), # Time increment used in the optimisation
                    'nop_CW' : int(0), # Number of points of the single pulse acquisition
                    'max_nopules' : int(0), # Maximum number of pulses for this instance
                    # Number of points of the optimisation of the single pulse: self.output['list_single_pulse_nutation'][stepcnt]
                    # Number of points of the SIMPSON simulation
                                        
                    # SIMPSON file
                    'head_1' : '', # Non-variable upper part 1 of the SIMPSON file
                    'head_2' : '', # Non-variable upper part 2 of the SIMPSON file
                    'tail' : '', # Non-variable lower part of the SIMPSON file
                    }
            
            self.data = {
                    # Bools
                    'success' : False, # Turns to True if the optimisation was successful
                    'timer' : float(0), # Stores the time required for the simulation 
                    # Main
                    'tnopules' : int(0), # Total number of pulses of FAM-N
                    'sim_dur' : [], # Duration of each simulation in number of points
                    'max_amp' : [], # Amplitudes of FAM-N at the points of maximum amplitude
                    'max_pos' : [], # Position of the maximum
                    'inv_pos' : [], # Position of the pulse inversion
                    'max_pos_abs' : [], # Absolute position of the maxima
                    'inv_pos_abs' : [], # Absolute position of the phase inversions
                    'selected_full_density' : [], # List of selected density matrix for each pulse
                    'end_density' : [], # Stores the density matrix at the final maximum
                    'density_last_inversion' : [], # Stores the density matrix at the point of last inversion
                    'TIR' : float(0), # Theoretical improvement ratio of the current step
                    # Final exploitable result
                    'pulse_durations_nop' : [], # List of pulse durations in number of points
                    'pulse_durations_time' : [], # List of pulse durations in number of points
                    # In case there is a final optimisation for this step (see the method 'final_simulation')
                    'final_full_density' : [], # Final full density matrix
                    'final_select_density' : [], # Final detect density matrix
                    'final_max_amp': int(0), # Amplitude of maximum after re-doing the optimisation
                    'final_max_pos': int(0), # Relavite position of the maximum
                    'final_max_pos_abs': int(0), # Absolute position of the maximum
                    }
            
            self.temp = { # Temporary variables
                    # ------------------------------
                    # Temp out arrays
                    'simpson_out_full' : [], # Read-in SIMPSON file as a complex numpy array with full density matrix
                    'simpson_out_detect' : [], # Read-in SIMPSON file containing only the detect element
                    # ------------------------------
                    # 
                    'pulse_counter' : int(0), # Counts the number of pulses
                    'inversion_counter' : int(0), # Counts the number of inversions
                    'suffix' : '', #
                    'set_filename' : None, # Function that sets the filename on every step
                    'set_filename_final' : None, # Function that sets the filename on every step for the final step of a multiple steps optimisation
                    'cur_name' : '', # Stores the name of the file currently written         
                    # -------------------------------------------------------------
                    'initpoint' : int(0), # First point of the plot with all simulations
                    # 'pulse_phase' : True, # True means that the phase of the pulse is '+x', False means it is '-x' # No used !
                    'car_phase' : True, # Alternates between False and True as the direction of the pulse changes in duration
                    'count+' : int(0), # Positive counter
                    'count−' : int(0), # Negative counter
                    'current_saved_pos+' : int(0), # Saved position of the phase inversion for positive values 
                    'current_saved_pos−' : int(0), # Saved position of the phase inversion for negative values 
                    'current_max' : float(0), # Stores the current maximum intensity 
                    'current_pos' : int(0), # Stores the current position of the phase
                                # inversion that allows this maximum intensity
                    'sim_dur' : float(0), # Duration of current pulse
                    'previous_pulse_dur' : float(0) , # Duration of saved pulse
                    'previous_max' : float(0), # Contains the amplitude of the previous maximum
                    'previous_maxpos' : int(0), # Contains the position of the previous maximum
                    'previous_density_full' : None, # Full density matrix evolution of the inversion that gave the maximum signel
                    'previous_density_detect' : None, # Same, but only containing the detect element
                    'idx_maxpos' : int(0), # Temporary variable storing the position of the maxpos among all of the inflexion points
                    }
            
            self.figures = { # Contains figure classes
                    'opt' : None, # Figure
                    'opt_all_ax' : None, # Axe containing all steps of the optimisation
                    'opt_select_ax' : None, # Axe containing all steps of the optimisation
                    # -------------------------
                    'final' : None, # 
                    'final_ax': None, # Axe containing the final optimisation if applicable
                    } 
            
            # Set the function that defines the file name
            self.define_filename_fun()
            
            # Initialise figures if applicable
            if self.FAMN_main.param['debug_figure']:
                self.figures['opt'], (self.figures['opt_all_ax'],
                self.figures['opt_select_ax']) = plt.subplots(2,1, figsize = (10,10))
                
                # Plot containing all simulation steps
                self.figures['opt_all_ax'].set_title('All optimisation steps')
                self.figures['opt_all_ax'].set_xlabel('Time /μs')
                
                # Plot containing only plots selected in the FAM-N optimisation
                self.figures['opt_select_ax'].set_title('Selected optimisation steps')
                self.figures['opt_select_ax'].set_xlabel('Time /μs')
                
        # ---------------------------------------------------------------------
        # ----- Filename methods
        def fixed_filename(self, truc = ''):
            ''' Always return the same filename. 'truc' is dummy'''
            return self.FAMN_main.param['simpson_basename']
        
        
        def variable_filename(self):
            ''' Return a filename that depends on the input parameters '''
            return (self.FAMN_main.param['simpson_basename'] +
                                    '_step' + str(self.opt['order_index']+1) +
                                    '_pulse' + str(self.temp['pulse_counter']+1) +
                                    '_inv' + str(self.temp['inversion_counter']+1))
                                    
        def allcrystalline_filename(self, cryst):
            ''' Return a filename that depends on the input parameters including all crystal files'''
            return (self.FAMN_main.param['simpson_basename'] +
                                    '_step' + str(self.opt['order_index']+1) +
                                    '_pulse' + str(self.temp['pulse_counter']+1) +
                                    '_inv' + str(self.temp['inversion_counter']+1) +
                                    '_cryst' + str(cryst + 1))
        
        def finalise_filename(self):
            ''' Return a filename that depends on the input parameters '''
            return self.FAMN_main.param['simpson_basename'] +\
                                    '_final_step' + str(self.opt['order_index']+1) +\
                                    '_pulse' + str(self.temp['pulse_counter']+1)               
        
        def define_filename_fun(self):
            ''' Define how are filename called during the optimisation '''
            
            if self.FAMN_main.param['keep_simpson_files']: # Name changes at every iteration
                if self.FAMN_main.param['do_every_cryst']:
                    self.temp['set_filename'] = self.allcrystalline_filename
                else:
                    self.temp['set_filename'] = self.variable_filename
            else: # Name is kept constant
                self.temp['set_filename'] = self.fixed_filename
                
                
        # -----------------------------------------------------------------
        # ----- Write Simpson methods
        
        def write_crystal_file(self,noc):
            ''' Write the 'noc' crystal orientation into the crystal file from the
            list of crystalites '''
            cryst_file = open(self.FAMN_main.param['crystfile'] ,'w')
            # ---------------------------
            cryst_file.write('1\n')
            if self.FAMN_main.param['only_beta']:
                cryst_file.write(' 0 ' + str((180./numpy.pi)*self.FAMN_main.output['crystal_list'][1][noc]) + ' 1\n')
            else:
                cryst_file.write( str((180./numpy.pi)*self.FAMN_main.output['crystal_list'][0][noc]) +
                                          ' ' + str((180./numpy.pi)*self.FAMN_main.output['crystal_list'][1][noc]) +
                                          ' 1\n')
            # ---------------------------
            cryst_file.close()
            
            
        def write_crystal_file_all(self):
            ''' Write all crystal orientation into the crystal file from the list of
            crystallites '''
            cryst_file = open(self.FAMN_main.param['crystfile'] ,'w')
            # ------------------------------------------------------
            cryst_file.write(str(self.FAMN_main.output['no_cryst']) + '\n')
            weight = str(1/self.FAMN_main.output['no_cryst'])
            
            if self.FAMN_main.param['only_beta']:
                for truc in range(self.FAMN_main.output['no_cryst']):
                    cryst_file.write(' 0 ' +str((180./numpy.pi)*self.FAMN_main.output['crystal_list'][1][truc]) + ' ' + weight + '\n')
            else:
                for truc in range(self.FAMN_main.output['no_cryst']):
                    cryst_file.write( str((180./numpy.pi)*self.FAMN_main.output['crystal_list'][0][truc]) +
                                          ' ' + str((180./numpy.pi)*self.FAMN_main.output['crystal_list'][1][truc]) +
                                          ' ' + weight + '\n')
            # ------------------------------------------------------
            cryst_file.close()
        
        
        def write_simpson(self, density_matrix, nop, dim = 0): # , positive_sign = True
            ''' Write SIMPSON input file according to the input parameters
            density_matrix : initial density matrix of the current optimisation
            nop : number of points of the pulse program
            dim: dimension of the density matrix
            positive_sign: True <=> pulse is phase 'x', else '-x'
            suffix : string that come after the name (optionnal)
            '''
            simpson_file = open('./' + self.FAMN_main.param['simpson_tempfolder'] + '/' +\
                                         (self.temp['cur_name']  + '.in'),'w')
            simpson_file.write(self.data['head_1'])
            simpson_file.write(str(nop))
            simpson_file.write(self.data['head_2'])
            # --------------------------------------------------
            simpson_file.write('matrix set start list ')
            simpson_file.write(fmc.array_to_simpson(density_matrix, dim ))
            simpson_file.write('reset' + '\n')
            simpson_file.write('offset ' + str(self.FAMN_main.output['offset']) + '\n' )
            simpson_file.write('for {set i 0} {$i < ' + str(nop) + ' } {incr i 1} {' + '\n')
            simpson_file.write('pulse $par(tsw) $par(RF) ')
            # --------------------------------------------------
            # Change signe of the phase
            simpson_file.write('x\n')
#            if positive_sign:
#                simpson_file.write('x\n')
#            else:
#                simpson_file.write('-x\n')
            # -------------------------------------------------
            simpson_file.write(self.data['tail'])
            simpson_file.close()
            
        def execute_simpson(self, initial_density, sim_duration):
            ''' Write, execute and read all SIMPSON files for all crystallites'''
            
            # Set file name
            self.temp['cur_name'] = self.temp['set_filename']()
            
            # Write SIMPSON input file for first pulse
            self.write_simpson(initial_density, sim_duration,\
				   self.FAMN_main.output['dim_density']) # , self.temp['pulse_phase']

            # SIMPSON Execution
            if not self.FAMN_main.param['no_SIMPSON_execution']:
                subprocess.check_call('simpson ' + self.temp['cur_name'] + '.in',
                                          cwd= './' + self.FAMN_main.param['simpson_tempfolder'])
                #subprocess.call('simpson ' + './' + self.FAMN_main.param['simpson_tempfolder'] +\
            	#	  '/' + self.temp['cur_name'] + '.in')
                
            # Read in the file, split it into a text file with one line par section
            # Parse it into complex128 numpy array
            self.temp['simpson_out_full'] = fmc.readin_simpson_out(('./' + self.FAMN_main.param['simpson_tempfolder'] +\
            	  '/' + self.temp['cur_name'] + '.out'), sim_duration,
                self.FAMN_main.output['dim_density'],
            	self.FAMN_main.output['size_density'],
            	)
            
            # Add the initial density matrix in the list of detect density
            self.temp['simpson_out_full'] = numpy.concatenate((
            		[numpy.conj(initial_density)],self.temp['simpson_out_full']))
            
            # Get only the detect element from this numpy array
            self.temp['simpson_out_detect'] = numpy.real(fmc.get_all_detect_element(
                     self.temp['simpson_out_full'],
                     self.FAMN_main.output['list_detect_density'][self.opt['order_index']],
                     sim_duration+1))
            
            
        def execute_simpson_loop(self, initial_density, sim_duration):
            ''' Write, execute and read SIMPSON for crystallites one by one'''
            
            # Pre-allocate output list
            self.temp['simpson_out_full'] = numpy.zeros((self.FAMN_main.output['no_cryst'],
                            sim_duration + 1,
                            self.FAMN_main.output['dim_density'],
                            self.FAMN_main.output['dim_density']),
                            dtype = 'complex128')
            
            # -----------------------------------------------------------------
            # Loop over crystallites
            for truc in range(self.FAMN_main.output['no_cryst']):
                
                # Set file name
                self.temp['cur_name'] = self.temp['set_filename'](truc)
                
                # If first pulse: Get the initial density which is the same for every step
                if self.temp['pulse_counter'] == 0:
                    self.write_simpson(initial_density, sim_duration,\
    				   self.FAMN_main.output['dim_density']) # , self.temp['pulse_phase']
                else: # If not first pulse: write the SIMPSON file containing the previous density
                    self.write_simpson(initial_density[truc], sim_duration,\
    				   self.FAMN_main.output['dim_density'])
                
                # Write crystal file
                self.write_crystal_file(truc)
                
                # SIMPSON Execution
                if not self.FAMN_main.param['no_SIMPSON_execution']:
                    subprocess.check_call('simpson ' + self.temp['cur_name'] + '.in',
                                          cwd= './' + self.FAMN_main.param['simpson_tempfolder'])
                	# subprocess.check_call('simpson ' + './' + self.FAMN_main.param['simpson_tempfolder'] +\
                	#	  '/' + self.temp['cur_name'] + '.in')
                    
                # Read in the file, split it into a text file with one line par section
                # Parse it into complex128 numpy array
                self.temp['simpson_out_full'][truc][1:] = fmc.readin_simpson_out(('./' +
                                     self.FAMN_main.param['simpson_tempfolder'] +
                	  '/' + self.temp['cur_name'] + '.out'), sim_duration,
                    self.FAMN_main.output['dim_density'],
                	self.FAMN_main.output['size_density'])
                
                # Add the initial density matrix in the list of detect density 
                if self.temp['pulse_counter'] == 0:
                    self.temp['simpson_out_full'][truc][0] = numpy.conj(initial_density)
                else: 
                    self.temp['simpson_out_full'][truc][0] = numpy.conj(initial_density[truc])
            
            # -----------------------------------------------------------------
            # Get only the detect element from this numpy array
            self.temp['simpson_out_detect'] = numpy.real(fmc.get_all_detect_element(
                     numpy.average(self.temp['simpson_out_full'], axis = 0),
                     self.FAMN_main.output['list_detect_density'][self.opt['order_index']],
                     sim_duration+1))
            
            # -----------------------------------------------------------------
            
            
        # ---------------------------------------------------------------------
        # ----- Main methods
        
        def run_singlepulse(self):
            ''' Start the optimisation loop'''

            # Get suffix if the debug option 'keep_simpson_files' is set to True
            self.temp['pulse_counter'] = int(0)
            self.data['sim_dur'].append(self.opt['nop_CW'])
            
            # -----------------------------------------------------------------
            print('Simulation for pulse n°1')
            
            # Set intial phase positive
            self.temp['pulse_phase'] = True
            
            # Write SIMPSON file, execute and read-in for first pulse
            if self.FAMN_main.param['do_every_cryst']: # If every crystallite is on
                self.execute_simpson_loop(
                    initial_density = self.FAMN_main.output['list_initial_density'][self.opt['order_index']],
                    sim_duration = self.opt['nop_CW'])
            else: # If every crystallite is off
                self.write_crystal_file_all() # Print crystal file
                self.execute_simpson(
                    initial_density = self.FAMN_main.output['list_initial_density'][self.opt['order_index']],
                    sim_duration = self.opt['nop_CW'])
            
            # -----------------------------
            # Get the maximum position and amplitude of the single pulse
            self.temp['max_pos'] = fmc.get_maxima_and_positions(self.temp['simpson_out_detect'])
            
            # Makes sure that the final point is not located at the end of the current simulation.
            # Remove the point if it is.
            if self.temp['max_pos'][-1] == self.opt['nop_CW']:
                self.temp['max_pos'] = self.temp['max_pos'][:-1]

            # Pick the one with the maximum amplitude
            self.temp['max_amp'] = self.temp['simpson_out_detect'][self.temp['max_pos']]
            
            # Get the position of this max amplitude, get the position and report it to the final list
            self.data['max_amp'].append(numpy.max(self.temp['max_amp']))
            self.data['max_pos'].append(fmc.find_nearest_idx(self.temp['simpson_out_detect'], self.data['max_amp'][0]))
            
            # Initialise the absolute position arrays
            self.data['max_pos_abs'].append(self.data['max_pos'][0])
            
            # Save the density matrix for next steps
            self.data['selected_full_density'].append(self.temp['simpson_out_full'])
            
            # -----------------------------------------------------------------
            # Update figures
            if self.FAMN_main.param['debug_figure']:
                self.figures['opt_all_ax'].plot(1.0e6 * self.opt['time_increment'] * numpy.arange(0,self.opt['nop_CW']+1),
                        self.temp['simpson_out_detect'],
                        color = fmc.colour_in_plot(self.temp['pulse_counter']))
    
                
                self.figures['opt_select_ax'].plot(1.0e6*self.opt['time_increment'] * numpy.arange(0,self.opt['nop_CW']+1),
                         self.temp['simpson_out_detect'],
                        color = fmc.colour_in_plot(self.temp['pulse_counter']))
                self.figures['opt_select_ax'].axvline(x = self.data['max_pos_abs'][0]*1.0e6*self.opt['time_increment'],
                            color = 'red', linestyle = '--')
            # -----------------------------------------------------------------
    
    
        def one_optimisation_famn(self, current_point):
            ''' Run one optimisation corresponding to the current_point FAM-N'''
            print('      Inversion: ' + fmc.display_plus(current_point))
            
            # Get pulse duration
            self.temp['sim_dur'] = current_point + self.opt['nop_CW']
            
            # Write SIMPSON file, execute and read-in for first pulse
            if self.FAMN_main.param['do_every_cryst']:
                self.execute_simpson_loop( initial_density =\
                         self.data['selected_full_density'][self.temp['pulse_counter']-1][ 0:self.FAMN_main.output['no_cryst'],
    				    current_point + self.data['max_pos'][self.temp['pulse_counter']-1] , ...],
                       sim_duration = self.temp['sim_dur'])
            else:
                self.execute_simpson(initial_density = self.data['selected_full_density'][self.temp['pulse_counter']-1]\
    				   [current_point + self.data['max_pos'][self.temp['pulse_counter']-1]],
                       sim_duration = self.temp['sim_dur'])
            
            # -----------------------------
            # Get the maximum position and amplitude of the current simulation
            self.temp['max_pos'] = fmc.get_maxima_and_positions(self.temp['simpson_out_detect'])
            
            if self.temp['max_pos'].size: # If no inflexion point
                # Makes sure that the final point is not located at the end of the current simulation.
                # Remove the point if it is.
                if self.temp['max_pos'][-1] == self.opt['nop_CW']:
                    self.temp['max_pos'] = self.temp['max_pos'][:-1]

                # Pick the one with the maximum amplitude
                self.temp['max_amp'] = self.temp['simpson_out_detect'][self.temp['max_pos']]
                self.temp['idx_maxpos'] = fmc.find_nearest_idx(self.temp['max_amp'], max(self.temp['max_amp']))
                self.temp['max_pos'] = self.temp['max_pos'][self.temp['idx_maxpos']]
                self.temp['max_amp'] = self.temp['max_amp'][self.temp['idx_maxpos']]                
                
            else: # If no inflexion point
                self.temp['max_amp'] = max(self.temp['simpson_out_detect'])
                self.temp['max_pos'] = 0

            # -----------------------------------------------------------------
            # Update figure with all simulations
            
            if self.temp['pulse_counter'] >= 2:
                self.temp['initpoint'] = 1+self.data['inv_pos_abs'][self.temp['pulse_counter']-2] +\
                self.data['max_pos'][self.temp['pulse_counter']-1]
            else:
                self.temp['initpoint'] = self.data['max_pos'][self.temp['pulse_counter']-1]
                
            # ----------------------------------------
            if self.FAMN_main.param['debug_figure']:
                self.figures['opt_all_ax'].plot(1.0e6 * self.opt['time_increment'] *\
                         numpy.arange(current_point + self.temp['initpoint']  \
                                      , 1+2 * current_point + self.temp['initpoint'] + self.opt['nop_CW']),      
                                      numpy.real(self.temp['simpson_out_detect']),
                        color = fmc.colour_in_plot(self.temp['pulse_counter']))
            # -----------------------------------------------------------------
            
            
        # ---------------------------------------------------------------------
        def add_one_pulse(self):
            ''' Add one more pulse and continue the optimisation'''
            
            # Get suffix if the debug option 'keep_simpson_files' is set to True
            self.temp['pulse_counter'] += 1
            # self.temp['pulse_phase'] = not self.temp['pulse_phase'] # Inverse the phase sign
            
            print('Simulations for pulse n°' + str(self.temp['pulse_counter'] + 1))
            
            # Get the density after the first inversion
            self.temp['inversion_counter'] = 0
            self.one_optimisation_famn(0)
            
            self.temp['previous_max'] = self.temp['max_amp']
            self.temp['previous_maxpos'] = 0
            self.temp['previous_density_full'] = None
            self.temp['previous_density_detect'] = None
            self.temp['previous_pulse_dur'] = 0
            self.temp['current_max'] = 0
            self.temp['current_invpos'] = 0
            self.temp['current_saved_pos+'] = 0
            self.temp['current_saved_pos−'] = 0
            
            # Reinitialise_the_counters
            self.temp['count+'] = 0
            self.temp['count−'] = 0
            # self.temp['max_count+'] = self.data['max_pos'][self.temp['pulse_counter']-1] - 1
            # self.temp['max_count−'] = self.opt['nop_CW'] - self.data['max_pos'][self.temp['pulse_counter']-1] - 1
            self.temp['car_phase'] = False
            
            # -----------------------------------------
            while self.temp['count+'] < self.FAMN_main.param['tol_count+'] \
                or self.temp['count−'] < self.FAMN_main.param['tol_count−'] :
                
                self.temp['inversion_counter'] += 1
                    
                # Inverse the direction of the duration change if the 
                # corresponding counter was not reached, and add or remove
                if not self.temp['car_phase']:
                    if self.temp['count+'] < self.FAMN_main.param['tol_count+']:
                    # Positive duration evolution
                        self.temp['car_phase'] = True
                        self.temp['current_invpos'] = self.temp['current_saved_pos+'] + 1
                    else:
                        # If counter reached for positive values
                        self.temp['current_invpos'] += -1
                    
                elif self.temp['car_phase']:
                    if self.temp['count−'] < self.FAMN_main.param['tol_count−']:
                    # Negative duration evolution
                        self.temp['car_phase'] = False
                        self.temp['current_invpos'] = self.temp['current_saved_pos−'] - 1
                    else:
                        # If counter reached for negative values
                        self.temp['current_invpos'] += 1
                
                # -------------------
                # Perform the optimisation
                self.one_optimisation_famn(self.temp['current_invpos'])
                self.temp['current_max'] = self.temp['max_amp']
                
                # Compare the old result and the new result
                if self.temp['current_max'] > self.temp['previous_max']:
                    self.temp['previous_max'] = self.temp['current_max']
                    self.temp['previous_maxpos'] = self.temp['max_pos']
                    self.temp['previous_invpos'] = self.temp['current_invpos']
                    self.temp['previous_density_full'] = self.temp['simpson_out_full']
                    self.temp['previous_density_detect'] = self.temp['simpson_out_detect']
                    self.temp['previous_pulse_dur'] = self.temp['sim_dur']
                    # Reset the counter if signal imporve
                    if self.temp['car_phase']:
                        self.temp['count+'] = 0
                    else:
                        self.temp['count−'] = 0
                else:
                    # Add one to the counter if no signal improve
                    if self.temp['car_phase']:
                        self.temp['count+'] += 1
                    else:
                        self.temp['count−'] += 1
            
                # Save the current point of phase inversion
                if self.temp['car_phase']:
                    self.temp['current_saved_pos+'] = self.temp['current_invpos']
                else:
                    self.temp['current_saved_pos−'] = self.temp['current_invpos']
            
                # Security checks: verify that the next point do not go over the point limit
                # or inside the previous pulse
                
                if self.temp['current_saved_pos+'] +1 == self.data['sim_dur'][self.temp['pulse_counter']-1] -\
                            self.data['max_pos'][self.temp['pulse_counter']-1]:
                    print('Warning: end of current pulse reached.')
                    self.temp['count+'] = self.FAMN_main.param['tol_count+']
                    
                if self.temp['current_saved_pos−'] +1 == self.data['max_pos'][self.temp['pulse_counter']-1] :
                    print('Warning: beginning of current pulse reached.')
                    self.temp['count−'] = self.FAMN_main.param['tol_count−']


        def optimise_famn(self):
            ''' Run the main optimisation loop of FAM-N'''
            
            # Run until the maximum allowed number of pulses is reached
            while self.temp['pulse_counter'] < self.opt['max_nopules']:
            
                # Add an additionnal pulse
                self.add_one_pulse()
                
                # Verify if the additionnal pulse allows a signal imporvement
                if self.temp['previous_max'] > self.FAMN_main.param['min_improvement'] *\
                            self.data['max_amp'][self.temp['pulse_counter']-1]:
                    
                    # Copy the current maximum into the stored data
                    self.data['max_amp'].append(self.temp['previous_max'])
                    self.data['max_pos'].append(self.temp['previous_maxpos'])
                    self.data['inv_pos'].append(self.temp['previous_invpos'])
                    self.data['sim_dur'].append(self.temp['previous_pulse_dur'])
                    self.data['selected_full_density'].append(self.temp['previous_density_full'])
                    
                    # Update the list containing the absolute positions of maxima and phase inversions
                    self.data['inv_pos_abs'].append(
                        self.data['inv_pos'][self.temp['pulse_counter']-1] +\
                        self.data['max_pos_abs'][self.temp['pulse_counter']-1])
    
                    self.data['max_pos_abs'].append(
                            self.data['inv_pos_abs'][self.temp['pulse_counter']-1] +\
                            self.data['max_pos'][self.temp['pulse_counter']])
                    
                    # Display the improvement for the current pulse on the console
                    self.data['TIR'] = self.data['max_amp'][self.temp['pulse_counter']]/\
                                       self.data['max_amp'][self.temp['pulse_counter']-1]
                    print('      Improvement: ×' + str(numpy.round(self.data['TIR'],2)))
                    
                    # Update the selected_simulation figure
                    if self.FAMN_main.param['debug_figure']:
                        self.figures['opt_select_ax'].plot(1.0e6 * self.opt['time_increment'] *\
                                 numpy.arange(self.data['inv_pos_abs'][self.temp['pulse_counter']-1],
                                              self.data['inv_pos_abs'][self.temp['pulse_counter']-1] +\
                                              self.data['sim_dur'][self.temp['pulse_counter']]+1,
                                              ),
                                self.temp['previous_density_detect'],
                                color = fmc.colour_in_plot(self.temp['pulse_counter']),)
                        # Plot vertical line at phase inversion
                        self.figures['opt_select_ax'].axvline(x = (self.data['inv_pos_abs'][self.temp['pulse_counter']-1])*\
                                    1.0e6*self.opt['time_increment'], color = 'black') 
                        # Plot vertical line at maximum
                        self.figures['opt_select_ax'].axvline(x = self.data['max_pos_abs'][self.temp['pulse_counter']]*\
                                    1.0e6*self.opt['time_increment'], color = 'red', linestyle = '--') 
                    
                else:
                    # Break if no further improvement
                    break
                
            # -----------------------------------------------------------------
            # Write output data
            self.data['tnopules'] = int(self.temp['pulse_counter'])
            
            print('Theoretical improvement ratio: ×'\
                  + str(numpy.round(self.data['max_amp'][self.data['tnopules']-1]/self.data['max_amp'][0],2)))
            
            # Write the density matrix at the end of the execution
            if self.FAMN_main.param['do_every_cryst']: # If every crystallite is on
                self.data['end_density'] = self.data['selected_full_density']\
                [self.data['tnopules']-1][0:self.FAMN_main.output['no_cryst'],
                 self.data['max_pos'][self.data['tnopules']-1], ...]                
            else: # If every crystallite is off
                self.data['end_density'] = self.data['selected_full_density']\
                [self.data['tnopules']-1][self.data['max_pos'][self.data['tnopules']-1]]
            
            # Write the density at the end of the final pulse
            if self.data['tnopules'] > 1: # If more than one pulse
                if self.FAMN_main.param['do_every_cryst']: # If every crystallite is on
                    self.data['density_last_inversion'] = self.data['selected_full_density']\
                    [self.data['tnopules']-1][0:self.FAMN_main.output['no_cryst'],
                     self.data['inv_pos'][self.data['tnopules']-2], ...]
                else:
                    self.data['density_last_inversion'] = self.data['selected_full_density']\
                    [self.data['tnopules']-1][self.data['inv_pos'][self.data['tnopules']-2]]
            else:
                # If the optimisation completed with only one pulse, just export the initial density matrix
                self.data['density_last_inversion'] =\
                self.FAMN_main.output['list_initial_density'][self.opt['order_index']-1]
            
            # Remove input files if required
            if not self.FAMN_main.param['keep_simpson_files']:
                os.remove('./' + self.FAMN_main.param['simpson_tempfolder'] + '/' +\
                                         (self.temp['cur_name']  + '.in'))
                os.remove('./' + self.FAMN_main.param['simpson_tempfolder'] + '/' +\
                                         (self.temp['cur_name']  + '.out'))
                
        # ---------------------------------------------------------------------
        # ----- Finalise methods
        
        def final_simulation(self):
            ''' Method executed only in case of multiple-step FAM-N.
            Re-simulate the optimised pulse of this step with the final density
            matrix of the previous pulse if applicable. See output['do_final_sim']
            in the class 'FAMN_main'.
            '''
            
            # Initialise some arrays
            self.data['final_full_density'] = []
            self.data['final_select_density'] = []
            
            # Set the function that defines the file name
            if self.FAMN_main.param['keep_simpson_files']: # Name changes at every iteration
                self.temp['set_filename_final'] = self.finalise_filename
            else: # Name is kept constant
                self.temp['set_filename_final'] = self.fixed_filename
            
            # Do the simulations for the single pulse
            print('Final simulation for pulse n°1')
            self.temp['pulse_counter'] = 0
            self.temp['cur_name'] = self.temp['set_filename_final']()
             
            # Write SIMPSON file, execute and read-in for first pulse
            self.execute_simpson(initial_density = self.FAMN_main.output['FAMNsteps']\
                                 [self.opt['order_index']-1].data['end_density']  ,
                   sim_duration = self.opt['nop_CW'] )

            # Save the current density evolution and add the initial density matrix
            # in the list of detect density
            self.data['final_select_density'].append(fmc.get_all_detect_element(
                    self.data['final_full_density'][0],
                    self.FAMN_main.output['list_detect_density'][self.opt['order_index']],
                    1+self.opt['nop_CW']))
            
            # Test figure
            if self.FAMN_main.param['debug_figure']: 
                self.figures['final'], self.figures['final_ax'] = plt.subplots(1,1, figsize = (10,5))
                self.figures['final_ax'].set_title('Final optimisation for step n°' + str(1+ self.opt['order_index']) )
                self.figures['final_ax'].set_xlabel('Time /μs')
                self.figures['final_ax'].plot(1.0e6 * self.opt['time_increment'] * numpy.arange(0,self.opt['nop_CW']+1),
                         numpy.real(self.data['final_select_density'][0]),
                         color = fmc.colour_in_plot(0),
                         )
            
            # -----------------------------------------------------------------
            # Loop over the number of pulses
            for self.temp['pulse_counter'] in range(1,self.data['tnopules']):
                
                print('Final simulation for pulse n°' + str(self.temp['pulse_counter']+1))
                # self.temp['pulse_phase'] = not self.temp['pulse_phase'] # Inverse phase # Not used
                
                # Change filename
                self.temp['cur_name'] = self.temp['set_filename_final']()
                
                # Write SIMPSON file, execute and read-in for first pulse
                self.execute_simpson(initial_density = self.data['final_full_density'][self.temp['pulse_counter']-1]\
                                   [self.data['max_pos'][self.temp['pulse_counter']-1]+\
                                    self.data['inv_pos'][self.temp['pulse_counter']-1]] ,
                       sim_duration = self.data['sim_dur'][self.temp['pulse_counter']] )
                
                # Save the current density evolution
                self.data['final_select_density'].append(fmc.get_all_detect_element(\
                         self.data['final_full_density'][self.temp['pulse_counter']],
                         self.FAMN_main.output['list_detect_density'][self.opt['order_index']],
                         1+self.data['sim_dur'][self.temp['pulse_counter']]))
            
                # Test figure
                if self.FAMN_main.param['debug_figure']:
                    self.figures['final_ax'].plot(1.0e6 * self.opt['time_increment'] *\
                                numpy.arange(self.data['inv_pos_abs'][self.temp['pulse_counter']-1],
                                          self.data['inv_pos_abs'][self.temp['pulse_counter']-1] +\
                                          self.data['sim_dur'][ self.temp['pulse_counter']] + 1),
                                            numpy.real(self.data['final_select_density'][self.temp['pulse_counter']]),
                                            color = fmc.colour_in_plot(self.temp['pulse_counter']),
                                            )
                    self.figures['final_ax'].axvline(x =  1.0e6 * self.opt['time_increment'] *\
                                (self.data['inv_pos_abs'][self.temp['pulse_counter']-1]),
                                color = 'black')
                            
            # -----------------------------------------------------------------
            # Get positions of maxima
            self.data['final_max_amp'] = numpy.max(self.data['final_select_density'][self.data['tnopules']-1])
            self.data['final_max_pos'] = fmc.find_nearest_idx(\
                     self.data['final_select_density'][self.data['tnopules']-1],
                     self.data['final_max_amp'])
            self.data['final_max_pos_abs'] = self.data['final_max_pos'] + self.data['inv_pos_abs'][self.data['tnopules']-2]
            
            # Plot vertical line at maximum
            if self.FAMN_main.param['debug_figure']:
                self.figures['final_ax'].axvline(x = self.data['final_max_pos_abs']*\
                            1.0e6*self.opt['time_increment'], color = 'red', linestyle = '--')
            
            # Remove the final file if required
            if not self.FAMN_main.param['keep_simpson_files']:
                # Remove input files if required
                os.remove('./' + self.FAMN_main.param['simpson_tempfolder'] + '/' +\
                                         (self.temp['cur_name']  + '.in'))
                os.remove('./' + self.FAMN_main.param['simpson_tempfolder'] + '/' +\
                                         (self.temp['cur_name']  + '.out'))
                os.remove(self.FAMN_main.param['crystfile'])
                
            
        def optimisation_complete(self):
            ''' Do a couple of things at the end of the optimisation if successful'''
            print('')
            print('FAM-N optimisation completed successfully')
            print('')
            
            # Get the pulse durations
            if self.FAMN_main.output['do_final_sim'][self.opt['order_index']-1]:
                # If a final optimisation was carried out
                temp = self.data['max_pos']
                temp[self.data['tnopules']-1] = self.data['final_max_pos']
                self.data['pulse_durations_nop'] = fmc.get_pulse_duration(
                        self.data['max_pos'],self.data['inv_pos'],self.data['tnopules'])
            else: # If not
                self.data['pulse_durations_nop'] = fmc.get_pulse_duration(
                        self.data['max_pos'],self.data['inv_pos'],self.data['tnopules'])
            
            # Get the durations in μs
            self.data['pulse_durations_time'] = self.opt['time_increment'] * self.data['pulse_durations_nop']
            
            # End
            # self.temp = {}
            self.data['success'] = True
            # -------------------------------
            
         # ----- Figure methods
        def save_opt_figure(self, filename, filepath = ''):
            ''' Save the optimisation figure '''
            
            if filepath == '':
                path = fmc.change_ext(filename,'pdf')
            else:
                path = os.path.normpath(filepath) + os.sep + fmc.change_ext(filename,'pdf')
            
            self.figures['opt'].figsave(path)
            
            
        def save_final_figure(self, filename, filepath = ''):
            ''' Save the optimisation figure '''

            if filepath == '':
                path = fmc.change_ext(filename,'pdf')
            else:
                path = os.path.normpath(filepath) + os.sep + fmc.change_ext(filename,'pdf')
            
            if self.figures['final'] != None :
                self.figures['final'].figsave(path)
                
    # -------------------------------------------------------------------------
    # End FAM-N step
    
    
    
    
    
    # -------------------------------------------------------------------------
    # ----- Utility functions
    
    def savestate(self, filename, filepath = ''):
        ''' Save the content of current instance into a filename located in filepath'''
        
        if filepath == '':
            outputfile = open(fmc.change_ext(filename,'pkl') ,'wb')
        else:
            outputfile = open(os.path.normpath(filepath) + os.sep + fmc.change_ext(filename,'pkl') ,'wb')
        # ------------------------------------
        pickle.dump(self.__dict__, outputfile, protocol=pickle.HIGHEST_PROTOCOL)
        outputfile.close()
        
        
    def loadstate(self, filename = '__previous_session__.pkl', filepath = ''):
        ''' Save the content of an instance into a filename located in filepath'''
        if filepath == '':
            inputfile = open(fmc.change_ext(filename,'pkl') ,'rb')
        else:
            inputfile = open(os.path.normpath(filepath) + os.sep + fmc.change_ext(filename,'pkl') ,'rb')
        # -------------------------------
        try:
            loadedfile = pickle.load(inputfile)
            inputfile.close()
            # Verifies that the loaded file correspond to the sort of instance
            if loadedfile.parameters['klupolnghsietd'] == 'vdrffuadbpt':
                self.__dict__ = loadedfile
                loadedfile.output['loading_error'] = False
            else:
                loadedfile.output['loading_error'] = True
                if self.param['console_warning']:
                    raise IOError('Loaded file not valid.')
        except:
            pass
        
    # -------------------------------------------------------------------------
    # ----- Main functions
    
    def initialise(self):
        ''' Fill in the FAMN_main class with the input parameters, check for errors etc. '''
        print('')
        print('--- FAM-N optimisation ---')
        print('')
        print('Initialising…')
        print('')
        
        # Verifies if SIMPSON is installed
        self.param['SIMPSON_installed'] = fmc.is_simpson_installed()
        
        if (not self.param['SIMPSON_installed']
            and not self.param['no_SIMPSON_execution']):
            if self.param['console_warning']:
                raise ImportError('SIMPSON is not intalled on this system')
            else:
                ImportWarning('Execution without SIMPSON installed')
                
        else: # If SIMPSON is installed
            
            # Set the current optimisation as "Not completed"
            self.output['completed'] = False
            
            # Start timer
            self.output['total_execution'] = fmc.start_timer()
            
            # Make temporary directory if not exist
            if not os.path.exists(self.param['simpson_tempfolder']):
                os.makedirs(self.param['simpson_tempfolder'])
                
            # Turn some variables into integers
            self.param['no_cryst'] = int(self.param['no_cryst'])
            self.param['no_γ'] = int(self.param['no_γ'])
            self.param['max_nopulses'] = int(self.param['max_nopulses'])
            self.param['nop_360'] = int(self.param['nop_360'])
            self.param['num_cores'] = int(self.param['num_cores'])
            
            # Try to read spin info from xmd file
            try:
                self.param['spin'] = fmc.read_spin_from_xml(fmc.chemical_symbol_first(self.param['nucleus']))
            except IndexError:
                # If non-existing nucleus
                self.output['non_existing_nucleus'] = True
                if self.param['console_warning']:
                    raise IndexError('Nucleus "' + self.param['nucleus'] + '" does not exist.')
                    
            else: # If existing nucleus
                self.output['non_existing_nucleus'] = False
                
                # Verifies if half-integer quadrupolar nucleus
                if self.param['spin'] not in ('3/2','5/2','7/2','9/2'):
                    self.output['non_quad_spin'] = True
                    if self.param['console_warning']:
                        raise ValueError('The nucleus ' + str(self.param['nucleus']) +\
                         ' is not a half-integer quadrupolar nucleus (I = ' + self.param['spin'] + ')')
                        
                else: # If half-integer quadrupolar nucleus
                    self.output['non_quad_spin'] = False
                    
                    try: # Verifies if the number of crystallite is integer
                        self.param['no_cryst'] = int(self.param['no_cryst'])
                    except:
                        self.output['crystallites_not_integer'] = True
                        if self.param['console_warning']:
                            raise ValueError('Number of crystallites provided ( ' + 
                                                  str(self.param['no_cryst']) + ' ) is not integer')
                            
                    else:
                        # Value given for crystallites is integer
                        self.output['crystallites_not_integer'] = False
                        
                        if self.param['no_cryst'] < 2: 
                            self.output['not_enough_crystallites'] = True
                            if self.param['console_warning']:
                                raise ValueError('Not enough crystallites ( ' +
                                                str(self.param['no_cryst']) + ' ) were given')
                                
                        else: # If enough crystallites
                            # ---------------------------------------------------------
                            # Make most of the preparation and verifications
                            
                            # Make sure the export filepath is a path, and empty it if not
                            try:
                                self.param['export_filepath'] = os.path.normpath(self.param['export_filepath'])
                            except:
                                self.param['export_filepath'] = os.path.normpath('')
                                
                            # Get Larmor frequency for current nucleus
                            self.output['Larmor_frequency'] = fmc.get_Larmor_frequency(\
                                       self.param['nucleus'],
                                       proton_larmor_frequency = self.param['proton_larmor_frequency'],
                                       )
                            
                            # Use only beta Euler angles in the assymetry parameter is low
                            if self.param['ηQ'] < 0.01:
                                self.param['only_beta'] = True
                                # Use only the square-root number of crystallites if only beta is used
                                self.output['no_cryst'] = int(numpy.round(numpy.sqrt(self.param['no_cryst'])))
                            else:
                                self.param['only_beta'] = False
                                self.output['no_cryst'] = int(self.param['no_cryst'])
                            
                            # Generate the crystal file
                            if self.param['cryst_gen_method'] == 'golden':
                                self.output['crystal_list'] = fmc.golden_sphere(self.output['no_cryst'])
                            elif self.param['cryst_gen_method'] == 'zcw':
                                self.output['crystal_list'] = fmc.zcw_angles(self.output['no_cryst'])
                            elif self.param['cryst_gen_method'] == 'random':
                                self.output['crystal_list'] = fmc.random_sphere_generator(self.output['no_cryst'])
                                
                            # Set the crystal file in for SIMPSON
                            self.param['crystfile'] = self.param['cryst_name'] 
                            #self.param['crystfile'] = ( './' + self.param['simpson_tempfolder'] +\
                            #               '/' + self.param['cryst_name'] )
                                
                            # Set the offset to QIS (Quadrupolar induced shift) if applicable, to 0 ejse
                            if self.param['set_offset_to_qis']:
                                self.output['offset'] = fmc.get_QIS(self.param['nucleus'],
                                    float(eval(self.param['spin'])),
                                    self.output['Larmor_frequency'],
                                    self.param['CQ'],
                                    self.param['ηQ'],
                                    )
                            else:
                                self.output['offset'] = 0
                            
                            # Set the size and dimension of the density matrix
                            self.output['dim_density'] = fmc.get_dim_from_spin(self.param['spin'])
                            self.output['size_density'] = self.output['dim_density']**2
                            
                            # Reinitialise error bools
                            self.output['partial_infos_complex'] = False
                            self.output['complex_parse_error'] = False
                            self.output['Non_consistent_number_of_famn_steps'] = False
                            self.output['nonexisting_transition'] = False
                            self.output['impossible_transition'] = False

                            # Get the required transition
                            # ---------------------------------------------------------
                            if ((self.param['complex_optimisation'] == '' and
                                self.param['single_pulse_nutation_list'] != '')
                                or (self.param['complex_optimisation'] != '' and
                                self.param['single_pulse_nutation_list'] == '')):
                                # If the input for complex optimisation is partial
                                
                                self.output['partial_infos_complex'] = True
                                if self.param['console_warning']:
                                    raise ValueError('Partial informations for informations for'
                                                     ' multiple FAM-N were provided:\n'
                                                     "parameters['complex_optimisation'] = " +
                                                     str(self.param['complex_optimisation']) + '\n' +
                                                     "parameters['single_pulse_nutation_list'] = " + 
                                                     str(self.param['single_pulse_nutation_list']) + '\n' +
                                                     'Please provide both information or empty both fields'
                                                     ' to carry out simple FAM-N optimisation.'
                                                     )
                            
                            # ---------------------------------------------------------
                            elif (self.param['complex_optimisation'] != '' and
                                self.param['single_pulse_nutation_list'] != '') : # Case the required optimisation is complex
                                
                                # Read the input string
                                self.output['list_single_pulse_nutation'] =\
                                    numpy.array([eval(self.param['single_pulse_nutation_list'])])
                                self.output['list_initial_density'], self.output['list_detect_density'] =\
                                    fmc.parse_complex_famn(self.param['complex_optimisation'],
                                                       self.output['dim_density'])
                                
                                # Get the number of steps
                                self.output['no_FAMN_steps'] = len(self.output['list_initial_density'])
                                
                                # Verify if no array is return with a 'nan',
                                # which means that there was a problem with parsing.
                                for truc in range(self.output['no_FAMN_steps']):
                                    try:
                                        if (numpy.any(numpy.isnan(self.output['list_initial_density'][truc])) or
                                            numpy.any(numpy.isnan(self.output['list_detect_density'][truc]))):
                                            self.output['complex_parse_error'] = True
                                            if self.param['console_warning']:
                                                raise ValueError('Error parsing the input for '
                                                    'complex FAM-N.\nPlease verify the input string, '
                                                    'if all transitions are possible for the given spin etc…')
                                    except TypeError:
                                        pass
                                        
                                        # self.output['list_initial_density'][truc].any() == numpy.nan
                                
                                if len(self.output['list_single_pulse_nutation']) == 1:
                                    # Fill the nutation list if only one value was added
                                    self.output['list_single_pulse_nutation'] =\
                                        numpy.full(self.output['no_FAMN_steps'],
                                        self.output['list_single_pulse_nutation'][0]
                                            )
                                
                                # Return an error if the number of steps is not consistent between
                                # the number of FAM-N steps and the number of nutation rates        
                                elif len(self.output['list_single_pulse_nutation']) != self.output['no_FAMN_steps']:
                                    self.output['Non_consistent_number_of_famn_steps'] = True
                                    if self.param['console_warning']:
                                                raise ValueError('Number of steps of this '
                                                    'optimisation (' + str(self.param['complex_optimisation']) + ' )'
                                                    'not consistent with the numbers of nutation rates ( ' +\
                                                    str(self.param['single_pulse_nutation_list']) + ' )')            
                                                
                                # ------------------------------
                                self.output['simple_opt'] = False
                                self.output['transfer_type'] = 'unknown'
                                
                            # ---------------------------------------------------------
                            else: # Case of a simple optimisation, i.e. one-step conversion pulse 
                                self.output['no_FAMN_steps'] = 1
                                self.output['transfer_type'] = 'conversion'
                                self.output['simple_opt'] = True
                                # --------------------------------
                                # Nb : all matrix elements are inverted relative to the convention
                                # of SIMPSON for detected elements
                                if self.param['transition'] == 3:
                                    if self.param['spin'] == '3/2':
                                        self.output['list_single_pulse_nutation'] = [120]
                                        self.output['list_initial_density'] =\
                                            [fmc.one_element_density(self.output['dim_density'],[1,4])]
                                        self.output['list_detect_density'] =\
                                            [fmc.one_element_density(self.output['dim_density'],[3,2])]
                                        
                                    if self.param['spin'] == '5/2':
                                        self.output['list_single_pulse_nutation'] = [75]
                                        self.output['list_initial_density'] =\
                                            [fmc.one_element_density(self.output['dim_density'],[2,5])]
                                        self.output['list_detect_density'] =\
                                            [fmc.one_element_density(self.output['dim_density'],[4,3])]
                                            
                                    if self.param['spin'] == '7/2':
                                        self.output['list_single_pulse_nutation'] = [55]
                                        self.output['list_initial_density'] =\
                                            [fmc.one_element_density(self.output['dim_density'],[3,6])]
                                        self.output['list_detect_density'] =\
                                            [fmc.one_element_density(self.output['dim_density'],[5,4])]
                                        
                                    if self.param['spin'] == '9/2':
                                        self.output['list_single_pulse_nutation'] = [40]
                                        self.output['list_initial_density'] =\
                                            [fmc.one_element_density(self.output['dim_density'],[4,7])]
                                        self.output['list_detect_density'] =\
                                            [fmc.one_element_density(self.output['dim_density'],[6,5])]
                                        
                                    
                                # --------------------------------
                                elif self.param['transition'] == 5:
                                    if self.param['spin'] == '3/2':
                                        self.output['impossible_transition'] = True
                                        
                                    if self.param['spin'] == '5/2':
                                        self.output['list_single_pulse_nutation'] = [120]
                                        self.output['list_initial_density'] =\
                                            [fmc.one_element_density(self.output['dim_density'],[1,6])]
                                        self.output['list_detect_density'] =\
                                            [fmc.one_element_density(self.output['dim_density'],[4,3])]
                                        
                                    if self.param['spin'] == '7/2':
                                        self.output['list_single_pulse_nutation'] = [75]
                                        self.output['list_initial_density'] =\
                                            [fmc.one_element_density(self.output['dim_density'],[2,7])]
                                        self.output['list_detect_density'] =\
                                            [fmc.one_element_density(self.output['dim_density'],[5,4])]
                                            
                                    if self.param['spin'] == '9/2':
                                        self.output['list_single_pulse_nutation'] = [55]
                                        self.output['list_initial_density'] =\
                                            [fmc.one_element_density(self.output['dim_density'],[3,8])]
                                        self.output['list_detect_density'] =\
                                            [fmc.one_element_density(self.output['dim_density'],[6,5])]
                                        
                                # --------------------------------
                                elif self.param['transition'] == 7:
                                    if self.param['spin'] == ('3/2' or '5/2'):
                                        self.output['impossible_transition'] = True
                                        
                                    if self.param['spin'] == '7/2':
                                        self.output['list_single_pulse_nutation'] = [120]
                                        self.output['list_initial_density'] =\
                                            [fmc.one_element_density(self.output['dim_density'],[1,8])]
                                        self.output['list_detect_density'] =\
                                            [fmc.one_element_density(self.output['dim_density'],[5,4])]
                                            
                                    if self.param['spin'] == '9/2':
                                        self.output['list_single_pulse_nutation'] = [75]
                                        self.output['list_initial_density'] =\
                                            [fmc.one_element_density(self.output['dim_density'],[2,9])]
                                        self.output['list_detect_density'] =\
                                            [fmc.one_element_density(self.output['dim_density'],[6,5])]
                                        
                                # --------------------------------
                                elif self.param['transition'] == 9:
                                    if self.param['spin'] == ('3/2' or '5/2' or '7/2'):
                                        self.output['impossible_transition'] = True
                                        
                                    if self.param['spin'] == '9/2':
                                        self.output['list_single_pulse_nutation'] = [120]
                                        self.output['list_initial_density'] =\
                                            [fmc.one_element_density(self.output['dim_density'],[1,4])]
                                        self.output['list_detect_density'] =\
                                            [fmc.one_element_density(self.output['dim_density'],[6,5])]
                                            
                                # --------------------------------
                                else:
                                    self.output['nonexisting_transition'] = True
                                    
                            # Returns error messages if applicable
                                if self.param['console_warning']:
                                    if self.output['nonexisting_transition']:
                                        raise ValueError('"' + str(self.param['transition']) +\
                                                         '" is not a valid transition.' )
                                    if self.output['impossible_transition']:
                                        raise ValueError('The transition ' + str(self.param['transition']) +\
                                                         ' is not possible for spin ' + self.param['spin'])
                                        
                        # Pre-allocate the array that determines if an optimisation for the final step must be redone
                        # or if there is a succession in subsequent pulses
                        self.output['is_succesion'] = numpy.zeros(self.output['no_FAMN_steps'], dtype = 'bool')
                        self.output['do_final_sim'] = numpy.zeros(self.output['no_FAMN_steps'], dtype = 'bool')
                
                # -------------------------------------------------------------------------------------------
        
    
    # ------------------------------------------------------------------------
    def initialise_step(self,stepcnt):
        ''' Initialise the class of a FAM-N optimisation
        "stepcnt" is the index of the step to be initialised
        '''

        # Create an instance of a FAM-N step optimisation
        self.FAMNsteps.append(self.FAMN_step(self))
        
        # Start timer
        self.FAMNsteps[stepcnt].data['timer'] = fmc.start_timer()
        
        # Set the order in the total FAM-N sequence of the newly introduced step
        self.FAMNsteps[stepcnt].opt['order_index'] = stepcnt

        # Initialise the optimisation
        self.FAMNsteps[stepcnt].opt['detect_density'] =\
            self.output['list_detect_density'][stepcnt]
        
        # Get maxdt
        self.FAMNsteps[stepcnt].data['maxdt'] =\
            1.0e6/(self.param['maxdt_ratio']*self.param['MAS'])
        
        self.temp['nop_calc'] = int(numpy.round((self.param['nop_360'])*\
            (self.output['list_single_pulse_nutation'][stepcnt]/360.)))
        
        # Compare the calculated time increment with the minimum time increment as 
        # specified in the input.

        if  self.param['min_time_increment'] != 0:
            self.temp['nop_min'] = int(numpy.round(self.output['list_single_pulse_nutation'][stepcnt]/\
                     (360.*self.param['RF']*self.param['min_time_increment'])))
            # Set the time increment to the minimum value between the parameter
            # famn.parameters['min_time_increment'] and the number of points self.param['nop_360']
            if self.temp['nop_min'] < self.temp['nop_calc']:
                # -------------------------------------------------------------
                self.output['using_min_timeinc'] = False
                # If time increment is below the minimum allowed
                # -------------------------------------------------------------
                self.FAMNsteps[stepcnt].opt['time_increment'] =\
                self.param['min_time_increment']
                
                # -------------------------------------------------------------
                self.FAMNsteps[stepcnt].opt['nop_CW'] =\
                int(numpy.round(((self.output['list_single_pulse_nutation'][stepcnt])/\
                    (360*self.param['RF']*self.param['min_time_increment']))))
            else:
                # If time increment is above the minimum
                self.FAMNsteps[stepcnt].opt['nop_CW'] =\
                    self.temp['nop_calc']
                self.FAMNsteps[stepcnt].opt['time_increment'] =\
                    (1./(self.param['RF']*self.param['nop_360']))
        
        else:
            # If the minimum time increment is not zero
            self.FAMNsteps[stepcnt].opt['nop_CW'] =\
                self.temp['nop_calc']
            self.FAMNsteps[stepcnt].opt['time_increment'] =\
                (1./(self.param['RF']*self.param['nop_360']))
            
        # Set maxdt
        self.FAMNsteps[stepcnt].data['maxdt'] =\
        1. / (self.param['MAS'] * self.param['maxdt_ratio'])
            
        # -----------------------------------------------------------
        # Set the invariable parts of the SIMPSON file
        # Head part 1:
        self.FAMNsteps[stepcnt].data['head_1'] =\
            ('spinsys {' + '\n'
            'channels ' + self.param['nucleus'] + '\n'
            'nuclei '  + self.param['nucleus'] + '\n'
            'quadrupole  1 2 ' + str(self.param['CQ']) + ' ' + str(self.param['ηQ']) +\
            ' 0 0 0  ' + '\n'
            '}' + '\n'
            'par {' + '\n'
            'num_cores ' + str(self.param['num_cores']) + '\n'
            'spin_rate ' + str(self.param['MAS']) + '\n'
            'method ' + str(self.param['opt_method']) + '\n'
            'variable tsw ' + str(1.0e6 * self.FAMNsteps[stepcnt].opt['time_increment']) + '\n'
            'variable RF ' + str(self.param['RF']) + '\n'
            'sw 1.0e6/tsw' + '\n'
            'crystal_file ' + self.param['crystfile'] + '\n'
            'gamma_angles ' + str(self.param['no_γ']) + '\n'
            'proton_frequency ' + str(self.param['proton_larmor_frequency']) + '\n'
            'verbose 0' + '\n'
            'np ' + str( self.output['size_density'] ) + '*'
            )
            
        # -----------------------------------------------------------
        # Head part 2:
        self.FAMNsteps[stepcnt].data['head_2'] =\
            ( '\n} \n'
            'proc pulseq {} {' + '\n'
            'global par' + '\n'
            'maxdt '+ str(1.0e6 * self.FAMNsteps[stepcnt].data['maxdt']) + '\n')

        # -----------------------------------------------------------
        # Tail
        self.FAMNsteps[stepcnt].data['tail'] = ''
        
        # ---
        gen = fmc.gen_detect_elements(self.output['dim_density'])
        for truc in range(self.output['size_density']):
            self.FAMNsteps[stepcnt].data['tail'] += (
            'matrix set detect elements {' + next(gen) + '}\n' + 'acq\n')
            
        # ---
        self.FAMNsteps[stepcnt].data['tail'] +=\
            (
            '}' + '\n'
            '}' + '\n'
            'proc main {} {' + '\n'
            'global par \n'
            'fsave [fsimpson] ' +
            '$par(name).out' + '\n' + '}' + '\n')
            # './' + self.param['simpson_tempfolder'] + '/'
    
    
    def FAMN_loop(self):
        ''' Runs the main loops of FAM-N optimisation '''
        
        for stepcnt in range(self.output['no_FAMN_steps']):
            
            # Print the current step index if superior to 1
            if self.output['no_FAMN_steps'] > 1:
                print('-----------------------')
                print('Current step : ' + str(stepcnt+1))
                print('')
            
            # Set the maximum number of pulses for each step
            if stepcnt != 0:
                self.temp['noprev_pulses'] = 0
                for truc in range(stepcnt-1):
                    self.temp['noprev_pulses'] += self.FAMNsteps[truc].opt['max_nopules']

            # Break this loop if there are no peak left for further optimisations
            if self.temp['noprev_pulses'] >= self.param['max_nopulses']:
                print('Maximum number of pulses reached for step n°' + str(stepcnt))
                break

            # Initialise class
            self.initialise_step(stepcnt)
            
            # Set the maximum number of pulses to the maximum for the first step …
            if stepcnt == 0:
                self.FAMNsteps[stepcnt].opt['max_nopules'] =\
                self.param['max_nopulses']
            else:
            # … or substract the number of pulses of the previous steps
                self.FAMNsteps[stepcnt].opt['max_nopules'] = \
                self.param['max_nopulses'] - self.temp['noprev_pulses']
                
            # If the list of initial density evaluates of the current step is 'None', 
            # get the final density matrix of previous pulse as initial density matrix of current
            if stepcnt > 0:
                if self.output['list_initial_density'][stepcnt].any() == None :
                    self.output['do_final_sim'][stepcnt] = False
                    self.output['is_succesion'][stepcnt] = True
                    self.output['list_initial_density'][stepcnt] =\
                    self.FAMNsteps[stepcnt-1].data['density_last_inversion']
                else:
                    self.output['do_final_sim'][stepcnt] = True
                
            # Run the method sequence
            self.FAMNsteps[stepcnt].run_singlepulse()
            self.FAMNsteps[stepcnt].optimise_famn()
            
            # If the current step was not propagated from the density matrix at the end of the
            # final pulse: redo the simulation with the previous optimised FAM-N pulse but with
            # the initial density of the previous pulse.
            if ( stepcnt > 0 and self.output['do_final_sim'][stepcnt] ):
                
                print('-----------------------')
                print('Proceeding to the final optimsation for step ' + str(stepcnt+ 1))
                print('')
                self.FAMNsteps[stepcnt].final_simulation()
                
            # Do a few things if the optimisation is successful
            self.FAMNsteps[stepcnt].optimisation_complete()

            # Remove the final pulse of the previous pulse if succession
            if (stepcnt > 0 and self.output['do_final_sim'][stepcnt]):
                numpy.delete(self.FAMNsteps[stepcnt-1].data['pulse_durations_nop'],-1)
                
            # Stop step timer and display if more than one step
            self.FAMNsteps[stepcnt].data['timer'] =\
                fmc.stop_timer(self.FAMNsteps[stepcnt].data['timer'])
            if self.output['no_FAMN_steps'] > 1:
                print('Total step time:', fmc.time_to_string(self.FAMNsteps[stepcnt].data['timer']) )
            
            
    def reopt(self):
        '''Function to redo a FAM-N simulation with experiment values different
        from the ones that it was optimised for. To do if someone asks me and if 
        I have the time.
        '''

    def finalise(self):
        ''' Finalise the FAM-N optimisation and organises the production of the output'''  
        
        # ----------------------------------------------------------
        # Final detect coherence
        self.output['final_coherence'] = self.output['list_detect_density'][  self.output['no_FAMN_steps']-1  ]
        
        # ----------------------------------------------------------
        # Number of pulses
        self.output['total_nopules'] = numpy.zeros(self.output['no_FAMN_steps'], dtype = 'int32')
        for stepcnt in range(self.output['no_FAMN_steps']): 
            self.output['total_nopules'][stepcnt] = int(self.FAMNsteps[stepcnt].data['tnopules'])
            if self.output['is_succesion'][stepcnt]:
                self.output['total_nopules'][stepcnt-1] -= 1
        
        # ----------------------------------------------------------
        # Select the final coherence evolution, pulse duration,
        # and amplitudes of maxima (unsorted)
        self.output['select_coherences'] = []
        self.output['pulse_durations'] = []
        self.output['max_amplitudes'] = []
        
        for stepcnt in range(self.output['no_FAMN_steps']):
            # Swipe all separated FAM-N of the optimisation
            
            temp_select_coherences = []
            temp_pulse_durations = []
            temp_max_amplitudes = []
            
            if self.output['do_final_sim'][stepcnt]:
                # 
                for FAMNpulses in range(self.FAMNsteps[stepcnt].data['tnopules']):
                    temp_select_coherences.append(fmc.get_all_detect_element(
                        self.FAMNsteps[stepcnt].data['final_full_density'][FAMNpulses],
                        self.output['final_coherence']))
                    temp_pulse_durations.append(len(temp_select_coherences[-1]))
                    temp_max_amplitudes.append(temp_select_coherences[-1][\
                            self.FAMNsteps[stepcnt].data['max_pos'][FAMNpulses]])
                    
                # -----------------------------------------------------
                self.output['max_amplitudes'][-1] = self.output['select_coherences'][-1][\
                           self.FAMNsteps[stepcnt].data['final_max_pos']]
            else:
                
                for FAMNpulses in range(self.FAMNsteps[stepcnt].data['tnopules']):
                    if self.param['do_every_cryst']:
                        # Get the powder average if crystallite were calculated separatly
                        temp_select_coherences.append(fmc.get_all_detect_element(
                            numpy.average( self.FAMNsteps[stepcnt].data['selected_full_density'][FAMNpulses], axis = 0  ),
                            self.output['final_coherence'])  )
                    else: # Get the powder average directly if not
                        temp_select_coherences.append(fmc.get_all_detect_element(
                            self.FAMNsteps[stepcnt].data['selected_full_density'][FAMNpulses],
                            self.output['final_coherence']))
                    
                    temp_pulse_durations.append(len(temp_select_coherences[-1]))
                    temp_max_amplitudes.append(temp_select_coherences[-1][\
                            self.FAMNsteps[stepcnt].data['max_pos'][FAMNpulses]])
                    
            # --------------------------------------------------------
            self.output['select_coherences'].append(numpy.array(temp_select_coherences))
            self.output['pulse_durations'].append(numpy.array(temp_pulse_durations))
            self.output['max_amplitudes'].append(numpy.array(temp_max_amplitudes))

        # -----------------------------------------------------
        # Final signal improvements    
        self.output['final_TIR'] = numpy.real(self.output['select_coherences'][-1][-1][ self.FAMNsteps[-1].data['max_pos'][-1] ])/\
                    numpy.real(self.output['select_coherences'][0][0][ self.FAMNsteps[0].data['max_pos'][0] ])
        
        # -----------------------------------------------------
        # Total number of points
        self.output['total_nop'] = numpy.full((self.output['no_FAMN_steps']),int(0))  
        for stepcnt in range(self.output['no_FAMN_steps']):
            self.output['total_nop'][stepcnt] =\
            numpy.sum(self.FAMNsteps[stepcnt].data['pulse_durations_nop'])
        
        # Total duration: self.output['total_nop']*self.FAMNsteps[stepcnt].opt['time_increment']
        
         # -----------------------------------------------------
        # First absolute time of each step
        self.output['initial_time'] = numpy.zeros(self.output['no_FAMN_steps'])
        for stepcnt in range(1,self.output['no_FAMN_steps']): 
            self.output['initial_time'][stepcnt] =\
            1.0e6 * self.FAMNsteps[stepcnt].opt['time_increment'] *\
            self.output['initial_time'][stepcnt-1] +\
            self.FAMNsteps[stepcnt].data['pulse_durations_nop'][-1]
        
        # -----------------------------------------------------
        # Text fractions of the final pulse for output in Topspin and Delta
        self.output['pulse_frac'] = []
        for stepcnt in range(self.output['no_FAMN_steps']): 
            # Pre-allocate
            temp = numpy.full( self.output['total_nopules'][stepcnt] , None )
            
            # Store in a temporary array
            for FAMNpulses in range(self.output['total_nopules'][stepcnt]):
                temp[FAMNpulses] = (str(self.FAMNsteps[stepcnt].data['pulse_durations_nop'][FAMNpulses]) +
                        '/' + str(self.output['total_nop'][stepcnt]))
            # -----------------
            self.output['pulse_frac'].append(temp)

        # -----------------------------------------------------
        self.output['completed'] = True
        if self.output['no_FAMN_steps'] > 1: 
            print('All FAM-N optimisation successful')
        
        # Automatically save current instance
        # self.savestate('__previous_session__')
            
        # Stop and show timer
        self.output['total_execution'] = fmc.stop_timer(self.output['total_execution'])
        print('')
        print('Total execution time:',fmc.time_to_string(self.output['total_execution']))
        print('-----------------------------------------')
    # -------------------------------------------------------------------------
    
    # -------------------------------------------------------------------------
    # ----- Output functions
    
    def screen_report(self):
        ''' Output a summary of the current execution in the console'''
        
        print('')
        print('---- Optimisation report -----')
        print('')
        print('Final Theoretical improvement: ×', numpy.round(self.output['final_TIR'],4))
        print('Pulse durations:', self.output['pulse_durations'])
        print('Total number of points:', self.output['total_nop'])
        # print(self.output['select_coherences'])
        print('')
        
        # ---------------------------------------------------------------------
        print('Positions of maxima:')
        print(str( numpy.real(self.output['max_amplitudes'] )))
        if self.output['no_FAMN_steps'] > 1:
            for stepcnt in range(self.output['no_FAMN_steps']):
                print('Step n°' + str(stepcnt+1)  + ':' +
                str(numpy.round(self.FAMNsteps[stepcnt].data['max_amp'],3)))

        # ---------------------------------------------------------------------
        print('Pulse durations (Number of points):')
        if self.output['no_FAMN_steps'] == 1:
            print(str(self.FAMNsteps[0].data['pulse_durations_nop']))
        else:
            for stepcnt in range(self.output['no_FAMN_steps']):
                print('Step n°' + str(stepcnt+1) + ':' +
                             self.FAMNsteps[stepcnt].data['pulse_durations_nop'])
                    
        # ---------------------------------------------------------------------
        print('Pulse durations (μs):')
        
        if self.output['no_FAMN_steps'] == 1:
            print(str(1.0e6*self.FAMNsteps[0].data['pulse_durations_time']))
        else:
            for stepcnt in range(self.output['no_FAMN_steps']):
                print('Step n°' + str(stepcnt+1) + ':' + 
                      str(numpy.round(1.0e6*self.FAMNsteps[stepcnt].data['pulse_durations_time'],3)) + '\n')


    def plot_allsteps_fig(self):
        ''' Plot a figure containing all steps of the optimisation '''
        
        # Pre-allocate
        self.figures['allsteps_fig'], self.figures['allsteps_ax'] = plt.subplots(1,1, figsize = (10,5))
        self.figures['allsteps_ax'].set_title('Final FAM-N pulse')
        self.figures['allsteps_ax'].set_xlabel('Time /μs')

        # ---------------------------------------------------------------------
        # Plot the figure
        for stepcnt in range(self.output['no_FAMN_steps']):
            # Plot the first pulse
            self.figures['allsteps_ax'].plot(self.output['initial_time'] +\
                    1.0e6 *self.FAMNsteps[stepcnt].opt['time_increment'] *\
                    numpy.arange(self.output['pulse_durations'][stepcnt][0]),
                self.output['select_coherences'][stepcnt][0],
                color = fmc.colour_in_plot(0)
                )
            # Plot subsequent pulses
            for FAMNpulses in range(1,self.FAMNsteps[stepcnt].data['tnopules']):
                # Plot curve
                self.figures['allsteps_ax'].plot(self.output['initial_time'] + 1.0e6 * \
                        self.FAMNsteps[stepcnt].opt['time_increment'] * numpy.arange(\
                       self.FAMNsteps[stepcnt].data['inv_pos_abs'][FAMNpulses-1],
                       self.FAMNsteps[stepcnt].data['inv_pos_abs'][FAMNpulses-1] +\
                       self.output['pulse_durations'][stepcnt][FAMNpulses]),
                       numpy.real(self.output['select_coherences'][stepcnt][FAMNpulses]),
                       color = fmc.colour_in_plot(FAMNpulses)
                       )
                # Trace point of phase inversion
                self.figures['allsteps_ax'].axvline(x = self.output['initial_time'][stepcnt] +\
                                1.0e6 *self.FAMNsteps[stepcnt].opt['time_increment'] *\
                                self.FAMNsteps[stepcnt].data['inv_pos_abs'][FAMNpulses-1]        
                                , color = 'black')

        # ---------------------------------------------------------------------
        # Plot position of first maximum
        self.figures['allsteps_ax'].axvline(1.0e6 *self.FAMNsteps[stepcnt].opt['time_increment'] *\
                                self.FAMNsteps[0].data['max_pos_abs'][0],
                color = 'red', linestyle = '--')
        # Plot position of last maximum
        self.figures['allsteps_ax'].axvline(self.output['initial_time'][-1] +\
                                1.0e6 *self.FAMNsteps[stepcnt].opt['time_increment'] *\
                                self.FAMNsteps[-1].data['max_pos_abs'][-1],
                color = 'red', linestyle = '--')
        
    
    def export_text_report(self, filename, filepath = ''):
        ''' Print a text report for all carried out optimisations'''
        
        print(os.path.normpath(filepath) + os.sep + fmc.change_ext(filename,'txt'))
        
        if filepath == '':
            report = open(fmc.change_ext(filename,'txt') ,'w', encoding='utf-8')
        else:
            report = open(os.path.normpath(filepath) + os.sep + fmc.change_ext(filename,'txt') ,'w', encoding='utf-8')
        
        # ---------------------------------------------------------------------
        report.write('------------------------------------------\n')
        report.write('------- FAM-N optimisation results -------\n')
        # --------------
        report.write('\n------- Input parameters-------\n')
        report.write('Instance name: '+ self.param['simpson_basename'] +'\n')
        # --------------
        report.write('Spin system:\n')
        report.write(self.param['tab'] + 'Nucleus: ' + self.param['nucleus'] + '\n')
        report.write(self.param['tab'] + 'Proton Larmor Frequency: '
                     + str(self.param['proton_larmor_frequency']/1.0e6) + ' MHz\n')
        report.write(self.param['tab'] + 'MAS rate: ' + str(self.param['MAS']) + ' Hz\n')
        report.write(self.param['tab'] + 'RF field: ' + str(self.param['RF']) + ' Hz\n')
        report.write(self.param['tab'] + 'Quadrupolar coupling constant (CQ): ' +
                     str(self.param['CQ']/1.0e6) + ' MHz\n')
        report.write(self.param['tab'] + 'Asymmetry parameters (ηQ): '
                     + str(self.param['ηQ']) + '\n')
        
        # --------------
        report.write('\nTransfer:\n')
        # --------------
        if self.output['simple_opt']:
            # If simple FAM-N optimisation
            report.write(self.param['tab'] + 'Simple optimisation\n ')
            report.write(self.param['tab'] + 'Transfer type: ' + str(self.output['transfer_type']) + '\n')
            report.write(self.param['tab'] + 'Transition: ' + str(self.param['transition']) + '\n')
            # Print density matrix
#            report.write(self.param['tab'] + 'Initial density matrice: \n')   
#            report.write(str(self.output['list_initial_density']))
#            report.write('\n' + self.param['tab'] + 'Detect density matrice: \n')
#            report.write(str(self.output['list_detect_density']) + '\n')
            
        else:
            # If complex transfer
            report.write(self.param['tab'] + 'Complex optimisation\n')
            report.write(self.param['tab'] + 'Number of steps:' + str(self.output['no_FAMN_steps']) + '\n')
            report.write(self.param['tab'] + 'Transfer type: ' + self.output['transfer_type'] + '\n')
            report.write(self.param['tab'] + 'Input optimisation steps: ' +
                         str( self.param['complex_optimisation']) + '\n')
            report.write(self.param['tab'] + 'Input nutation list: ' + 
                         str(self.param['single_pulse_nutation_list']) + '\n')
            # Print density matrices
#            report.write(self.param['tab'] + 'Initial density matrices: \n')   
#            report.write(str(self.output['list_initial_density']))
#            report.write('\n' + self.param['tab'] + 'Detect density matrices: \n')
#            report.write(str(self.output['list_detect_density']) + '\n')

        # --------------
        report.write('\nAdvanced parameters:\n')
        # --------------
        report.write(self.param['tab'] + 'Minimum improvement: ' + str(self.param['min_improvement']) + '\n')
        report.write(self.param['tab'] + 'Maximum number of pulses: ' + str(self.param['max_nopulses']) + '\n')
        report.write(self.param['tab'] + 'Number of points per 360° nutation: ' + str(self.param['nop_360']) + '\n')
        report.write(self.param['tab'] + 'Minimum time increment (μs): '+ str(numpy.round(1.0e6*self.param['min_time_increment'],3)) + '\n')
        report.write(self.param['tab'] + 'maxdt ratio: ' + str(self.param['maxdt_ratio']) + '\n')
        report.write(self.param['tab'] + 'Crystal file: ' + str(self.param['crystfile'])    + '\n')
        report.write(self.param['tab'] + 'Number of γ: ' + str(self.param['no_γ']) + '\n')
        report.write(self.param['tab'] + 'num_cores: ' + str(self.param['num_cores']) + '\n')
        # ---------------------------------------------------------------------
        
        report.write('\n------- Results -------\n')
        
        # ---------------------------------------------------------------------
        report.write('Theoretical improvement ratio: ')
        report.write('×' +str(numpy.round(self.output['final_TIR'],2)))
        if self.output['no_FAMN_steps'] > 1:
            for stepcnt in range(self.output['no_FAMN_steps']):
                report.write(self.param['tab'] + '\nStep n°' + str(stepcnt+1)  + ': ')
                report.write('×' + str(numpy.round(self.FAMNsteps[stepcnt].data['TIR'],2))  + '')
            report.write('\n')
        report.write('\n')
        
        # ---------------------------------------------------------------------
        report.write('Number of pulses: ')
        report.write(str( numpy.sum(self.output['total_nopules']) ))
        if self.output['no_FAMN_steps'] > 1:
            for stepcnt in range(self.output['no_FAMN_steps']):
                report.write(self.param['tab'] + '\nStep n°' + str(stepcnt+1)  + ': ')
                report.write(str(self.FAMNsteps[stepcnt].data['tnopules'])  + '')
            report.write('\n')
        report.write('\n')
        
        # ---------------------------------------------------------------------
        report.write('Time increment (μs): ')
        if self.output['no_FAMN_steps'] == 1:
            report.write(str(numpy.round(1.0e6*self.FAMNsteps[0].opt['time_increment'],3)))
        else:
            for stepcnt in range(self.output['no_FAMN_steps']):
                report.write(self.param['tab'] + '\nStep n°' + str(stepcnt+1)  + ': ')
                report.write(str(numpy.round(1.0e6*self.FAMNsteps[stepcnt].opt['time_increment'],3))  + '')
            report.write('\n')
        report.write('\n')
        
        # --------------
        report.write('Total duration (μs): ')
        if self.output['no_FAMN_steps'] == 1:
            report.write(str(numpy.round(1.0e6*self.output['total_nop'][0]*\
                                         self.FAMNsteps[0].opt['time_increment'] ,3)))
        else:
            for stepcnt in range(self.output['no_FAMN_steps']):
                report.write(self.param['tab'] + '\nStep n°' + str(stepcnt+1)  + ': ')
                report.write(str(numpy.round(1.0e6*self.output['total_nop'][stepcnt]*\
                                    self.FAMNsteps[stepcnt].opt['time_increment'],5))  + '')
            report.write('\n')
        report.write('\n')
        
        # ---------------------------------------------------------------------
        #
        # ---------------------------------------------------------------------
        report.write('''\n– Points of phase inversion (Number of points)\n''')
        
        if self.output['no_FAMN_steps'] == 1:
            report.write(self.param['tab'] + str(self.FAMNsteps[0].data['inv_pos_abs'])  + '\n')
        else:
            for stepcnt in range(self.output['no_FAMN_steps']):
                report.write('Step n°' + str(stepcnt+1)  + ':\n')
                report.write(self.param['tab'] + str(self.FAMNsteps[stepcnt].data['inv_pos_abs'])  + '\n')
        

        # ---------------------------------------------------------------------
        report.write('''\n– Position of maxima (Number of points)\n''')
        
        if self.output['no_FAMN_steps'] == 1:
            report.write(self.param['tab'] + str(self.FAMNsteps[0].data['max_pos_abs'])  + '\n')
        else:
            for stepcnt in range(self.output['no_FAMN_steps']):
                report.write('Step n°' + str(stepcnt+1)  + ':\n')
                report.write(self.param['tab'] + str(self.FAMNsteps[stepcnt].data['max_pos_abs'])  + '\n')
                    
                # self.data['final_max_pos']
                
        # ---------------------------------------------------------------------
        report.write('''\n– Amplitude of maxima\n''')
        
        report.write(self.param['tab'] + str( numpy.real(self.output['max_amplitudes'] )))
        if self.output['no_FAMN_steps'] > 1:
            for stepcnt in range(self.output['no_FAMN_steps']):
                report.write('Step n°' + str(stepcnt+1)  + ':\n')
                report.write(self.param['tab'] + str(numpy.round(self.FAMNsteps[stepcnt].data['max_amp'],3)) + '\n')
            report.write('\n')
        report.write('\n')
                    
        # ---------------------------------------------------------------------
        report.write('''\n– Pulse durations (Number of points)\n''')  # self.data['pulse_durations_nop']
        
        if self.output['no_FAMN_steps'] == 1:
            report.write(self.param['tab'] + str(self.FAMNsteps[0].data['pulse_durations_nop']) + '\n')
        else:
            for stepcnt in range(self.output['no_FAMN_steps']):
                report.write('Step n°' + str(stepcnt+1)  + ':\n')
                report.write(self.param['tab'] +
                             str(self.param['tab'] +
                            self.FAMNsteps[stepcnt].data['pulse_durations_nop']) + '\n')
                    
        # ---------------------------------------------------------------------
        report.write('''\n– Pulse durations (μs)\n''')
        
        if self.output['no_FAMN_steps'] == 1:
            report.write(self.param['tab'] + str(1.0e6*self.FAMNsteps[0].data['pulse_durations_time']) + '\n')
        else:
            for stepcnt in range(self.output['no_FAMN_steps']):
                report.write('Step n°' + str(stepcnt+1) + ':\n')
                report.write(self.param['tab'] +
                             str(self.param['tab'] +\
                        numpy.round(1.0e6*self.FAMNsteps[stepcnt].data['pulse_durations_time'],3)) + '\n')
        
        # --------------
        report.close()
        # --------------

    # -------------------------------------------------------------------------
    # ----- Output pulse programs
    # ---------------------------------------
    #  Bruker pulse programs common
    
    def gen_pulseidx(self, first):
        ''' Generate pulse indexes starting from "first" '''
        for truc in range(numpy.sum(self.output['total_nopules'])):
            yield str(int(first + truc + 1))

    
    def get_final_pp_comments(self,prg):
        ''' Return the comment written at the bottom of bruker pulse sequences.
        prg = 'jeol' or 'bruker'. It switches the comment symbol.
        '''
        
        # Set comment
        if prg == 'bruker':
            com = ';'
        elif 'jeol':
            com = '-- '
        
        # Set the comment for the transfer
        if self.output['simple_opt']:
            # If simple FAM-N optimisation
            transfer = (com + 'Simple optimisation\n' + com + ' Transfer type: ' +
                        str(self.output['transfer_type']) + 
            '\n ' + com + 'Transition: ' + str(self.param['transition']) + '\n')
            
        else:
            # If complex transfer
            transfer = (
            com + 'Complex optimisation\n' +
            com + 'Number of steps:' + str(self.output['no_FAMN_steps']) + '\n' +
            com + 'Transfer type: ' + self.output['transfer_type'] + '\n' +
            com + 'Input optimisation steps: ' +
                         str( self.param['complex_optimisation']) + '\n' +
            com + 'Input nutation list: ' + 
                         str(self.param['single_pulse_nutation_list']) + '\n')

        # Add to the rest
        return (com + '-------- FAM-N Opt Details --------\n' + 
        com + 'Instance name: ' + self.param['simpson_basename'] + '\n' +
        com + 'Spin system: \n' +
        com + 'Nucleus: ' + self.param['nucleus'] +
        '\n' + com + 'Proton Larmor Frequency: '
                     + str(self.param['proton_larmor_frequency']/1.0e6) +
        ' MHz\n' + com + 'MAS rate: ' + str(self.param['MAS']) +
        ' Hz\n' + com + 'RF field: ' + str(self.param['RF']) +
        ' Hz\n' + com + 'Quadrupolar coupling constant (CQ): ' +
                     str(self.param['CQ']/1.0e6) +
        ' MHz\n' + com + 'Asymmetry parameters (ηQ): '
                     + str(self.param['ηQ']) +
        '\n\n' + com + 'Transfer:\n' + transfer)
    
    # ---------------------------------------------------------------------
    # ----- Export Bruker
    
    def bruker_get_variable_parts_pp(self, first):
        ''' Get variable components of Bruker pulse programs '''
        
        # Pulse list
        pulse_list = ''
        pidx = self.gen_pulseidx(first)
        for stepcnt in range(self.output['no_FAMN_steps']):
            for FAMNpulses in range(self.FAMNsteps[stepcnt].data['tnopules']):
                pulse_list += ('"p' +  pidx.__next__()  + '=(1u*' +\
                    self.output['pulse_frac'][stepcnt][FAMNpulses]  +\
                    ")*cnst" + str(stepcnt+1) + '"\n')
                
        # Insert conversion pulse with phase inversions
        conversion_pulses = ''
        pidx = self.gen_pulseidx(first)
        positive_phase = True
        for stepcnt in range(self.output['no_FAMN_steps']):
            for FAMNpulses in range(self.FAMNsteps[stepcnt].data['tnopules']):
                conversion_pulses += ('(p' + pidx.__next__() + ' ph' +\
                                fmc.get_bruker_phase_alternation(positive_phase) + '):f1\n')
                positive_phase = not positive_phase
                
                
        # Total duration
        total_duration = ''
        for stepcnt in range(self.output['no_FAMN_steps']):
            total_duration += (';cnst' + str(stepcnt+1) +\
                                ' : FAM-N duration (should be ' +\
                                str(numpy.round(1.0e6*self.output['total_nop'][stepcnt]*\
                                    self.FAMNsteps[stepcnt].opt['time_increment'],5))
                                + ' us) \n')
            
        return pulse_list, conversion_pulses, total_duration
    
    def export_topspin_pulselist(self,filename, filepath = '', first = 1):
        ''' Export a pulse list for Bruker input'''
        
        if filepath == '':
            topspin_pulselist = open(fmc.change_ext(filename,'txt') ,'w', encoding='utf-8')
        else:
            topspin_pulselist = open(os.path.normpath(filepath) + os.sep + fmc.change_ext(filename,'txt') ,'w', encoding='utf-8')
        # ---------------------------------------------------------------------
        # Get the elements
        pulse_list, conversion_pulses, total_duration = self.bruker_get_variable_parts_pp(first)
        
        # print the list
        topspin_pulselist.write('------- Pulse list -------\n\n\n')
        topspin_pulselist.write(pulse_list)
        topspin_pulselist.write('\n')
        topspin_pulselist.write(conversion_pulses)
        topspin_pulselist.write('\n')
        topspin_pulselist.write(total_duration)
    
        # ---------------------------------------------------------------------
        print('Pulse list for Topspin successfully exported')
        topspin_pulselist.close()
    
    def export_topspin_MQF(self, filename, filepath = ''):
        ''' Writes a multiple-quantum filtered NMR pulse program for Bruker Topspin'''
        
        if filepath == '':
            topspin_MQF = open(fmc.change_ext(filename,'') ,'w', encoding='utf-8')
        else:
            topspin_MQF = open(os.path.normpath(filepath) + os.sep + fmc.change_ext(filename,'') ,'w', encoding='utf-8')
        # ---------------------------------------------------------------------
        # Get some elements
        phase_cycling = fmc.get_bruker_mqf_phasecycling(self.param['transition'])
        
        pulse_list, conversion_pulses, total_duration = self.bruker_get_variable_parts_pp(1)
        
        # Introduction
        topspin_MQF.write(';Multiple-quantum filtered MAS experiment\n' +
                            ';transition: ' + str(self.param['transition']) + 'Q\n\n'
                            )
        
        # Pulse list
        topspin_MQF.write(pulse_list)
        
        # Beginning of pulse program, split-t1 constants and excitation pulse
        topspin_MQF.write('\n1 ze \n' +
                            '2 1u pl1:f1 \n' +
                            '  d1 \n' +
                            '  (p1 ph1):f1 \n' +
                            '  d2 pl2:f1 \n' +
                            '3 ')
        
        # Insert conversion pulse with phase inverson
        topspin_MQF.write(conversion_pulses)        
        
        # End of the pulse program
        topspin_MQF.write('go=2 ph31 \n' +
                            '  wr #0 \n' +
                            '  exit \n\n')
        
        # Phase cycling
        topspin_MQF.write(phase_cycling + '\n\n')
        
        # Total duration
        topspin_MQF.write(total_duration)
        
        # Other comments
        topspin_MQF.write(';pl1 : Excitation Power Level \n' +
            ';p1 : Excitation Pulse Length \n' +
            ';pl2 : Conversion Power Level \n' +
            ';d1 : Recycle Delay \n' +
            ';d0 : Initial Interpulse Delay \n\n')
        
        # Optimisation comments
        topspin_MQF.write(self.get_final_pp_comments('bruker'))
        # ---------------------------------------------------------------------
        print('MQF pulse program for Topspin successfully exported')
        topspin_MQF.close()
        
        
        
    def export_topspin_MQMAS(self, filename, filepath = ''):
        ''' Writes a Multiple-Quantum Magic-Angle-Spinning NMR pulse program
        for Bruker Topspin'''
        
        if filepath == '':
            topspin_MQMAS = open(fmc.change_ext(filename,'') ,'w', encoding='utf-8')
        else:
            topspin_MQMAS = open(os.path.normpath(filepath) + os.sep + fmc.change_ext(filename,'') ,'w', encoding='utf-8')
        # ---------------------------------------------------------------------
        # Get some elements
        splitt1 = fmc.get_splitt1(self.param['spin'],self.param['transition'])
        phase_cycling = fmc.get_bruker_mqmas_phasecycling(self.param['transition'])
        pulse_list, conversion_pulses, total_duration = self.bruker_get_variable_parts_pp(2)
        
        # Introduction
        topspin_MQMAS.write(';MQMAS split-t1 with FAM-N conversion\n' +
                            ';for spin ' + self.param['spin'] + ' nuclei\n' +
                            ';transition: ' + str(self.param['transition']) + 'Q\n\n' +
                          'define delay sqd \n' + 
                          '"inf1=in0*' + splitt1[4]
                          + '"\n\n'
                          )
        
        # Pulse list
        topspin_MQMAS.write(pulse_list)
        
        # Beginning of pulse program, split-t1 constants and excitation pulse
        topspin_MQMAS.write('\n  ze\n' +
                            '  "sqd=d0*' + splitt1[2] + '"\n'
                            '1 d1\n' +
                            '3 10u pl1:f1\n' +
                            '(p1 ph1):f1\n' +
                            'd0 pl2:f1\n' +
                            '4 ')
        
        # Insert conversion pulse with phase inverson
        topspin_MQMAS.write(conversion_pulses)        
        
        # Delay
        topspin_MQMAS.write(' d6 pl13:f1 \n')
        
        # Change position of the single-quantum delay if required
        if fmc.get_splitt1_pos(self.param['spin'],self.param['transition']):
            topspin_MQMAS.write(' sqd \n (p2 ph10):f1\n')
        else:
            topspin_MQMAS.write(' (p2 ph10):f1\n sqd \n')
                            
        # End of the pulse program
        topspin_MQMAS.write(' go=1 ph31 \n' +
                            ' d1 mc #0 to 1 F1QF(id0) \n' +
                            ' exit \n\n')
        
        # Phase cycling
        topspin_MQMAS.write(phase_cycling + '\n\n')
        
        # Total duration
        topspin_MQMAS.write(total_duration)
        
        # Other comments
        topspin_MQMAS.write(
        ';pl1 : Excitation Power Level \n' +
        ';p1 : Excitation Pulse Length \n' +
        ';pl2 : Conversion Power Level \n' +
        ';p2 : 180 CT Selective Pulse Length \n' +
        ';d0 : Initial Interpulse Delay \n' +
        ';d1 : Recycle Delay \n' +
        ';d6 : Echo Delay \n\n')
        
        # Optimisation comments
        topspin_MQMAS.write(self.get_final_pp_comments('bruker'))
        # ---------------------------------------------------
        print('MQMAS pulse program for Topspin successfully exported')
        topspin_MQMAS.close()
        
    
    def export_topspin_MQCPMGMAS(self, filename, filepath = ''):
        ''' Writes a multiple-quantum Carr-Purcell-Meiboom-Gill Magic-Angle-Spinning
        NMR pulse program for Bruker Topspin'''
    
        if filepath == '':
            topspin_MQCPMGMAS = open(fmc.change_ext(filename,'') ,'w', encoding='utf-8')
        else:
            topspin_MQCPMGMAS = open(os.path.normpath(filepath) + os.sep + fmc.change_ext(filename,'') ,'w', encoding='utf-8')
        # ---------------------------------------------------------------------
        # Get some elements
        phase_cycling = fmc.get_bruker_mqcpmgmas_phasecycling(self.param['transition'])
        splitt1 = fmc.get_splitt1(self.param['spin'],self.param['transition'])
        pulse_list, conversion_pulses, total_duration = self.bruker_get_variable_parts_pp(3)
        
        # Introduction
        topspin_MQCPMGMAS.write(';CPMG detected MQMAS split-t1 with FAM-N conversion\n' +
                            ';for spin ' + self.param['spin'] + ' nuclei\n' +
                            ';transition: ' + str(self.param['transition']) + 'Q\n\n' +
                            
                            '; CPMG part written by Stefan Steuernagel\n' +
                            '; pulse program for MQCPMGMAS sequence\n' +
                            '; samples continuously, including ALL pulses and ringdown delays\n' +
                            '; may be used with digmod digital\n' +
                            '; important: only runs via SGU in channel 3\n\n' +
        
                            ';$CLASS=Solids\n' +
                            ';$DIM=2D\n' +
                            ';$TYPE=half integer quadrupoles\n' +
                            ';$SUBTYPE=simple 2D\n' +
                            ';$OWNER=Bruker\n' +
                            '#include <Avancesolids.incl>\n' +
                            ';#include <Delayssolids.incl> ;This worked for topspin 2.1, but not 3.1.6\n\n' +
        
                            'define delay del3\n' +
                            '"del3=d3-1u"\n' +
                            '"cnst30=((d6*2+d3*2+p3)*l22+d6+d3)/dw"\n' +
                            'define delay rest\n' +
                            '"rest=aq-(cnst30*dw)"\n' +
                            '"d6=((l12*1s/cnst31)-d3*2-p3)/2"\n' +
                            '"d2=(((l13*1s/cnst31)-p3)/2)"\n' +
                            '"d22=d2-1u"\n' +
                            ';cnst11 : to adjust t=0 for acquisition, if digmod = baseopt\n' +
                            '"acqt0=1u*cnst11"\n' +
                            'define delay echospacing\n' +
                            '"echospacing=(2*d6)+(2*del3)+2u+p3"\n' +
                            'define delay peakspacing\n' +
                            '"peakspacing = 1/echospacing"\n\n' +
                            
                            'define delay sqd \n' + 
                            '"inf1=in0*"' + splitt1[4]
                            + '\n\n'
                            )
        
        # Pulse list
        topspin_MQCPMGMAS.write(pulse_list)
        
        # Beginning of pulse program, split-t1 constants and excitation pulse
        topspin_MQCPMGMAS.write('\n1 ze \n' +
                            '2 1u pl1:f1 \n' +
                            '   "sqd=d0*' + splitt1[2] + '"\n'
                            '  d1 \n' +
                            '  (p1 ph1):f1 \n' +
                            '  d2 pl2:f1 \n' +
                            '3 ')
        
        # Insert conversion pulse with phase inverson
        topspin_MQCPMGMAS.write(conversion_pulses)        
        topspin_MQCPMGMAS.write('  10u pl22:f1\n')
                                
        # Position of the single-quantum delay depending on the spin and transition
        if fmc.get_splitt1_pos(self.param['spin'],self.param['transition']):
            topspin_MQCPMGMAS.write('  sqd\n (p3  ph2):f1\n')
        else:
            topspin_MQCPMGMAS.write(' (p2 ph10):f1\n sqd\n') 
                                
        # CPMG part and end of the pulse program
        topspin_MQCPMGMAS.write('4 d2\n' +
                                '  (p3  ph3)\n' +
                                '  d2\n' +
                                '  lo to 4 times l21\n' +
                                '  d2\n' +
                                '  (p3  ph3)\n' +
                                '  d22\n' +
                                '  1u DWL_CLK_ON\n' +
                                '5 d6 RG_ON  ;;;;;;\n' +
                                '  1u RG_OFF \n' +
                                '  del3 \n' +
                                '  (p3  ph4):f1\n' +
                                '  1u \n' +
                                '  del3 \n' +
                                '  d6 RG_ON \n' +
                                '  lo to 4 times l22 ;;;;;\n' +
                                '  d6\n' +
                                '  del3 RG_OFF \n' +
                                '  rest\n' +
                                '  1u DWL_CLK_OFF\n' +
                                '  rcyc=2\n' +
                                '  10m mc #0 to 2 F1PH(ip1,id0)\n' +
                                'exit \n\n')
        
        # Phase cycling
        topspin_MQCPMGMAS.write(phase_cycling + '\n\n')
        
        # Total duration
        topspin_MQCPMGMAS.write(total_duration)
        
        # Other comments
        topspin_MQCPMGMAS.write('\n' +
                                ';cnst30 : set td to number of acquired complex data points\n' +
                                ';cnst31 : MAS frequency in Hz\n' +
                                ';ns : 96 * n\n' +
                                ';d0 : Initial 3QF delay\n' +
                                ';d1 : recycle delay\n' +
                                ';d2 : delay after first 180 (no acq)\n' +
                                ';d3 : time to allow pulse ringdown, 10 to 100 us\n' +
                                ';d6 : half-echo duration\n' +
                                ';pl1 : =120 dB, not used\n' +
                                ';pl20: MQ Excitation power level\n' +
                                ';pl21 : MQ Conversion RF power level\n' +
                                ';pl22 : 180 RF power level\n' +
                                ';p1 : MQ Excitation pulse\n' +
                                ';p3 : CT selective 180 pulse\n' +
                                ';p4 : MQ Conversion pulse\n' +
                                ';l12 : number of rotor periods for acquisition echoes\n' +
                                ';l13 : number of rotor periods for first echo\n' +
                                ';l21 : number of non-acquired CPMG echoes\n' +
                                ';l22 : number of acquired CPMG echoes\n')
        
        # Optimisation comments
        topspin_MQCPMGMAS.write(self.get_final_pp_comments('bruker'))
        # ---------------------------------------------------------------------
        print('MQ-CPMG-MAS pulse program for Topspin successfully exported')
        topspin_MQCPMGMAS.close()
        
    # -------------------------------------
    #
    #
    # -------------------------------------
    # Jeol pulse programs
    
    # ----- Export delta
    def jeol_get_variable_parts_pp(self):
        ''' Get variable components of Jeol pulse programs '''
        
        # Pulse list
        pulse_list = ''
        for stepcnt in range(self.output['no_FAMN_steps']):
            for FAMNpulses in range(self.FAMNsteps[stepcnt].data['tnopules']):
                pulse_list += ('    conv_mq1q_s' + str(stepcnt+1) + '_pw' + str(FAMNpulses+1) +
                               ' =? famn_ttd' + str(stepcnt + 1) +
                               ' * ' + self.output['pulse_frac'][stepcnt][FAMNpulses]
                               + ', help "FAMN pulse";\n')
                
        # Insert conversion pulse with phase inversions
        conversion_pulses = ''
        positive_phase = True
        for stepcnt in range(self.output['no_FAMN_steps']):
            for FAMNpulses in range(self.FAMNsteps[stepcnt].data['tnopules']):
                conversion_pulses += ('    conv_mq1q_s' + str(stepcnt+1) + '_pw' +
                                str(FAMNpulses+1) + ', (obs.gate, obs.phs.obs_phs_MQconv'
                        + fmc.get_jeol_phase_alternation(positive_phase) + ', obs.amp.x_amp_convpulse, obs.atn.x_atn);\n')
                positive_phase = not positive_phase
                
        # Total duration
        total_duration = ''
        total_duration_comment = ''
        for stepcnt in range(self.output['no_FAMN_steps']):
            duration_str = str(numpy.round(1.0e6*self.output['total_nop'][stepcnt]*\
                        self.FAMNsteps[stepcnt].opt['time_increment'],5))
            total_duration += ('	famn_ttd' + str(stepcnt+1) +' 			=> ' + duration_str
                         + ' [us], 0[us]->1[ms]:10[ns],' +
                        ' help "Total duration of the FAM-N pulse, should be ' +
                        duration_str + ' us";\n')
            total_duration_comment += ('--  famn_ttd' + str(stepcnt+1) +
                                       ' : total duration of FAM-N (should be ' + duration_str + ' us) \n')
            
        return pulse_list, conversion_pulses, total_duration, total_duration_comment
    
    
    def export_delta_pulselist(self,filename, filepath = ''):
        ''' Export a pulse list for Bruker input'''
    
        if filepath == '':
            delta_pulselist = open(fmc.change_ext(filename,'txt') ,'w', encoding='utf-8')
        else:
            delta_pulselist = open(os.path.normpath(filepath) + os.sep + fmc.change_ext(filename,'txt') ,'w', encoding='utf-8')
        # ---------------------------------------------------------------------
        # Get the elements
        pulse_list, conversion_pulses, total_duration, _ = self.jeol_get_variable_parts_pp()
        
        # print the list
        delta_pulselist.write('------- Pulse list -------\n\n\n')
        delta_pulselist.write(pulse_list)
        delta_pulselist.write('\n')
        delta_pulselist.write(conversion_pulses)
        delta_pulselist.write('\n')
        delta_pulselist.write(total_duration)
    
        # ---------------------------------------------------------------------
        print('Pulse list for Delta successfully exported')
        delta_pulselist.close()
    
    
    def export_delta_MQF(self, filename, filepath = ''):
        ''' Function that writes the report of the optimisation for Jeol Delta'''
        
        if filepath == '':
            delta_MQF = open(fmc.change_ext(filename,'jxp') ,'w', encoding='utf-8')
        else:
            delta_MQF = open(os.path.normpath(filepath) + os.sep + fmc.change_ext(filename,'jxp') ,'w', encoding='utf-8')
        # ---------------------------------------------------------------------
        phase_cycling, required_scan = fmc.get_jeol_mqf_phasecycling(self.param['transition'])
        pulse_list, conversion_pulses, total_duration, total_duration_comment = self.jeol_get_variable_parts_pp()
        
        # ---------------------------------------------------
        # Get some elements
        delta_MQF.write('-- HELP.eng: ' + str(self.param['transition']) + 'QF-MAS pulse programme\n'
        '-- Category: mqmas, solids\n'
        '-- File name : ' + str(self.param['transition']) + 'Q-filtered MAS experiment, for optimisation purposes\n'
        '--\n'
        '-- Sequence name :\n'
        '--\n'
        '-- Parameter :\n'
        '--  excite_mq_pw    : ' + str(self.param['transition']) + 'Q excite pulse width\n'
        '--  x_amp_pulse     : amplifier for the excitation pulse\n'
        '--  x_atn           : attenuator of x_pulse\n'
        '--\n' +
        '\n')
        
        # Write pulse lengths in the comments
        delta_MQF.write(total_duration_comment)
        
        # Write end of the introduction
        delta_MQF.write('\n'
        '--  conv_mq1q_pw : conversion pulse\n'
        '--  x_amp_convpulse     : ampliper of conversion pulse\n'
        '--\n'
        '--  relaxation_delay    : inter-pulse delay\n'
        '--  repetition_time     : pulse repetition_time (= relaxation_delay+x_acq_time)\n'
        '--\n'
        '-- Note :\n'
        '--  scans =' + str(required_scan) + '*n\n'
        '--\n'
        '-- END HELP\n\n')
         
        # Set file metadata
        delta_MQF.write('header\n'
        '	filename			=> "";\n'
        '	sample_id			=> "";\n'
        '	comment				=> "";\n'
        '	include "header_solid";\n'
        '	round_points		=> false;\n'
        '    process     		=> "1d_solid.list";\n'
        'end header;\n'
        '\n'
        'instrument\n'
        '	include "instrument_solid";\n'
        '	recvr_phase			=> 0[deg], help "Adjust to put all signal in Real Channel";\n'
        'end instrument;\n\n')
        
        # Set acquisition
        delta_MQF.write('acquisition \n'
        '	x_domain			=> "' +
        fmc.read_deltaname_from_xml(fmc.chemical_symbol_first(self.param['nucleus'])) + '" ;\n'
        '	x_offset			=> 0[ppm];\n'
        '	x_sweep				=> 15[kHz];\n'
        '	x_points			=> 512;\n'
        '	req_scans			=>' + str(required_scan) + ', help "Number of scans";\n'
        '	scans				=?' + str(required_scan) + '* Upper (req_scans / 12), help "Upper scans to nearest multiple of phase cycle";\n'
        '	x_prescans			=> 0;\n'
        '	mod_return			=> 1;\n'
        '\n'
        '	include "acquisition_solid";\n'
        '\n'
        'end acquisition; \n\n')
        
        # Print the first part of the pulse definition
        delta_MQF.write('pulse\n'
        '	collect COMPLEX,OBS;\n'
        '	include "pulse_solid";\n'
        '	initial_wait			=  10.0[ms];\n'
        '	excite				=?"#Setup Excitation Pulse#";\n'
        '	excite_mq_pw			=> 3[us], 0[us]->1[ms]:10[ns], help "3Q excite pulse width";\n'
        '	x_amp_pulse			=> 100[%],  0[%]->100[%]:0.01[%], help "set amplifier for 3Q excitation pulse";\n'
        '	conv				=?"#FAM-N pulse Pulse#";\n')
        
        # Print pulse duration and pulse list
        delta_MQF.write(total_duration)
        delta_MQF.write('	x_amp_convpulse		=> 100[%],  0[%]->100[%]:0.01[%], help "set amplifier for 3Q conversion pulse";\n')
        delta_MQF.write(pulse_list)
        
        # Write phase cycling
        delta_MQF.write('	mqevolution			=?"#Multiple-quantum evolution delay#";\n'
        '	MQ_delay			=> 0.0[us], help "3Q evolution delay";\n'
        '	recycle_Setup			=? "#Setup Recycle Times#";\n'
        '	relaxation_delay		=> 5.0[s], help "inter-pulse delay";\n'
        '	repetition_time			=? relaxation_delay + x_acq_time, help "relaxation_delay+x_acq_time";\n'
        '	include "dec_solid_optional";\n'
        '	atn_Setup         		=? "#Experiment Attenuator Settings#";\n'
        '	x_atn				=> xatn, help "set attenuator for x pulse";\n'
        '	When Decoupling do\n'
        '		irr_atn			=>	irratn, help "attenuator for irr";\n'
        '	end when;\n'
        )
        
        # Write end of the pulse program
        delta_MQF.write(phase_cycling)
        
        # Pulse list
        delta_MQF.write('	obs_phs_MQconv0				=  {0};\n'
        '	obs_phs_MQconv180			=  {180};\n'
        '	obs_phs_phase_acq			=  {0,90,180,270};\n\n'
        '	module_config = "solid_sample"; \n\n')
        
        # Position of the single-quantum delay depending on the spin and the transition
        delta_MQF.write('begin\n'
        '	initial_wait;\n'
        '	relaxation_delay;\n'
        '	when Decoupling do\n'
        '		on        (irr.gate, irr.phs.0, irr.amp.irr_amp_dec, irr.atn.irr_atn, irr.noise.irr_noise);\n'
        '	end when;\n'
        '	excite_mq_pw, (obs.gate, obs.phs.obs_phs_MQexc, obs.amp.x_amp_pulse, obs.atn.x_atn);\n'
        '	MQ_delay;\n')

        # Pulse list
        delta_MQF.write(conversion_pulses)

        # End of the pulse program
        delta_MQF.write('	acq (dead_time, delay, obs_phs_phase_acq);\n'
        '	when Decoupling do\n'
        '		off (irr.gate);\n'
        '	end when;\n'
        'end pulse;\n\n')
        
        # Write optimisation comments
        delta_MQF.write(self.get_final_pp_comments('jeol'))
        
        # ---------------------------------------------------------------------
        print('MQF pulse program for Delta successfully exported')
        delta_MQF.close()
        
        
        
    def export_delta_MQMAS(self, filename, filepath = ''):
        ''' Writes a Multiple-Quantum Magic-Angle-Spinning NMR pulse program
        for Jeol Delta'''
        
        if filepath == '':
            delta_MQMAS = open(fmc.change_ext(filename,'jxp') ,'w', encoding='utf-8')
        else:
            delta_MQMAS = open(os.path.normpath(filepath) + os.sep + fmc.change_ext(filename,'jxp') ,'w', encoding='utf-8')
        # --------------------------------------------------------------------- 
        # Get some elements
        phase_cycling, required_scan = fmc.get_jeol_mqmas_phasecycling(self.param['transition'])
        splitt1 = fmc.get_splitt1(self.param['spin'],self.param['transition'])
        pulse_list, conversion_pulses, total_duration, total_duration_comment = self.jeol_get_variable_parts_pp()
        
        # Write header
        delta_MQMAS.write('-- HELP.eng: ' + str(self.param['transition']) + 'QMAS pulse programme\n'
        '-- Category: mqmas, solids\n'
        '-- File name : shifted-echo ' + str(self.param['transition']) +
        'QMAS experiment with split-t1 for spins ' + str(self.param['spin']) + '\n'
        '-- Sequence name : Please adjust pulse widths for ' + str(self.param['transition']) + 'Q excitation and echo-delay\n'
        '-- \n'
        '-- Parameter:\n'
        '--  excite_mq_pw    : 3Q excite pulse width\n'
        '--  x_amp_pulse     : ampliper of excite pulse\n'
        '--  x_atn           : attenuator of x_pulse\n'
        '-- \n')
        
        # Write pulse lengths in the comments
        delta_MQMAS.write(total_duration_comment)
        
        # Write end of the introduction
        delta_MQMAS.write('--  conv_mq1q_pw : conversion pulses \n'
        '--  x_amp_convpulse     : ampliper of conversion pulse \n'
        '-- \n'
        '--  relaxation_delay    : inter-pulse delay \n'
        '--  repetition_time     : pulse repetition_time (= relaxation_delay+x_acq_time) \n'
        '-- \n'
        '-- tau             : Echo delay, normally 1/2 acq time \n'
        '-- \n'
        '-- Note : \n'
        '--  scans = ' + str(required_scan) + '*n \n'
        '-- \n'
        '-- END HELP \n\n')
        
        # Set file metadata
        delta_MQMAS.write('header \n' +
        '	filename			=> "";\n' +
        '	sample_id			=> "";\n' +
        '	comment				=> "";\n' +
        '	include "header_solid";\n' +
        '	round_points		=> false;\n' +
        '    process     		=> "1d_solid.list";\n' +
        'end header;\n' +
        '\n' +
        'instrument\n' +
        '	include "instrument_solid";\n' +
        '	recvr_phase			=> 0[deg], help "Adjust to put all signal in Real Channel";\n' +
        'end instrument;\n\n')
        
        
        # Set acquisition
        delta_MQMAS.write('acquisition \n' +
        '	x_domain			=> "' + fmc.read_deltaname_from_xml(fmc.chemical_symbol_first(self.param['nucleus'])) + '"; \n' +
        '	x_offset			=> 0[ppm]; \n' +
        '	x_sweep				=> 15[kHz]; \n' +
        '	x_points			=> 512; \n' +
        '	req_scans			=> ' + str(required_scan) +', help "Number of scans"; \n' +
        '	scans				=? ' + str(required_scan) + ' * Upper (req_scans /' + str(required_scan) +
                '), help "Upper scans to nearest multiple of phase cycle";' +
        '	x_prescans			=> 0;\n' +
        '	mod_return			=> 1;\n' +
        '\n' +
        '	y_domain			=  x_domain;\n' +
        '	y_offset			=  x_offset;\n' +
        '	mqd_rotorperiod		=> 50[us], help "MQ evolution delay, multiple of divisor of the MAS rotor period";\n' +
        '   sqd		=? ' + splitt1[2] + ' * mqd_rotorperiod, help "SQ evolution delay";\n' +
        '	y_sweep				=? mqd_rotorperiod + sqd ;\n' +
        '	y_dwell				=? 1/y_sweep;\n' +
        '	y_points			=> 64;\n' +
        ' \n' +
        '	include "acquisition_solid";\n' +
        '	include "acquisition_2d" \n' +
        'end acquisition; \n\n')
        
        # Print the first part of the pulse definition
        delta_MQMAS.write('pulse\n' +
        '	collect COMPLEX,OBS REAL,OBS,NTYPE;\n' +
        '	include "pulse_solid";\n' +
        '	initial_wait			=  10.0[ms];\n' +
        '\n' +
        '	eexcite				=?"#Setup Excitation Pulse#";\n' +
        '	excite_mq_pw			=> 3[us], 0[us]->1[ms]:10[ns], help "3Q excite pulse width";\n' +
        '	x_amp_pulse			=> 100[%],  0[%]->100[%]:0.01[%], help "set amplifier for 3Q excitation pulse";\n' +
        '\n' +
        '	conv				=?"#FAM-N pulse Pulse#";\n')
        
        # Print pulse duration and pulse list
        delta_MQMAS.write(total_duration)
        delta_MQMAS.write(' 	x_amp_convpulse		=> 100[%],  0[%]->100[%]:0.01[%],' +
                          ' help "set amplifier for 3Q conversion pulse";\n')
        delta_MQMAS.write(conversion_pulses)
        delta_MQMAS.write(' 	CTsel				=?"#CT-selective 180degree pulse#";\n'  +
        '	sel_180_pulse			=> 10 [us], help "Duration for the CT-selective 180degree pulse";\n'  +
        '	sel_amp_pulse			=> 10[%],  help "set amplifier for CT-selective 180degree pulse";\n'  +
        '\n'  +
        '	mqevolution			=?"#Delays#";\n'  +
        '	initial_t1			=> 0.0[us], help "Initial 3Q evolution delay";\n'  +
        '	tau				=> x_acq_time / 2, help "Echo delay";\n'  +
        '	Phase_1ocorr =? 360[deg] * (7/9) * tau / x_dwell, help "full-echo first-order phase correction";\n'  + 
        '\n'  +
        '	recycle_Setup			=? "#Setup Recycle Times#";\n'  +
        '	relaxation_delay		=> 5.0[s], help "inter-pulse delay";\n'  +
        '	repetition_time			=? relaxation_delay + x_acq_time, help "relaxation_delay+x_acq_time";\n'  +
        '\n'  +
        '	include "dec_solid_optional";\n'  +
        '\n'  +
        '	atn_Setup           		=? "#Experiment Attenuator Settings#";\n'  +
        '	x_atn				=> xatn, help "set attenuator for x pulse";\n'  +
        '	When Decoupling do\n'  +
        '		irr_atn			=>	irratn,										help "attenuator for irr";\n'  +
        '	end when;\n'  +
        '\n')
        
        # Write phase cycling
        delta_MQMAS.write('	obs_phs_MQconv0		=  {0};\n' +
        '	obs_phs_MQconv180	=  {180};\n')
        delta_MQMAS.write(phase_cycling)
        delta_MQMAS.write('module_config = "solid_sample";')
        
        # Write end of the pulse program
        delta_MQMAS.write('begin\n' +
        '	initial_wait;\n' +
        '	relaxation_delay;\n' +
        '	when Decoupling do\n' +
        '		on        (irr.gate, irr.phs.0, irr.amp.irr_amp_dec, irr.atn.irr_atn, irr.noise.irr_noise);\n' +
        '	end when;\n' +
        '	excite_mq_pw, (obs.gate, obs.phs.obs_phs_MQexc, obs.amp.x_amp_pulse, obs.atn.x_atn);\n' +
        '	(' + splitt1[1] + ' * initial_t1) ystep (mqd_rotorperiod);\n')
            
        # Pulse list
        delta_MQMAS.write(conversion_pulses)
        delta_MQMAS.write('	tau;' + '\n')
        
        # Position of the single-quantum delay depending on the spin and the transition
        if fmc.get_splitt1_pos(self.param['spin'],self.param['transition']):
            delta_MQMAS.write('	 (' + splitt1[2] + '* initial_t1) ystep (sqd);' + '\n' + 
        '	sel_180_pulse, (obs.gate, obs.phs.obs_phs_CTsel, obs.amp.sel_amp_pulse, obs.atn.x_atn);' + '\n') 
        else:
            delta_MQMAS.write('	sel_180_pulse, (obs.gate, obs.phs.obs_phs_CTsel, obs.amp.sel_amp_pulse, obs.atn.x_atn);' + '\n' +
        '	(' + splitt1[2] + '* initial_t1) ystep (sqd);' + '\n') 
        
        # End of the pulse program
        delta_MQMAS.write('	acq (dead_time, delay, obs_phs_phase_acq);' + '\n' +
        '	when Decoupling do ' + '\n' +
        '	off (irr.gate);' + '\n' +
        '	end when; ' + '\n' +
        'end pulse; ')
        
        # Optimisation comments
        delta_MQMAS.write(self.get_final_pp_comments('jeol'))
        
        # ---------------------------------------------------------------------
        print('MQMAS pulse program for Delta successfully exported')
        delta_MQMAS.close()
        
    # -------------------------------------------------------------------------
    # ----- Export GUI
    
    def export_all(self):
        '''Select folder for the output of things and export'''
        self.param['export_filepath'] = os.path.normpath(fmc.select_directory('Select output folder')[0])
        self.export_all_files()
        
    def export_all_files(self, exportfilepath = None):
        ''' Save all possible outputs in the output folder'''
        
        # Set the output file if no input argument is given
        if exportfilepath == None:
            exportfilepath = self.param['export_filepath']
        
        # Control and export
        if (exportfilepath != None and os.path.exists(exportfilepath)):
            print('')
            print('---- Exporting results ----')
            self.output['existing_allexport_filepath'] = True 
            
            # General things
            self.export_text_report((self.param['simpson_basename'] + '_text_report'),
                                      exportfilepath)
            self.plot_allsteps_fig()
            self.figures['allsteps_fig'].savefig(exportfilepath +
                        os.sep + self.param['simpson_basename'] + '_figure.pdf')
            # Export Bruker
            self.export_topspin_MQMAS((self.param['simpson_basename'] + '_topspin_MQMAS'),
                                      exportfilepath)
            self.export_topspin_pulselist((self.param['simpson_basename'] + '_topspin_pulselist'),
                                          exportfilepath)
            self.export_topspin_MQF((self.param['simpson_basename'] + '_topspin_MQF'),
                                          exportfilepath)
            self.export_topspin_MQCPMGMAS((self.param['simpson_basename'] + '_topspin_MQCPMGMAS'),
                                          exportfilepath)
            # Export Jeol
            self.export_delta_MQMAS((self.param['simpson_basename'] + '_jeol_MQMAS'),
                                      exportfilepath)
            self.export_delta_MQF((self.param['simpson_basename'] + '_jeol_MQF'),
                                      exportfilepath)
            self.export_delta_pulselist((self.param['simpson_basename'] + '_jeol_pulselist'),
                                      exportfilepath)
            
            # Export save file (.pkl)
            self.savestate((self.param['simpson_basename'] + '.pkl'),exportfilepath)
            
            # Open folder at the end if executing
            try:
                fmc.openFolder(os.path.normpath(exportfilepath))
            except:
                pass
            
        else:
            self.output['existing_allexport_filepath'] = False 
            if self.param['console_warning']:
                print('Warning: Empty or invalid export filepath. No file exported.')
            
    
    def save(self):
        ''' Save the current state of this class instance as a serialised version of the class instance'''
        # Get file path
        inputpath = fmc.save_dialog(Caption = 'Saving current optimisation',
                    Filter_tupple = (".pkl file (*.pkl)","All files (*)"))
        
        # Write file
        if inputpath != None:
            filepath, filename = os.path.split(inputpath[0])
            self.savestate(filename, filepath)
            print('Instance', filename ,'successfully saved as', filepath)
        
        
    def load(self):
        ''' Load a saved state of this class instance '''
        inputpath = fmc.open_dialog(Caption = 'Loading saved optimisation',
                    Filter_tupple = (".pkl file (*.pkl)","All files (*)"))
        
        if inputpath != None:
            filepath, filename = os.path.split(inputpath[0])
            self.loadstate(filename, filepath)
            print('File', filename ,'successfully loaded')
            
            
    # ----------------------------
    def text_report(self):
        ''' Export the report of the optimisation '''
        # Get file path
        inputpath = fmc.save_dialog(Caption = 'Save optimisation report',
                    Filter_tupple = ("Text file (*.txt)","All files (*)"))
        if inputpath != None:
            # Write file
            filepath, filename = os.path.split(inputpath[0])
            self.export_text_report(filename, filepath)
            print('File', filename ,'successfully created in', filepath)
    
    
    def all_steps_figure(self):
        ''' Export the optimisation picture '''
        # Plot figure
        self.plot_allsteps_fig()
        # Get file path
        inputpath = fmc.save_dialog(Caption = 'Save optimisation report',
                    Filter_tupple = fmc.fig_tuple() )
        # Write file
        if inputpath != None:
            filepath, filename = os.path.split(inputpath[0])
            self.figures['allsteps_fig'].savefig(inputpath[0])
            print('File', filename ,'successfully created in', filepath)
    
    
    def step_figure(self,nos = 1):
        ''' Saves the figure containing all optimisations from the step "nos" with the name
        'filename' (with extension) at the location 'filepath'.
        '''
        
        if self.param['debug_figure'] and (nos-1) < self.output['no_FAMN_steps']:
            # Get file path
            inputpath = fmc.save_dialog(Caption = 'Save optimisation report',
                        Filter_tupple = fmc.fig_tuple() )
            if inputpath != None:
                filepath, filename = os.path.split(inputpath[0])
                self.FAMNsteps[nos-1].figures['opt'].savefig(inputpath[0])
                print('File', filename ,'successfully created in', filepath)
                
        
        # plt.savefig()
    
        
    # ----------------------------
    def topspin_pulselist(self):
        ''' Open dialog and save pulse list for Bruker Topspin '''
        # Get file path
        inputpath = fmc.save_dialog(Caption = 'Save pulse list for Bruker Topspin',
                    Filter_tupple = ("Text file (*.txt)","All files (*)"))
        if inputpath != None:
            # Write file
            filepath, filename = os.path.split(inputpath[0])
            self.export_topspin_pulselist(filename, filepath)
            print('File', filename ,'successfully created in', filepath)
    
    def topspin_mqf(self):
        ''' Open dialog and save MQF pulse program for Bruker Topspin '''
        # Get file path
        inputpath = fmc.save_dialog(Caption = 'Save MQF pulse program for Bruker Topspin',
                    Filter_tupple = ("All files (*)",))
        if inputpath != None:
            # Write file
            filepath, filename = os.path.split(inputpath[0])
            self.export_topspin_MQF(filename, filepath)
            print('File', filename ,'successfully created in', filepath)
    
    def topspin_mqmas(self):
        ''' Open dialog and save MQMAS pulse program for Bruker Topspin '''
        # Get file path
        inputpath = fmc.save_dialog(Caption = 'Save MQMAS pulse program for Bruker Topspin',
                    Filter_tupple = ("All files (*)",))
        if inputpath != None:
            # Write file
            filepath, filename = os.path.split(inputpath[0])
            self.export_topspin_MQF(filename, filepath)
            print('File', filename ,'successfully created in', filepath)
        
    def topspin_mqcpmgmas(self):
        ''' Open dialog and save MQ-CPMG-MAS pulse program for Bruker Topspin '''
        # Get file path
        inputpath = fmc.save_dialog(Caption = 'Save MQ-CPMG-MAS pulse program for Bruker Topspin',
                    Filter_tupple = ("All files (*)",))
        if inputpath != None:
            # Write file
            filepath, filename = os.path.split(inputpath[0])
            self.export_topspin_MQCPMGMAS(filename, filepath)
            print('File', filename ,'successfully created in', filepath)
    
    # ----------------------------
    def delta_pulselist(self):
        ''' Open dialog and save pulse list for Jeol Delta '''
        # Get file path
        inputpath = fmc.save_dialog(Caption = 'Save pulse list for Bruker Topspin',
                    Filter_tupple = ("Text file (*.txt)","All files (*)"))
        if inputpath != None:
            # Write file
            filepath, filename = os.path.split(inputpath[0])
            self.export_topspin_pulselist(filename, filepath)
            print('File', filename ,'successfully created in', filepath)
    
    def delta_mqf(self):
        ''' Open dialog and save MQF pulse program for Jeol Delta '''
        # Get file path
        inputpath = fmc.save_dialog(Caption = 'Save MQF pulse program for Bruker Topspin',
                    Filter_tupple = ("All files (*)",))
        if inputpath != None:
            # Write file
            filepath, filename = os.path.split(inputpath[0])
            self.export_topspin_MQF(filename, filepath)
            print('File', filename ,'successfully created in', filepath)
    
    def delta_mqmas(self):
        ''' Open dialog and save MQMAS pulse program for Jeol Delta '''
        # Get file path
        inputpath = fmc.save_dialog(Caption = 'Save MQMAS pulse program FAM-N optimisation',
                    Filter_tupple = ("All files (*)",))
        if inputpath != None:
            # Write file
            filepath, filename = os.path.split(inputpath[0])
            self.export_topspin_MQF(filename, filepath)
            print('File', filename ,'successfully created in', filepath)

# -----------------------------------------------------------------------------
# End of class FAMN_main





#import gui_main

# -----------------------------------------------------------------------------
#class FAMN_gui(QtCore, QtGui, QtWidgets ,gui_main.Ui_MainWindow):  
#    'GUI for the FAM-N optimisation class'
#    
#    def __init__(self, parent=None):
#        super(FAMN_gui, self).__init__(parent)
#        self.setupUi

# -----------------------------------------------------------------------------
#def show_gui():
#    app = QtGui.QApplication(sys.argv)
#    form = FAMN_gui()
#    form.show()
#    app.exec_()


if __name__ == '__main__':
    ''' Start GUI if this function is executed'''
    # show_gui()
    pass