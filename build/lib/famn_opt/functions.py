# -*- coding: utf-8 -*-
"""
Created on Tue Apr 16 17:40:44 2019

@author: Henri
"""

# À faire: add dynamic pulse: increase pulse lenght if maximum is reached

# Import packages
import sys, os
import subprocess
import numpy 
import lxml.etree as xmlp # To read xml files
import matplotlib.pyplot as plt # To plot curves
import pickle # To save classes

# GUI packages
from PyQt5 import QtCore, QtGui, QtWidgets
# -----------------------------------------------------------------------------
#
# Commun functions
#
# -----------------------------------------------------------------------------

if sys.platform == 'darwin':
    def openFolder(path):
        subprocess.check_call(['open', '--', path])
elif sys.platform == 'linux2':
    def openFolder(path):
        subprocess.check_call(['xdg-open', '--', path])
elif sys.platform == 'win32':
    def openFolder(path):
        subprocess.check_call(['explorer', path])
        

def write_famn_console_input(filename, filepath):
    ''' Write a generic console input file for famn_opt '''
 
    if filepath == '':
        inputfile = open(change_ext(filename,'py') ,'w', encoding='utf-8')
    else:
        inputfile = open(os.path.normpath(filepath) + os.sep + change_ext(filename,'py') ,'w', encoding='utf-8')
    
    filenamebase , _ = os.path.splitext(filename) # Remove extension for the inside of the file
    inputfile.write('''# -*- coding: utf-8 -*-
"""
Console mode execution of the FAM-N optimisation program.
Fill in the prameters execute this script (console)
or press the run button/F5 (Spyder/pyCharm/VSstudio…)
"""

# Initialise the script
from famn_opt import FAMN_main
''' + filenamebase + ''' = None
''' + filenamebase + ''' = FAMN_main()
''' + filenamebase + '''.param['simpson_basename'] = "''' + filenamebase + '''"

# ----------------------------------------------------------------------------
# (Optionnal) Path where all output files are exported
''' + filenamebase + '''.param['export_filepath'] = ''

# ----------------------------------------------------------------------------
# Input parameters

# Spin system input parameters
''' + filenamebase + '''.param['nucleus'] = '23Na' # Nucleus, as written in Simpson
''' + filenamebase + '''.param['proton_larmor_frequency'] = 600e6 # Proton Larmor frequency (Hz) 
''' + filenamebase + '''.param['MAS'] = 12500 # Magic angle spinning rate (Hz) 
''' + filenamebase + '''.param['RF'] = 100000 # Radio-frequency field ν1 (Hz)  
''' + filenamebase + '''.param['CQ'] = 1.0e6 # Quadrupolar coupling constant CQ (Hz) 
''' + filenamebase + '''.param['ηQ'] = 0 # Asymmetry parameter EtaQ 

# ----------------------------------------------------------------------------
# Transfer parameters

# To fill for a standard FAM-N optimisation
''' + filenamebase + '''.param['transition'] = 3 # Transition of the FAM-N method.

# To fill in for complex FAM-N optimisation
''' + filenamebase + '''.param['complex_optimisation'] = ''
''' + filenamebase + '''.param['single_pulse_nutation_list'] = '' # If complex optimisation, string list of the single pulse nutation angles for each step

# ---------------------------------------------
# Crystal file
''' + filenamebase + '''.param['cryst_gen_method'] = 'golden' # Method used to genenerate the crystal file. Can be 'golden' or 'zcw'
''' + filenamebase + '''.param['no_cryst'] = 66 # Number of crystallites
''' + filenamebase + '''.param['no_γ'] = 4 # Number of γ-angles for the SIMPSON optimisation

# ---------------------------------------------
# Advanced user parameters.

''' + filenamebase + '''.param['min_improvement'] = 1.1 # Minimum improvement between two subsequent pulses
                                    # for the simulation to be allowed to add additionnal pulses
''' + filenamebase + '''.param['max_nopulses'] = 30 # Maximum number of pulses that the program is allowed to output.
''' + filenamebase + '''.param['nop_360'] = 180 # Number of points par 360° nutation.
''' + filenamebase + '''.param['min_time_increment'] = 1e-8 # Minimum time increment (seconds)
''' + filenamebase + '''.param['maxdt_ratio'] = 20 # Ratio between MAS period and maxdt in SIMPSON.
''' + filenamebase + '''.param['num_cores'] = 0 # Number of cores on which the optimisation is programmed

# ----------------------------------------------------------------------------
# Execution of the code
''' + filenamebase + '''.initialise()
''' + filenamebase + '''.FAMN_loop()
''' + filenamebase + '''.finalise()
''' + filenamebase + '''.export_all_files() # Export automatically if the output path is specified

# ----------------------------------------------------------------------------
''')
    inputfile.close()

# ----------------------------------
# --------------- Utility functions

def is_simpson_installed():
    ''' Verifies if simpson is properly installed. '''
    try:
        subprocess.call('simpson')
    except FileNotFoundError:
        return False
    else:
        return True

    
def chemical_symbol_first(chaine):
    ''' Makes sure that the element symbol comes before the number.
    This is only to comply with the xml format that prohibits string starting
    with a number.
    '''
    valeur_num = ''
    valeur_car = ''
    for truc in range(len(chaine)):
        if chaine[truc] in '0123456789':
           valeur_num += chaine[truc]
        else:
           valeur_car += chaine[truc]
            
    return valeur_car + valeur_num

def get_modulepath():
    ''' Returns the filepath of the module'''
    return os.path.dirname(os.path.abspath(__file__))

def read_spin_from_xml(nucleus):
    ''' Nb: element number must go first: e.g. : 'H1', 'Na23' '''
    return xmlp.parse(get_modulepath() + os.sep + 'spininfo.xml').xpath(nucleus.lower() + '/spin')[0].text

def colour_in_plot(integer):
    ''' Return a colour of the plot as a function of the modulo of the input integer  '''
    mod_integer = numpy.mod(integer,6)
    if mod_integer == 0:
        return 'black'
    elif mod_integer == 1:
        return '#2a7fff' # Some blue-green
    elif mod_integer == 2:
        return '#b10000' # Some dark-red
    elif mod_integer == 3:
        return '#55d400' # Some green
    elif mod_integer == 4:
        return '#be00be' # Some purple
    elif mod_integer == 5:
        return '#d59265' # Some orange


def find_nearest_idx(array, value):
    '''Find the index in 'array' that is not the first value of array.'''
    array = numpy.asarray(array)
    idx = (numpy.abs(array - value)).argmin()
#    print(idx)
    return int(idx)


def display_plus(number):
    ''' Display number as string with visible plus sign'''
    if number < 0:
        return ('−' + str(numpy.abs(number)))
    elif number > 0:
        return ('+' + str(numpy.abs(number)))
    else:
        return '0'
    
def to_full(array):
    ''' Fill a numpy array with uneven values with numpy.nan values'''
    output = numpy.full((len(array), max(map(len, array))), numpy.nan)
    for truc, row in enumerate(array):
        output[truc, :len(row)] = row
    return output

# --------------------------------------------------------------------------
# --------------- Load and save classes  
    
def change_ext(in_filename,in_ext):
    ''' Change the extension of filename to 'ext' if filename had no extension '''
    ext = in_ext.split('.')[-1] # Remove the point if any; do nothing if not.
    path,file = os.path.split(in_filename) # Remove path if any
    filename, extension = os.path.splitext(file) # Separate filename from extension
    if extension == '': # If no extension
        return (str(path) + str(filename) + '.' + str(ext))
    else:
        return in_filename


def save_pkl(saved_object, filename, filepath = ''):
    ''' Save class filename located in filepath'''
    outputfile = open(filepath + os.path.basename(filename) + '.pkl', 'wb')
    pickle.dump(saved_object, outputfile, protocol=pickle.HIGHEST_PROTOCOL)
    outputfile.close()
    
    
def load_pkl(filename, filepath = ''):
    ''' Load filename (no extension) located in 'filepath' '''
    inputfile = open(filepath + os.path.basename(filename) + '.pkl', 'rb')
    return pickle.load(inputfile)
    inputfile.close()
    
    
# -------------------------------------------------------------------------
# --------------- NMR/SIMPSON functions

def one_element_density(dim_density,matrix_element,amplitude = 1):
    ''' Get an initial density matrix that is all zeroes except one matrix element with
    coordinates matrix_element[0] and matrix_element[1] and amplitute 'amplitude' (default = 1)
    '''
    density = numpy.zeros((dim_density,dim_density), dtype = 'complex128')
    density[matrix_element[0]-1,matrix_element[1]-1] = amplitude
    return density


def get_dim_from_spin(spin):
    ''' Get the dimension of the density matrix from the spin number written as a string.'''
    if spin == '3/2':
        return 4
    elif spin =='5/2':
        return 6
    elif spin =='7/2':
        return 8
    elif spin =='9/2':
        return 10

def get_element_from_coherence(coherence,dim_density):
    ''' Returns the matrix element corresponding to a coherence for a given spin.'''
    if dim_density == 4:
        if coherence == 1:
            return [2,3]
        elif coherence ==-1:
            return [3,2]
        elif coherence == 3:
            return [1,4]
        elif coherence == -3:
            return [4,1]
    elif dim_density == 6:
        if coherence == 1:
            return [3,4]
        elif coherence ==-1:
            return [4,3]
        elif coherence == 3:
            return [2,5]
        elif coherence == -3:
            return [5,2]
        elif coherence == 5:
            return [1,6]
        elif coherence == -5:
            return [6,1]
    elif dim_density == 8:
        if coherence == 1:
            return [4,5]
        elif coherence ==-1:
            return [5,4]
        elif coherence == 3:
            return [3,6]
        elif coherence == -3:
            return [6,3]
        elif coherence == 5:
            return [2,7]
        elif coherence == -5:
            return [7,2]
        elif coherence == 7:
            return [1,8]
        elif coherence == -7:
            return [8,1]
    elif dim_density == 10:
        if coherence == 1:
            return [5,6]
        elif coherence ==-1:
            return [6,5]
        elif coherence == 3:
            return [4,7]
        elif coherence == -3:
            return [7,4]
        elif coherence == 5:
            return [3,8]
        elif coherence == -5:
            return [8,3]
        elif coherence == 7:
            return [2,9]
        elif coherence == -7:
            return [9,2]
        elif coherence == 9:
            return [1,10]
        elif coherence == -9:
            return [10,1]
    else: return None
        
            
            
def pauli_matrix(dim_density,direction):
    ''' Get Pauli's spin matrix for each spin/direction '''
    if dim_density == 2:
        if direction == 'x':
            return numpy.array([[0,1],
                    [1,0]],
                    dtype = 'complex128')
        if direction == 'y': 
            return numpy.array([[0,-1j],
                    [1j,0]],
                    dtype = 'complex128')
        if direction == 'z':
            return numpy.array([[1,0],
                    [0,-1]],
                    dtype = 'complex128')
            
    if dim_density == 4:
        if direction == 'x':
            return numpy.array([[0,numpy.sqrt(3),0,0],
                    [numpy.sqrt(3),0,2,0],
                    [0,2,0,numpy.sqrt(3)],
                    [0,0,numpy.sqrt(3),0]],
                    dtype = 'complex128')
        if direction == 'y': 
            return numpy.array([[0,-1j*numpy.sqrt(3),0,0],
                    [1j*numpy.sqrt(3),0,2j,0],
                    [0,2j,0,-1j*numpy.sqrt(3)],
                    [0,0,1j*numpy.sqrt(3),0]]
                    , dtype = 'complex128')
        if direction == 'z':
            return numpy.array([[3,0,0,0],
                    [0,1,0,0],
                    [0,0,-1,0],
                    [0,0,0,-3]]
                    , dtype = 'complex128')
    if dim_density == 6:
        if direction == 'x':
            return numpy.array([[0,numpy.sqrt(5),0,0,0,0],
                    [numpy.sqrt(5),0,numpy.sqrt(8),0,0,0],
                    [0,numpy.sqrt(8),0,3,0,0],
                    [0,0,3,0,numpy.sqrt(8),0],
                    [0,0,0,numpy.sqrt(8),0,numpy.sqrt(5)],
                    [0,0,0,0,numpy.sqrt(5),0]]
                    , dtype = 'complex128')
        if direction == 'y': 
            return numpy.array([[0,-1j*numpy.sqrt(5),0,0,0,0],
                    [1j*numpy.sqrt(5),0,-1j*numpy.sqrt(8),0,0,0],
                    [0,1j*numpy.sqrt(8),0,-1j*3,0,0],
                    [0,0,3,0,1j*numpy.sqrt(8),0],
                    [0,0,0,1j*numpy.sqrt(8),0,1j*numpy.sqrt(5)],
                    [0,0,0,0,1j*numpy.sqrt(5),0]]
                    , dtype = 'complex128')
        if direction == 'z':
            return numpy.array([[5,0,0,0,0,0],
                    [0,3,0,0,0,0],
                    [0,0,1,0,0,0],
                    [0,0,0,-1,0,0],
                    [0,0,0,0,-3,0],
                    [0,0,0,0,0,-5]]
                    , dtype = 'complex128')
    if dim_density == 8:
        if direction == 'x':
            return numpy.array([[0,numpy.sqrt(7),0,0,0,0,0,0],
                    [numpy.sqrt(7),0,2*numpy.sqrt(3),0,0,0,0,0],
                    [0,2*numpy.sqrt(3),0,numpy.sqrt(15),0,0,0,0],
                    [0,0,numpy.sqrt(15),0,4,0,0,0],
                    [0,0,0,4,0,numpy.sqrt(15),0,0],
                    [0,0,0,0,numpy.sqrt(15),0,2*numpy.sqrt(3),0],
                    [0,0,0,0,0,2*numpy.sqrt(3),0,numpy.sqrt(7)],
                    [0,0,0,0,0,0,numpy.sqrt(7),0]]
                    , dtype = 'complex128')
        if direction == 'y': 
            return numpy.array([[0,-1j*numpy.sqrt(7),0,0,0,0,0,0],
                    [1j*numpy.sqrt(7),0,-1j*2*numpy.sqrt(3),0,0,0,0,0],
                    [0,1j*2*numpy.sqrt(3),0,-1j*numpy.sqrt(15),0,0,0,0],
                    [0,0,1j*numpy.sqrt(15),0,-1j*4,0,0,0],
                    [0,0,0,1j*4,0,-1j*numpy.sqrt(15),0,0],
                    [0,0,0,0,1j*numpy.sqrt(15),0,-1j*2*numpy.sqrt(3),0],
                    [0,0,0,0,0,1j*2*numpy.sqrt(3),0,-1j*numpy.sqrt(7)],
                    [0,0,0,0,0,0,1j*numpy.sqrt(7),0]]
                    , dtype = 'complex128')
        if direction == 'z':
            return numpy.array([[7,0,0,0,0,0,0,0],
                    [0,5,0,0,0,0,0,0],
                    [0,0,3,0,0,0,0,0],
                    [0,0,0,1,0,0,0,0],
                    [0,0,0,0,-1,0,0,0],
                    [0,0,0,0,0,-3,0,0],
                    [0,0,0,0,0,0,-5,0],
                    [0,0,0,0,0,0,0,-7]]
                    , dtype = 'complex128')
    if dim_density == 10:
        if direction == 'x':
            return numpy.array([[0,3,0,0,0,0,0,0,0,0],
                    [3,0,4,0,0,0,0,0,0,0],
                    [0,4,0,numpy.sqrt(21),0,0,0,0,0,0],
                    [0,0,numpy.sqrt(21),0,2*numpy.sqrt(6),0,0,0,0,0],
                    [0,0,0,2*numpy.sqrt(6),0,5,0,0,0,0],
                    [0,0,0,0,5,0,2*numpy.sqrt(6),0,0,0],
                    [0,0,0,0,0,2*numpy.sqrt(6),0,numpy.sqrt(21),0,0],
                    [0,0,0,0,0,0,2*numpy.sqrt(21),0,4,0],
                    [0,0,0,0,0,0,0,4,0,3],
                    [0,0,0,0,0,0,0,0,3,0]]
                    , dtype = 'complex128')
        if direction == 'y': 
            return numpy.array([[0,-1j*3,0,0,0,0,0,0,0,0],
                    [1j*3,0,-1j*4,0,0,0,0,0,0,0],
                    [0,1j*4,0,-1j*numpy.sqrt(21),0,0,0,0,0,0],
                    [0,0,1j*numpy.sqrt(21),0,-1j*2*numpy.sqrt(6),0,0,0,0,0],
                    [0,0,0,1j*2*numpy.sqrt(6),0,-1j*5,0,0,0,0],
                    [0,0,0,0,1j*5,0,-1j*2*numpy.sqrt(6),0,0,0],
                    [0,0,0,0,0,1j*2*numpy.sqrt(6),0,-1j*numpy.sqrt(21),0,0],
                    [0,0,0,0,0,0,1j*2*numpy.sqrt(21),0,-1j*4,0],
                    [0,0,0,0,0,0,0,1j*4,0,-1j*3],
                    [0,0,0,0,0,0,0,0,1j*3,0]]
                    , dtype = 'complex128')
        if direction == 'z':
            return numpy.array([[9,0,0,0,0,0,0,0,0,0],
                    [0,7,0,0,0,0,0,0,0,0],
                    [0,0,5,0,0,0,0,0,0,0],
                    [0,0,0,3,0,0,0,0,0,0],
                    [0,0,0,0,1,0,0,0,0,0],
                    [0,0,0,0,0,-1,0,0,0,0],
                    [0,0,0,0,0,0,-3,0,0,0],
                    [0,0,0,0,0,0,0,-5,0,0],
                    [0,0,0,0,0,0,0,0,-7,0],
                    [0,0,0,0,0,0,0,0,0,-9]]
                    , dtype = 'complex128')


def get_Larmor_frequency(nucleus,proton_larmor_frequency = 0, static_magnetic_field = 0 ):
    ''' Get the Larmor frequency from the xml file 'spininfo.xml' '''

    if proton_larmor_frequency == 0 and static_magnetic_field == 0:
        raise ValueError('Larmor frequency and the static magnetic field provided '
                         'must not be both zero.')
    gyromag = 1.0e6 * float(xmlp.parse(get_modulepath() + os.sep + 'spininfo.xml').xpath(chemical_symbol_first(nucleus).lower() + '/gamma')[0].text)
    if proton_larmor_frequency != 0:
        return float(proton_larmor_frequency * gyromag / 42.57747892e6)
    elif static_magnetic_field != 0:
        return float(gyromag * static_magnetic_field)


def get_PQ(CQ,etaQ):
    ''' Get the quadrupolar product from the parameters '''
    return CQ*numpy.sqrt(1 + ( (etaQ**2)/(3) ) )


def get_QIS(nucleus,spin,nucleus_larmor_frequency,CQ,etaQ):
    ''' Function to retrieve the quadruporal-induced shift in the simulations'''
    return (-1)*(3/10)*( ((spin*(spin+1)-(3/4))) / (2*spin*(2*spin-1))**2 )*\
            ( ( pow(get_PQ(CQ,etaQ),2) ) / ( nucleus_larmor_frequency ) )


def get_density_complex_famn(element,dim_density):
    try: # Case this thing is integer (i.e. is a coherence order)
        
        element = get_element_from_coherence(int(element),dim_density)
        
    except:
        # Case it is a direction operator
        if element in ['x','y','z']:
            return pauli_matrix(dim_density,element)
        # Case it is an operator
        else:
            return one_element_density(dim_density,eval(element))
        
    else:
        # Verifies that the transition is possible for the input spin
        if element != None:
            return one_element_density(dim_density,element)
        else: # Returns nan
            return  numpy.full((dim_density,dim_density),numpy.nan)


def get_detect_element(idx, in_density, detect_density):
    ''' Returns the detect element defined by 'detect_density' from the matrix element
    'in_density' with position 'idx'. Detect density most be the transconjugate of the
    element that you actually want to detect.
    '''
    return numpy.trace(numpy.dot(in_density[idx],detect_density), dtype = 'complex128')


def get_all_detect_element(in_density, detect_density , nop = 0):
    ''' Returns a complex128 numpy array with the first 'nop' detect elements of
    matrix 'in_density' '''
    if nop == 0:
        nop = len(in_density)
    out_array = numpy.zeros(nop, dtype = 'complex128')  
    for truc in range(nop):
        out_array[truc] = get_detect_element(truc, in_density, detect_density)
    return out_array


def parse_complex_famn(string,dim_density):
    ''' Function that converts the text input corresponding to complex FAM-N to
    a list of start and detect density matrices
    '''
    
    # Pre-allocate the output list
    input_list = []
    output_list = []
    
    # Separate with ';'
    temp1 = string.split(';')
    for truc1 in range(len(temp1)):
        
        temp2 = temp1[truc1].split('>')
        # ---------------
        len_temp2 = len(temp2)
        if len_temp2 > 1:
            # --------------------------------
            # Set the first element
            input_list.append(get_density_complex_famn(temp2[0],dim_density))
            # --------------------------------
            # Set middle elements if any
            for truc2 in range(1,(len_temp2 - 1)):
                output_list.append(numpy.transpose(get_density_complex_famn(temp2[truc2],dim_density)))
                # input_list.append(None)
                input_list.append(numpy.full((dim_density,dim_density),None)) # None means "take final density matrix of last fam-n"
            # --------------------------------
            # Set the end element
            output_list.append(numpy.transpose(get_density_complex_famn(temp2[len_temp2-1],dim_density)))
        else:
            raise ValueError('The input string ' + temp2 + ' contains only one element')
        
    # --------------
    return input_list, output_list
        

def array_to_simpson(matrix, dim = 0):
    ''' Turn a numpy array into a SIMPSON matrix'''
    
    # Get the dimension of the input array
    if dim == 0:
        dim = int(numpy.sqrt(len(matrix)))
    
    # Sort the array
    string = '{\n'
    for i in range(dim):
        string += '{'
        for j in range(dim):
            string += '{' + str(numpy.real(matrix[ i , j ])) + ' ' +\
                        str(numpy.imag(matrix[ i , j  ])) + '}' + ' '
        string += '}' + '\n'
    # ---------------------------
    return (string + '}'  + '\n')


def readin_simpson_out(filepath,nop,dim_density,size_density):
    ''' Reads in the content of SIMPSON output (.out) file'''

    # Open .out simpson file
    simpson_out_str = open(filepath).read().splitlines()
    
    # Parse the text output file of SIMPSON into a numpy array.
    out_array = numpy.zeros((nop,dim_density,dim_density), dtype = 'complex128')
    
    for truc1 in range(nop):
        for truc2 in range(dim_density):
            for truc3 in range(dim_density):
                temp = simpson_out_str[ 5+truc1*size_density+truc2*dim_density + truc3 ].split(' ')
                out_array[truc1,truc2,truc3] = float(temp[0]) + 1j * float(temp[1])
    return out_array

def gen_detect_elements(noel):
    ''' Generate all the detect elements of a spin system.
    'noel' is the dimension of the density matrix '''
    for truc1 in range(noel):
        for truc2 in range(noel):
            yield '{' + str(truc2+1) + ' ' + str(truc1+1) + '}'



def get_maxima_and_positions(data):
    ''' Get the position of the local maxima of "data" '''
    return (numpy.diff(numpy.sign(numpy.diff(data))) < 0).nonzero()[0] + 1


def get_pulse_duration(max_array,inv_array,nopulses):
    ''' Returns an array containing the pulse durations (in number of points) from
    the output arrays of FAMN_steps.'''
    out_array = numpy.zeros(nopulses, dtype = 'int32')
    for truc in range(nopulses-1):
        out_array[truc] = max_array[truc] + inv_array[truc]
    out_array[nopulses-1] = max_array[nopulses-1]
    return out_array

 # -------------------------------------------------------------------------
 # --------------- Timing functions
 
def start_timer():
    ''' Module that calculates the total time of the execution of the program.
    NB: Optionnal. If the module 'time' in not found, this will be skipped.'''
    try:
        global time
        import time
        return time.time()
    except ModuleNotFoundError:
        print('The module "time" was not found. The executation time of this program will not be displayed.')
        return 0


def stop_timer(ref_time):
    ''' Second part, to use after set_timer.  '''
    try:
        ref_time -= time.time()
        return numpy.abs(ref_time)
    except NameError:
        return 0


def convert_into_smhd(time):
    '''A function that turns a given time in seconds into seconds, minutes,
    hours and days'''
    seconds = numpy.mod(time,60)
    minutes = numpy.mod((time - seconds)/60,60)
    hours = numpy.mod((time - minutes*60 - seconds)/(3600), 24)
    days = (time - 3600*hours - minutes*60 - seconds)/(86400)
    return (seconds,minutes,hours,days)


def time_to_string(time):
    '''Display execution time '''
    time_tupple = convert_into_smhd(time)
    time_string = str(numpy.round(time_tupple[0],2)) + ' s '
    #
    if time_tupple[1] != 0:
        time_string = str(time_tupple[1]) + ' min ' + time_string
        if time_tupple[2] != 0:
            time_string = str(time_tupple[2]) + ' hours ' + time_string
            if time_tupple[3] != 0:
                time_string = str(time_tupple[3]) + ' days ' + time_string
    return time_string 

# -----------------------------------------------------------------------------
# --------------- Monte-Carlo functions
    
def golden_sphere(nop, step = (1 + numpy.sqrt(5)), offset = 0.5, cartesians = False ):
    ''' So-called "golden spiral" method to genenertate 'nop' points roughly evenly separated on a sphere.
    step: step used in theta
    offset (between 0 and 1, can be above)
    cartesian = True return cartesian coordinates
    
    Source:
    https://stackoverflow.com/questions/9600801/evenly-distributing-n-points-on-a-sphere
    '''
    indices = numpy.arange(0, nop, dtype=float) + offset
    
    phi = numpy.arccos(1 - 2*indices/nop)
    theta = numpy.pi * step * indices
    # weight  = numpy.ones(samples)/samples
    
    if cartesians:
        return numpy.cos(theta) * numpy.sin(phi), numpy.sin(theta) * numpy.sin(phi), numpy.cos(phi) #, weight
    else:
        return numpy.mod(theta,2*numpy.pi) , phi #, weight
    

def fib(n):
        """
        From http://enumpy.literateprograms.org/Fibonacci_numbers_(Python)
        Returns (n+2)-th and n-th fibonacci number
        """
        if n > 90:
                raise ValueError("Can't represent higher numbers!")
        start = numpy.array([[1,1],[1,0]],dtype='int64')
        temp = start[:]
        for i in range(n):
                temp = numpy.dot(start,temp)
        return temp[0,0] ,temp[0,1],temp[1,1]


def zcw_angles(m):
        """
        Getting N=fibonacchi(m) numbers of ZCW angles
        
        Source:
        https://www.mail-archive.com/numpy-discussion@scipy.org/msg19081.html
        http://dx.doi.org/10.1002/cmr.a.10065
        Computer simulations in solid-state NMR. III. Powder averaging
        """

        samples, fib_1, fib_2 = fib(m)
        #fib_2 = fib_1

        js = numpy.arange(samples, dtype='float64')/samples
        
        # c_full = (1.,2.,1.)
        # c_hemi = (-1.,1.,1.)
        # c_oct  = (2.,1.,8.)

        c = (1.,2.,1.)
        
        j_samples = fib_2 * js
        # alpha
        phi = 2*numpy.pi/c[2] * numpy.mod(j_samples,1.0)
        # beta
        theta = numpy.arccos(c[0]*(c[1]*numpy.mod(js,1.0) - 1))
        # weights
        weight = numpy.ones(samples)/samples
        return phi,theta,weight
    
    
def naive_random_sphere_generator(nop):
    ''' Naive random angle generator on a sphere '''
    theta = numpy.zeros(nop)
    phi = numpy.zeros(nop)
    
    for truc in range(nop):
        norm = float(2)
        # ----------------
        while norm > 1.0 :
        
            x = 2*numpy.random.random()-1
            y = 2*numpy.random.random()-1
            z = 2*numpy.random.random()-1
            
            norm = numpy.sqrt(x**2+y**2+z**2)
        # ----------------
        try:
            theta[truc] = numpy.arctan(y/x)
        except ZeroDivisionError:
            theta[truc] = numpy.random.choice([-1,1]) * numpy.pi/2
        finally:
            theta[truc] +=  numpy.random.choice([0,numpy.pi])
        # ----------------
        try:
            phi[truc] = numpy.arctan(numpy.sqrt(x**2+y**2)/z)
        except ZeroDivisionError:
            phi[truc] = numpy.random.choice([-1,1]) * numpy.pi/2
    # -----------------------------------------------------------
    theta = numpy.mod(theta, 2 * numpy.pi)
    phi = numpy.mod(phi, numpy.pi)
    # weight = numpy.ones(nop)/nop
    
    return theta, phi #, weight

def random_sphere_generator(ndp):
    ''' Random angles sphere generator'''
    test = numpy.zeros((2,ndp))
    for truc in range(ndp):
        test[0][truc] = 2*numpy.pi*numpy.random.random()
        test[1][truc] = numpy.arccos(2*numpy.random.random()-1)
    return test

# -----------------------------------------------------------------------------
# --------------- Functions for pulse program export
    
def get_bruker_phase_alternation(test):
    ''' Return '2' if test is True, 3 if not'''
    if test:
        return '2'
    else:
        return '3'
    
def get_splitt1_pos(spin,transition):
    ''' Returns True if the single-quantum delay is located after the CT 180° '''
    
    if transition == 3 and  spin == '3/2':
        return True
    elif transition == 5 and spin == '5/2':
        return True
    elif transition == 7 and spin == '7/2':
        return True
    elif transition == 9 and spin == '9/2':
        return True
    else:
        return False
        

def get_splitt1(spin,transition):
    ''' Returns the split-t1 ratio as a function of the transition and the spin '''
    
    splitt1 = []
    
    if transition == 3:
        if spin == '3/2':
            splitt1.append('7/16')
            splitt1.append('9/16')
            splitt1.append('7/9')
            splitt1.append('16/7')
            splitt1.append('16/9')
            
        elif spin == '5/2':
            splitt1.append('19/31')
            splitt1.append('12/31')
            splitt1.append('19/12')
            splitt1.append('31/19')
            splitt1.append('31/12')
            
        elif spin == '7/2':
            splitt1.append('101/146')
            splitt1.append('45/146')
            splitt1.append('101/45')
            splitt1.append('146/101')
            splitt1.append('146/101')
            
        elif spin == '9/2':
            splitt1.append('91/127')
            splitt1.append('36/127')
            splitt1.append('91/36')
            splitt1.append('127/91')
            splitt1.append('127/36')
            
    if transition == 5:
        if spin == '5/2':
            splitt1.append('25/37')
            splitt1.append('12/37')
            splitt1.append('25/12')
            splitt1.append('37/25')
            splitt1.append('37/12')
            
        elif spin == '7/2':
            splitt1.append('11/20')
            splitt1.append('9/20')
            splitt1.append('19/11')
            splitt1.append('20/11')
            splitt1.append('20/9')
            
        elif spin == '9/2':
            splitt1.append('96/131')
            splitt1.append('36/131')
            splitt1.append('95/36')
            splitt1.append('131/96')
            splitt1.append('131/36')
            
    if transition == 7:
        if spin == '7/2':
            splitt1.append('161/206')
            splitt1.append('45/206')
            splitt1.append('161/45')
            splitt1.append('206/161')
            splitt1.append('206/45')
            
        elif spin == '9/2':
            splitt1.append('7/25')
            splitt1.append('18/25')
            splitt1.append('7/18')
            splitt1.append('25/7')
            splitt1.append('25/18')
            
    if transition == 9:
        splitt1.append('31/37')
        splitt1.append('6/37')
        splitt1.append('31/6')
        splitt1.append('37/31')
        splitt1.append('37/6')
    # ----------------------------------------------------------------
    return splitt1


def get_bruker_mqmas_pp_splitt1(spin,transition):
    ''' Returns the split-t1 string as a function of the transition and the spin '''
    if transition == 3:
        if spin == '3/2':
            splitt1 = '''"sqd=d0*7/16"\n"tqd=d0*9/16"\n"inf1=in0*16/9"'''
        elif spin == '5/2':
            splitt1 = '''"sqd=d0*19/31"\n"tqd=d0*12/31"\n"inf1=in0*31/12"'''
        elif spin == '7/2':
            splitt1 = '''"sqd=d0*101/146"\n"tqd=d0*45/146"\n"inf1=in0*146/45"'''
        elif spin == '9/2':
            splitt1 = '''"sqd=d0*91/127"\n"tqd=d0*36/127"\n"inf1=in0*127/36"'''
    if transition == 5:
        if spin == '5/2':
            splitt1 = '''"sqd=d0*25/37"\n"tqd=d0*12/37"\n"inf1=in0*37/12"'''
        elif spin == '7/2':
            splitt1 = '''"sqd=d0*11/20"\n"tqd=d0*9/20"\n"inf1=in0*20/9"'''
        elif spin == '9/2':
            splitt1 = '''"sqd=d0*96/131"\n"tqd=d0*36/131"\n"inf1=in0*131/36"'''
    if transition == 7:
        if spin == '7/2':
            splitt1 = '''"sqd=d0*161/206"\n"tqd=d0*45/206"\n"inf1=in0*45/206"'''
        elif spin == '9/2':
            splitt1 = '''"sqd=d0*7/25"\n"tqd=d0*18/25"\n"inf1=in0*25/18"'''
    if transition == 9:
        splitt1 = '''"sqd=d0*31/37"\n"tqd=d0*6/37"\n"inf1=in0*37/6"'''
    # ----------------------------------------------------------------
    return splitt1


def get_bruker_mqf_phasecycling(transition):
    ''' Returns string blocks for the output of Topspin MQF pulse programs.
    phase_cycling = the phase cycling for the corresponding transition.'''
    if transition == 3:
        phase_cycling = '''ph1=(12) 0 1 2 3 4 5 6 7 8 9 10 11'''
    elif transition == 5:
        phase_cycling = '''ph1=(20) 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19'''
    elif transition == 7:
        phase_cycling = '''ph1=(28) 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27'''
    elif transition == 9:
        phase_cycling = '''ph1=(36) 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35'''
    return (phase_cycling + '''\nph2=(4) 0 0 0 0\nph3=(4) 2 2 2 2\nph31=(4) 0 1 2 3''' )
    

def get_bruker_mqmas_phasecycling(transition):
    ''' Returns string blocks for the output of Topspin MQMAS pulse programs.
    phase_cycling = the phase cycling for the corresponding transition.
    '''
    # -------------------------------------------------------------------------
    if transition == 3:
        phase_cycling = '''ph1=(12) 0 1 2 3 4 5 6 7 8 9 10 11
ph2=(4) 0 0 0 0
ph3=(4) 2 2 2 2
ph10= (8) {0}*12 {1}*12 {2}*12 {3}*12 {4}*12 {5}*12 {6}*12 {7}*12
ph31=(4) 0 3 2 1 0 3 2 1 0 3 2 1
         1 0 3 2 1 0 3 2 1 0 3 2
         2 1 0 3 2 1 0 3 2 1 0 3
         3 2 1 0 3 2 1 0 3 2 1 0'''

    # -------------------------------------------------------------------------
    if transition == 5:
        phase_cycling = '''ph1=(20) 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19
ph2=(4) 0 0 0 0
ph3=(4) 2 2 2 2
ph10= (8) {0}*20 {1}*20 {2}*20 {3}*20 {4}*20 {5}*20 {6}*20 {7}*20
ph31=(4) 0 3 2 1 0 3 2 1 0 3 2 1 0 3 2 1 0 3 2 1
         1 0 3 2 1 0 3 2 1 0 3 2 1 0 3 2 1 0 3 2
         2 1 0 3 2 1 0 3 2 1 0 3 2 1 0 3 2 1 0 3
         3 2 1 0 3 2 1 0 3 2 1 0 3 2 1 0 3 2 1 0 '''

    # -------------------------------------------------------------------------  
    elif transition == 7:
        phase_cycling = '''ph1=(28) 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27
ph2=(4) 0 0 0 0');
ph3=(4) 2 2 2 2');
ph10= (8) {0}*28 {1}*28 {2}*28 {3}*28 {4}*28 {5}*28 {6}*28 {7}*28')
ph31=(4) 0 3 2 1 0 3 2 1 0 3 2 1 0 3 2 1 0 3 2 1 0 3 2 1 0 3 2 1')
         1 0 3 2 1 0 3 2 1 0 3 2 1 0 3 2 1 0 3 2 1 0 3 2 1 0 3 2')
         2 1 0 3 2 1 0 3 2 1 0 3 2 1 0 3 2 1 0 3 2 1 0 3 2 1 0 3')
         3 2 1 0 3 2 1 0 3 2 1 0 3 2 1 0 3 2 1 0 3 2 1 0 3 2 1 0')'''

    # -------------------------------------------------------------------------
    elif transition == 9:
        phase_cycling = '''ph1=(36) 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35
ph2=(4) 0 0 0 0
ph3=(4) 2 2 2 2
ph10= (8) {0}*36 {1}*36 {2}*36 {3}*36 {4}*36 {5}*36 {6}*36 {7}*36
ph31=(4) 0 3 2 1 0 3 2 1 0 3 2 1 0 3 2 1 0 3 2 1 0 3 2 1 0 3 2 1 0 3 2 1 0 3 2 1
         1 0 3 2 1 0 3 2 1 0 3 2 1 0 3 2 1 0 3 2 1 0 3 2 1 0 3 2 1 0 3 2 1 0 3 2
         2 1 0 3 2 1 0 3 2 1 0 3 2 1 0 3 2 1 0 3 2 1 0 3 2 1 0 3 2 1 0 3 2 1 0 3
         3 2 1 0 3 2 1 0 3 2 1 0 3 2 1 0 3 2 1 0 3 2 1 0 3 2 1 0 3 2 1 0 3 2 1 0'''

    return phase_cycling


def get_bruker_mqcpmgmas_phasecycling(transition):
    ''' Returns string blocks for the output of Topspin MQMAS pulse programs.
    phase_cycling = the phase cycling for the corresponding transition.'''
    # -------------------------------------------------------------------------
    if transition == 3:
        phase_cycling = '''ph0=0
ph1=(12) 0 2 4 6 8 10
ph2=(4) {0}*24 {2}*24
ph3=(4) {1}*48
ph4=(4) {1}*48  {3}*48
ph5=(4) {0}*6 {1}*6 {2}*6 {3}*6
ph6=(4) {2}*6 {3}*6 {0}*6 {1}*6
ph30=0
ph31={0 2}*3 {2 0}*3'''

    # -------------------------------------------------------------------------
    if transition == 5:
        phase_cycling = ''' '''

    # -------------------------------------------------------------------------  
    elif transition == 7:
        phase_cycling = ''' '''

    # -------------------------------------------------------------------------
    elif transition == 9:
        phase_cycling = ''' '''

    return phase_cycling

# -----------------------------------------------------------------------------
# Functions for jeol delta pulse program export

def read_deltaname_from_xml(nucleus):
    ''' Nb: element number must go first: e.g. : 'H1', 'Na23' '''
    return xmlp.parse(get_modulepath() + os.sep + 'spininfo.xml').xpath(nucleus.lower() + '/delta_fullname')[0].text

def get_jeol_phase_alternation(test):
    ''' Return '2' if test is True, 3 if not'''
    if test:
        return '0'
    else:
        return '180'
    

def get_jeol_mqf_phasecycling(transition):
    ''' Returns MQF phase cycling corresponding to the given transition.
    '''

    if transition == 3:
        phase_cycling = ('	obs_phs_MQexc		=  {0,30,60,90,120,150,180,210,240,270,300,330};\n')
        required_scan = int(12)
    
    elif transition == 5:
        phase_cycling = ('	obs_phs_MQexc		=  {0,18,36,54,72,90,108,126,144,162,180,198,216,234,252,270,288,306,324,342};\n')
        required_scan = int(20)
    
    elif transition == 7:
        phase_cycling = '	obs_phs_MQexc		=  {0,90/7,2*90/7,3*90/7,4*90/7,5*90/7,6*90/7,90,8*90/7,9*90/7,10*90/7,11*90/7,12*90/7,13*90/7,180,15*90/7,16*90/7,17*90/7,18*90/7,19*90/7,20*90/7,270,22*90/7,23*90/7,24*90/7,25*90/7,26*90/7,27*90/7};\n'
        required_scan = int(28)
        
    elif transition == 9:
        phase_cycling = '	obs_phs_MQexc		=  {0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,210,220,230,240,250,260,270,280,290,300,310,320,330,340,350};\n'
        required_scan = int(36)

    return phase_cycling, required_scan


def get_jeol_mqmas_phasecycling(transition):
    ''' Returns string for MQMAS phase cycling corresponding to the given transition.
    '''
    
    if transition == 3:
        phase_cycling = '''    obs_phs_MQexc = {0,30,60,90,120,150,180,210,240,270,300,330};
	obs_phs_CTsel = {12(0),12(45) ,12(90),12 (135),12(180),12(225), 12(270),12(315)};
	obs_phs_phase_acq = {3(0,270,180,90),3(90,0,270,180),3(180,90,0,270),3(270,180,90,0)};'''
        
        required_scan = int(96)
    
    elif transition == 5:
        phase_cycling = '''    obs_phs_MQexc = {0,18,36,54,72,90,108,126,144,162,180,198,216,234,252,270,288,306,324,342};
	obs_phs_CTsel = {20(0),20(45) ,20(90),20 (135),20(180),20(225), 20(270),20(315)};
	obs_phs_phase_acq = {5(0,270,180,90),5(90,0,270,180),splitt1 = numpy.full(3,'')5(180,90,0,270),5(270,180,90,0)};'''

        required_scan = int(160)
    
    elif transition == 7:
        phase_cycling = '''    obs_phs_MQexc =  {0,90/7,2*90/7,3*90/7,4*90/7,5*90/7,6*90/7,90,8*90/7,9*90/7,10*90/7,11*90/7,12*90/7,13*90/7,180,15*90/7,16*90/7,17*90/7,18*90/7,19*90/7,20*90/7,270,22*90/7,23*90/7,24*90/7,25*90/7,26*90/7,27*90/7};
	obs_phs_CTsel = {28(0),28(45) ,28(90),28(135),28(180),28(225),28(270),28(315)};
	obs_phs_phase_acq = {7(0,270,180,90),7(90,0,270,180),7(180,90,0,270),7(270,180,90,0)};
'''
        
        required_scan = int(224)

    elif transition == 9:
        phase_cycling = '''    obs_phs_MQexc = {0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,210,220,230,240,250,260,270,280,290,300,310,320,330,340,350};
	obs_phs_CTsel = {36(0),36(45) ,36(90),36(135),36(180),36(225),36(270),36(315)};
	obs_phs_phase_acq = {9(0,270,180,90),9(90,0,270,180),9(180,90,0,270),9(270,180,90,0)};
'''

        required_scan = int(288)

    return phase_cycling, required_scan

# ----------------------------------------------------------------------------
# ------- GUI functions
    
def save_dialog(Caption = "", Filter_tupple = (("All files (*)",))):
    ''' Open save dialog box and return filepath to save a file'''
    app = QtWidgets.QApplication.instance()
    if app is None:
        app = QtWidgets.QApplication([])
        QtCore.QTimer.singleShot(100, app.quit)
        dlg = QtWidgets.QFileDialog(caption = Caption)
        dlg.setFileMode(QtWidgets.QFileDialog.AnyFile)
        dlg.setAcceptMode(QtWidgets.QFileDialog.AcceptSave)
        dlg.setNameFilters(Filter_tupple)
        if dlg.exec_():
            filenames = dlg.selectedFiles()
            return filenames
    else:
        print(app)
        

def open_dialog(Caption = "", Filter_tupple = (("All files (*)",))):
    ''' Open save dialog box to open a file and return filepath '''
    app = QtWidgets.QApplication.instance()
    if app is None:
        app = QtWidgets.QApplication([])
        QtCore.QTimer.singleShot(100, app.quit)
        dlg = QtWidgets.QFileDialog(caption = Caption)
        dlg.setFileMode(QtWidgets.QFileDialog.AnyFile)
        dlg.setAcceptMode(QtWidgets.QFileDialog.AcceptOpen)
        dlg.setNameFilters(Filter_tupple)
        if dlg.exec_():
            filenames = dlg.selectedFiles()
            return filenames
    else:
        print(app)


def select_directory(Caption = ""):
    ''' Open save dialog box to return a directory and return filepath '''
    app = QtWidgets.QApplication.instance()
    if app is None:
        app = QtWidgets.QApplication([])
        QtCore.QTimer.singleShot(100, app.quit)
        dlg = QtWidgets.QFileDialog(caption = Caption)
        dlg.setFileMode(QtWidgets.QFileDialog.Directory)
        dlg.setAcceptMode(QtWidgets.QFileDialog.AcceptOpen)
        if dlg.exec_():
            filenames = dlg.selectedFiles()
            return filenames
    else:
        print(app)


def inputfile():
    ''' Open dialog and save input file '''
    # Get file path
    inputpath = save_dialog(Caption = 'Save input file for FAM-N optimisation',
                Filter_tupple = ("Python file (*.py)","All files(*)"))
    # Write file
    if inputpath != None:
        filepath, filename = os.path.split(inputpath[0])
        write_famn_console_input(filename, filepath)
        print('File', filename ,'successfully created in', filepath)
        # Try to open the folder
        try:
            openFolder(os.path.normpath(filepath))
        except:
            pass
    
def fig_tuple(): # eps, jpeg, jpg, pdf, pgf, png, ps, raw, rgba, svg, svgz, tif, tiff
    return ("All supported picture formats (*.pdf , *.png , *.ps , *.eps , *.svg, *.jpg, *.tiff)",
                                     "Portable Document Format (*.pdf)",
                                     "Portable Network Graphics (*.png)",
                                     "Postscript (*.ps)",
                                     "Encapsulated PostScript (*.eps)",
                                     "Scalable Vector Graphics (*.svg)",
                                     "Joint Photographic Experts Group (*.jpg)",
                                     "Tagged Image File Format (*.tiff)",
                                     "All files (*)")

# ----------------------------------------------------------------------------