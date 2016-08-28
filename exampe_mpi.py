#!/usr/bin/env python
'''

################################################################################
# NB. Russell.
# THis is a fushion of code from LFPy and Traub.
# Traub hoc code is called somewhere about half way down, its sourced from a different file (init.hoc)




# An LFPy example file showing how cells can be run in parallel using MPI.
# To run using MPI with 4 cpu cores, issue in terminal
# openmpirun -np 4 python example_mpi.py
#
# The example uses mpi4py with openmpi, and do not rely on NEURON's MPI.
################################################################################
'''

import os
from os.path import join
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
import LFPy
import sys
#if sys.version < '3':
#    from urllib2 import urlopen
#else:    
#    from urllib.request import urlopen
import zipfile
from mpi4py import MPI

#MPI stuff we're using
COMM = MPI.COMM_WORLD
SIZE = COMM.Get_size()
RANK = COMM.Get_rank()


if RANK == 0:
    #compile mod files
    os.system('''
              cd cells
              nrnivmodl
              ''')

    #sync threads
COMM.Barrier()


LFPy.cell.neuron.load_mechanisms('cells')    


from neuron import h

# Call the initialisation of Traubs model, using HOC native NEURON code the contents of the
# h object will reflect NEURON objects in the Python name space.

h.xopen('init.hoc')

# Try dir(h) to check out the hoc variables existing in the python name space.

#LFPy.cell.neuron.load_mechanisms('cells')




# Create a grid of measurement locations, in (mum)
X, Z = np.mgrid[-500:501:20, -400:1201:40]
Y = np.zeros(X.shape)

# Define electrode parameters
electrode_parameters = {
    'sigma' : 0.3,      # extracellular conductivity
    'x' : X.flatten(),  # electrode requires 1d vector of positions
    'y' : Y.flatten(),
    'z' : Z.flatten()
}

# Create electrode object
electrode = LFPy.RecExtElectrode(**electrode_parameters)
    
        
    
        


'''
def distribute_cellsims(self):
	#start unique cell simulation on every RANK,
	#and store the electrode and cell objects in dicts indexed by cellindex
	results = {}
	for cellindex in range(self.POPULATION_SIZE):
	    if divmod(cellindex, SIZE)[1] == RANK:
		results.update({cellindex : self.cellsim(cellindex)})
	return results
'''

def cellsim(cellindex):
    '''main cell- and LFP simulation procedure'''
    #create extracellular electrode object
    electrode = LFPy.RecExtElectrode(**self.electrodeParameters)
    #perform NEURON simulation, results saved as attributes in cell
    cell.simulate(electrode = electrode)
    
    #return dict with primary results from simulation
    return {'LFP' : electrode.LFP, 'somav' : cell.somav}

    



def run(self):
    '''execute the proper simulation and collect simulation results'''
    #produce simulation results on each RANK
    self.results = self.distribute_cellsims()
    
    #superimpose local LFPs on every RANK, then sum using MPI to RANK 0
    self.LFP = []
    for key, value in list(self.results.items()):
        self.LFP.append(value['LFP'])
    self.LFP = np.array(self.LFP).sum(axis=0)
    self.LFP = COMM.reduce(self.LFP)        #LFP is None on all but RANK 0
    
    #collect all simulation results on RANK 0, including single cell LFP
    if RANK == 0:
        for i in range(1, SIZE):
            result = COMM.recv(source=MPI.ANY_SOURCE) #receive from ANY rank
            self.results.update(result)     #collect
    else:
        COMM.send(self.results, dest=0)     #send to RANK 0
        self.results = None                 #results only exist on RANK 0
        
    COMM.Barrier()  #sync MPI threads
