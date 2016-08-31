#!/usr/bin/env python
'''
################################################################################
# An LFPy example file showing how cells can be run in parallel using MPI.
# To run using MPI with 4 cpu cores, issue in terminal
# openmpirun -np 4 python example_mpi.py
#
# The example uses mpi4py with openmpi, and do not rely on NEURON's MPI.
################################################################################


However, there is one important caveat: The NEURON extension module does not initialize MPI itself, but rather delegates this job to Python. To initialize MPI in Python, one must import a Python MPI module, such as “MPI for Python” (mpi4py) (Dalcín et al., 2008), prior to importing neuron:


from mpi4py import MPI
from neuron import h
pc = h.ParallelContext()
s = "mpi4py thinks I am %d of %d,\
 NEURON thinks I am %d of %d\n"
cw = MPI.COMM_WORLD
print s % (cw.rank, cw.size, \
           pc.id(),pc.nhost())
pc.done()


'''
import pdb
import os
#from os.path import join
import numpy as np
#import matplotlib.pyplot as plt
#from matplotlib.collections import PolyCollection
import LFPy
#import sys
#if sys.version < '3':
#    from urllib2 import urlopen
#else:    
#    from urllib.request import urlopen
#import zipfile
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
    
        
h.xopen('init.hoc')
h.xopen('simulation_run.hoc')
#dir(electrode)
electrode.cal_lfp()
#pdb.set_trace()    
        
#electrode_dic = {'LFP' : electrode.LFP }#, 'somav' : cell.somav}


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

#def cellsim(cellindex):
'''main cell- and LFP simulation procedure'''
#create extracellular electrode object
#electrode = LFPy.RecExtElectrode(electrodeParameters)
#perform NEURON simulation, results saved as attributes in cell
#cell.simulate(electrode = electrode)
    
#return dict with primary results from simulation

    


'''
def run(self):
    execute the proper simulation and collect simulation results
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
run()
'''
