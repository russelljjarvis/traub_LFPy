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
'''





import pdb
import os
import numpy as np
from neuron import h
pc = h.ParallelContext()

#import LFPy
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
h.xopen('init.hoc')

#h('forall{ { psection() } }')
h('forall{ { insert xtra } }')
h('grindaway()') 
h('forall{ { insert extracellular} }')
h('xopen("interpxyz.hoc")')

#POINTER im, ex2, ex

h('forall{ if (ismembrane("xtra")){ { setpointer im_xtra(x), i_membrane(x) } } }')    
h('forall{ if (ismembrane("xtra")){ { ex_xtra(x), e_extracellular(x) } } }')    

exec_hoc='''
//https://www.neuron.yale.edu/phpBB/viewtopic.php?f=31&t=3083
proc fieldrec() { local sum
  sum = 0
  forall {
    if (ismembrane("xtra")) {
       sum += er_xtra(x)
    }
  }
  //return sum
}
'''
h(exec_hoc)
local_sum=h.fieldrec()
print local_sum
h.xopen('simulation_run.hoc')



#LFPy.cell.neuron.load_mechanisms('cells')    


#from neuron import h

# Call the initialisation of Traubs model, using HOC native NEURON code the contents of the
# h object will reflect NEURON objects in the Python name space.


# Try dir(h) to check out the hoc variables existing in the python name space.

#LFPy.cell.neuron.load_mechanisms('cells')



'''
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
grid_electrode = LFPy.RecExtElectrode(**electrode_parameters)
    

# Define electrode parameters
point_electrode_parameters = {
    'sigma' : 0.3,  # extracellular conductivity
    'x' : np.array([-130., -220.]),
    'y' : np.array([   0.,    0.]),
    'z' : np.array([   0.,  700.]),
}
        

'''
'''
For email:

Above, the morphology and specifications of the biophysical properties were given as keyword arguments to LFPy.Cell. Models existing in memory can in principle be executed by supplying the keyword arguments 
morphology=None, delete_sections=False in addition to the above cell_parameters, e.g., if a model is defined via the NEURON graphical or command line interface. Defining scripts with the full specification of the model loaded with LFPy.Cell, is, however, in most cases more tractable.


While the present version 1.0 of LFPy is focused on the calculation of single-neuron contributions to the extracellular potentials, the computational scheme generalizes directly to the calculation of signals from populations of neurons. This was illustrated in Example 3 in section 3.3 where also parallelization of the computational scheme by means of MPI was employed, however without communication between units. At present, LFPy is however less suitable for the investigation of extracellular potentials generated in genuine network models that require parallelization of the network activity. At present, this is a limitation in the current version of the software mainly in that the simulation control is incorporated as an LFPy.Cell class method, and that the class LFPy.TemplateCell (which allows for multiple simultaneous cell representations) is not using the capabilities of NEURON for assigning each cell to different MPI ranks. However, as simulation of extracellular signals from network activity likely will become increasingly important, we aim to implement solutions to these limitations in future versions of LFPy.

#dir(LFPy.Cell._run_custom_codes)
'''

'''
#pdb.set_trace()
cells_py=h.cells
for c in cells_py:
	h.psection() in c.Soma.Section() 
	
h.psection() in cells_py[1].Soma.Section() 

for c in cells_py:
     print c     
     c.Soma.Section().push()
     #h.secname()
     #h.psection()
     h.cas()
     h.pop_section()



# Define cell parameters
for c in cells_py:
	Need to get rm,cm e_pas, and v_init out of existing cell object hoc objects stored in memory. How I don't know.
	print c.Soma.Section().Ra
	print c.Soma.Section().L
	

	cell_parameters = {          # various cell parameters,
		'morphology' : None, 
		'delete_sections' : False,
		'rm' : 30000.,      # membrane resistance
		'cm' : 1.0,         # membrane capacitance
		'Ra' : c.Soma.Section().Ra,        # axial resistance
		'v_init' : -65.,    # initial crossmembrane potential
		'e_pas' : -65.,     # reversal potential passive mechs
		'passive' : True,   # switch on passive mechs
		'nsegs_method' : 'lambda_f',
		'lambda_f' : 100.,
		'timeres_NEURON' : 2.**-3,   # [ms] dt's should be in powers of 2 for both,
		'timeres_python' : 2.**-3,   # need binary representation
		'tstartms' : 0.,    # start time of simulation, recorders start at t=0
		'tstopms' : 100.,   # stop simulation at 200 ms. These can be overridden
		                    # by setting these arguments i cell.simulation()
	}

	# Create cell
	cell = LFPy.Cell(**cell_parameters)
'''	
	

