# To build NEURON+PYTHON+MPI
Try this experimental docker file at
https://github.com/russelljjarvis/PyNeuron-Toolbox

# To test if the docker container works:
git clone https://github.com/russelljjarvis/traub_LFPy.git
cd traub_LFPy
and run test0.py

Then try
$mpiexec -np 4 ipython -i example_mpi.py

# traub_LFPy
There are two approaches to recording the LFP in Traub. The first approach is to use LFPy. The second approach is to use the mod files: .xtra.mod and extracellular.mod, and to construct a grid of point processors corresponding to recording electrodes. The HOC/MOD way is usually more complex/tedious and less portable/reproducible so the LFPy way is probably more warranted.

This approach might be better, the nice thing about LFPy


## 0. 
Inside the Traub model
Traubs model uses NEURONs MPI communicator object, however LFPy uses mpi4py. There are some NEURON examples that show how you can make mpi4py.COMM a proxy for NEURONs MPI communicator, such that we can traubs NEURON mpi communication and mpi4py can coexist, and also objects from the NEURON name space can be proxied by mpi4py instead, when we do other MPI communication operations later.

##1. 
Move all the mod files from their own special directory to the trunk Neuron directory (up a level of tree).

##2. 
compile there.

##3. 
Make an init.py that simply imports the neuron module and then uses that to open init.hoc. Does this drag all of the HOC objects into the python space? Yes. Good.

##4. 
Comment out the part of init.hoc that runs the traub simulation, as our code has to initialise LFP recording electrodes in a mesh properly fir st.

##5. 
Launch the script with $ ipython -i init.py # an interactive mode.

##7. 
$ run dir(h) # To check what is in the Python name space. Everything is there.

##9. 
ipython -c "import LFPy"
NEURON -- Release 7.4 (1370:16a7055d4a86) 2015-11-09
Duke, Yale, and the BlueBrain Project -- Copyright 1984-2015
See http://www.neuron.yale.edu/neuron/credits


Need LFP-mpi initialised code before Traub simulation is started, LFP is recorded on dots of grid line intersections on a sqaure surface (the LFPs recorded at each point can all be summed togethor later).

# Create a grid of measurement locations, in (mum)
Currently this grid only makes sense relative to the LFPY-MPI example, however we want it to make be a sensible surface location relative to the traub model.

If we are trying to match up our data with experimental eCOG intracranial recordings from epilepsy. Then our surface location is probably the top surface of Traubs cortex's outer most layer (if it is layered).

X, Z = np.mgrid[-500:501:20, -400:1201:40]
Y = np.zeros(X.shape)

# Define electrode parameters
electrode_parameters = {
    'sigma' : 0.3,      # extracellular conductivity
    'x' : X.flatten(),  # electrode requires 1d vector of positions
    'y' : Y.flatten(),
    'z' : Z.flatten()
}

electrode = LFPy.RecExtElectrode(**electrode_parameters)


##11. 
## First significant problem. 

LFPs electrode results are temporarily stored inside LFPy.cell objects (why?).

THe constructor for an LFPy.cell object requires a morphology file (swc). Traubs model defines morphologies by HOC code not swc. This is a bummer.

Need to find out if its really necessary to store LFP vectors inside cell objects, or if they can be stored as a dictionary for each host instead.


def cellsim(self, cellindex):
    '''main cell- and LFP simulation procedure'''
    #create extracellular electrode object
    electrode = LFPy.RecExtElectrode(**self.electrodeParameters)
    
    #Initialize cell instance, using the LFPy.Cell class
    cell = LFPy.Cell(**self.cellParameters)
    #set the position of midpoint in soma
    cell.set_pos(xpos = self.cellPositions[cellindex, 0],
                 ypos = self.cellPositions[cellindex, 1],
                 zpos = self.cellPositions[cellindex, 2])
    #rotate the morphology
    cell.set_rotation(z = self.cellRotations[cellindex])
    
    #attach synapse with parameters and set spike time
    synapse = LFPy.Synapse(cell, **self.synapseParameters)
    synapse.set_spike_times(self.synapseTimes[cellindex])
    
    #perform NEURON simulation, results saved as attributes in cell
    cell.simulate(electrode = electrode)
    
    #return dict with primary results from simulation
    return {'LFP' : electrode.LFP, 'somav' : cell.somav}


The LFPy project requires LFPy.cell class 
Help on module LFPy.cell in LFPy:

NAME
    LFPy.cell - Copyright (C) 2012 Computational Neuroscience Group, UMB.

FILE
    /usr/local/lib/python2.7/dist-packages/LFPy/cell.py

DESCRIPTION
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

CLASSES
    __builtin__.object
        Cell
    
    class Cell(__builtin__.object)
     |  The main cell class used in LFPy.
     |  
     |  Arguments:
     |  ::
     |      
     |      morphology : [str]: path/to/morphology/file
     |  
     |      v_init: [-65.]: initial potential
     |      passive: [True]/False: passive mechs are initialized if True
     |      Ra: [150.]: axial resistance
     |      rm: [30000]: membrane resistivity
     |      cm: [1.0]: membrane capacitance
     |      e_pas: [-65.]: passive mechanism reversal potential
     |      extracellular: [True]/False: switch for NEURON's extracellular mechanism
     |  
     |      timeres_NEURON: [0.1]: internal dt for NEURON simulation
     |      timeres_python: [0.1]: overall dt for python simulation
     |  
     |      tstartms: [0.]:  initialization time for simulation <= 0 ms
     |      tstopms: [100.]: stop time for simulation > 0 ms
     |  
     |      nsegs_method: ['lambda100']/'lambda_f'/'fixed_length': nseg rule
     |      max_nsegs_length: [None]: max segment length for method 'fixed_length'
     |      lambda_f: [100]: AC frequency for method 'lambda_f'
     |      d_lambda: [0.1]: parameter for d_lambda rule
     |      
:
