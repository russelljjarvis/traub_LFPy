# To build NEURON+PYTHON+MPI
Try this experimental docker file at
https://github.com/russelljjarvis/PyNeuron-Toolbox

# To test if the docker container works:
git clone https://github.com/russelljjarvis/traub_LFPy.git
cd traub_LFPy
and run test0.py

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
NEURON version of Traub et al J Neurophysiol. 2005 Apr;93(4):1829-30.
Single-column thalamocortical network model exhibiting gamma oscillations
sleep spindles, and epileptogenic bursts.

See: http://senselab.med.yale.edu/senselab/modeldb/ShowModel.asp?model=45539

Prepare for running with
nrnivmod mod

See the comparison of each cell type with the fortran output with
nrngui onecell.hoc

Selecting the "Exact traub style ri" which forces the NEURON connection
coefficients to be exactly the same as computed by the Traub's equivalent
circuit style from the FORTRAN demonstrates that cell and channel properties
are reliably represented in this NEURON translation of the FORTRAN code.


Reliability:

Since the NEURON and Fortran models do not produce quantitatively
identical results there is always some question as to whether simulation
differences are due to substantive parameter translation errors or
can be attributed to different numerical methods.
It must be realized that our experience has
been that every test into a new runtime domain has exhibited
discprepancies that were ultimately resolved by fixing
translation errors.

And also that our comparison
tests are only between NEURON and an already significantly modified
FORTRAN code. The bulk of the FORTRAN modifications are toward more
generic FORTRAN syntax to allow the original ifc compatible FORTRAN to
run under g77. Most of the execution places where the g77 version
differs from the ifc version are straight forward transformations of
bulk array assignment into equivalent elementwise assignment via do loops.
Did we get them all? Did we assign over ALL the elements in each array?
We did manually review all ifc to g77 editing changes but
a few cases involved our judgement with regard to whether there was a bug
in the original ifc fortran version.  The modified FORTRAN used for the
NEURON comparisons is available from this model=45539 page. As is, the
g77 FORTRAN model can only be run as 14 processes, one for each cell type
and a full model run takes 20 hours or so. Simplifying to 1/10 the number
of cells gives a model that takes approximate 1.5 hours for 100 ms of
simulation time. Our last network bug, based on significant spike raster
plot discrepancies, was found using a 10 ms run.

We consider the translation of the 14 individual cell types to be quite
reliable based on the quantitative similarity of the g77 and NEURON
isolated 100 ms cell trajectories at the spike detection location for
0 and 0.3 nA constant current stimulation into the soma. Note that
quantitative similarity demands compartment coupling of exactly the
same values used by the FORTRAN version algorithms (imitated in NEURON
using the "traub_exact()" algorithm where some branch points had the
form of "wye" equivalent circuits, some had the "delta" form, and all had
a different view of how resistance from child to parent should be computed.)

Network topology and chemical synapse parameter reliability is limited to the
diagnostic power of our specific tests. For quantitative comparisons we
printed the precise FORTRAN network topology to files and used that information
to define the NEURON network connections. For 10 ms with a 1/10 size network
we focused on quantitative similarity of the spike raster plots. The FORTRAN
version has a spike resolution time of 0.1 ms and all synaptic conductance
trajectories are step functions with that resolution (the underying
dt is 50 times smaller, dt = 0.002). We prepared a special version of the
NEURON executable to force spike threshold detection on 0.1 ms boundaries
to allow convenient comparison of spike rasters. For the first 10 ms
we judged whether spike discrepancies were due to FORTRAN-NEURON spikes
straddling the 0.1 ms boundaries or whether the discrepancy was likely to
be due to a topology or synaptic parameter error. The judgement was based
on the details of the voltage trajectory at the spike detector compartment.
We believe that careful analysis of the first
10 ms of the  100 ms spike raster overlap plot for the
FORTRAN (fat red marks) and NEURON (thin black marks) in combination with
the spike trajectory sensitivity of suppyrRS cells with respect to number
of spikes in their burst after the first spike will convince that we
have gone as far as possible with quantitative spike location similarity
as a diagnostic technique. Further diagnostics will likely have to be
based on specific questions in regard to qualitative discrepancies and
a focus on the NEURON model itself as the tool for exploration in terms
of certain properties added or subtracted from successive runs. Unfortunately
that can only be done in response to a specific suspicion on the part of the
user. Enumerated below are the known discrepancies between the representations
of the Traub model in FORTRAN and NEURON:

The NEURON nmda saturation is turned off. See the NMDA_saturation_fact in
the FORTRAN groucho.f file and the nrntraub/mod/traub_nmda.mod file.
Warning: NEURON will not mimic the FORTRAN
merely by setting the factor to 80.

The groucho.f axon_refrac_time is normally set to 1.5. Our quantitative
tests temporarily set this parameter to 0.5. The NEURON spike detection
algorithm defines a spike as a positive going transition past the
trigger value.

Enumerated below are those major components (of which we
are aware) that are in the model but have not been tested in terms of their
quantitative equivalence to FORTRAN:
gap junctions.
long term nmda properties and effects.
ectopic spikes
random current stimulation

The bottom line: The spike rasters for a
full g77 FORTRAN run (gap junctions and current stimulation present) and
a full NEURON run with its own independent random variables
(no "traub_exact" connection coefficients, dt = 0.025 (ms),
secondorder = 2, and spike detection with dt resolution) is presented.
