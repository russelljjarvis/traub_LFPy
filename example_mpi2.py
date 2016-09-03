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
from mpi4py import MPI

pc = h.ParallelContext()

#import LFPy

#MPI stuff we're using
COMM = MPI.COMM_WORLD#h.ParallelContext()#MPI.COMM_WORLD
SIZE = COMM.Get_size()#int(pc.nhost)#
RANK = COMM.Get_rank()#int(pc.id)#


if RANK == 0:
    #compile mod files
    os.system('''
              cd cells
              nrnivmodl
              ''')

    #sync threads
COMM.Barrier()
h.xopen('init.hoc')
h.xopen('simulation_run.hoc')
h.define_shape()
cells=h.cells
ncell=len(cells)
celldict={}#construct a dictionary. Keys: Global IDentifiers. Values: cells.
# This way can use syntax. If GID in celldict    
for i in xrange(0,SIZE*ncell):
	if pc.gid_exists(i):
		celldict[i]=pc.gid2cell(i) 
		
#for k, v in celldict.iteritems():
#	print k,v

#seglist= iter( (seg, sec, cell) for sec in cell for seg in sec )     





#itercell= ( (i,t) for i,t in self.celldict.iteritems() if i in self.celldict.keys() if int(t.gid1) != int(k['gid']) )       
#suspect lambda rule needs to be called before x is a thing.
hoc_string='''
forall {
   //for (x, 0) {
      insert xtra
   //}
}
'''
h(hoc_string)
h('forall{ for(x,0){ insert xtra }}')
h('forall{ for(x,0){ insert extracellular}}')    
h('xopen("interpxyz.hoc")')
h('grindaway()')  
pdb.set_trace()
print 'got here?'
coordinates=[]
for cell in celldict.values():
    #seglist= iter( (seg, sec, self.celldict[j]) for sec in self.celldict[j].spk_trig_ls for seg in sec )     
    #for (seg,sec, cellc) in seglist:	   
	seglist= iter((seg, sec) for sec in cell.all for seg in sec)
	for seg, sec in seglist:
		sec.push()
		h('print x_xtra()')#, y_xtra(), z_xtra()')
		h('print psection()')
		'''h('objref coords')
		h('coords = new Vector()')
		get_cox = str('coords.x[0]=x_xtra('
		              + str(0.5) + ')')
		h(get_cox)                   
		get_coy = str('coords.x[1]=y_xtra('
		              + str(0.5) + ')')
		h(get_coy)
		get_coz = str('coords.x[2]=z_xtra('
		              + str(0.5) + ')')
		h(get_coz)
		coordinates.append(h.coords.to_python())
		'''
		h.pop_section()

seglist= [ seg for cell in cells for seclist in cell.all for sec in seclist for seg in sec ]
#for seg in seglist:
#	print seg
	#h.insert('xtra')

h('forall{ insert xtra }')
h('forall{ insert extracellular }')

#h('forall{ { psection() } }')
#seglist= iter( (seg, sec, cell) for sec in cell for seg in sec ) 
seglist= iter( (seg, sec, cell) for cell in cells for sec in cell for seg in sec ) 

#for cell in cells:
#	for seg in cell.Soma.Section():
#		print seg.x

#h('forall{ for(0,x) { insert xtra  } }')
print 'got here'
h('grindaway()') 

h('xopen("interpxyz.hoc")')

#POINTER im, ex2, ex

h('forall{ if (ismembrane("xtra")){ for(x,0) { setpointer im_xtra(x), i_membrane(x) } } }')    
h('forall{ if (ismembrane("xtra")){ for(x,0) { ex_xtra(x), e_extracellular(x) } } }')    

exec_hoc='''
//https://www.neuron.yale.edu/phpBB/viewtopic.php?f=31&t=3083
proc fieldrec() { local sum
  sum = 0
  forall {
    for(x,0) {

		if (ismembrane("xtra")) {
		   sum += er_xtra(x)
		}
	}
  }
  //return sum
}
'''
h(exec_hoc)
local_sum=h.fieldrec()
print local_sum

#############################################################
########    Simulation control    ##########################
{pc.set_maxstep(10)}
h.stdinit()
{pc.psolve(h.tstop)}
###############################################################
        
def prun(self,tstop):
    h=self.h    
    pc=h.ParallelContext()
    NCELL=self.NCELL
    SIZE=self.SIZE
    COMM = self.COMM
    RANK=self.RANK
    checkpoint_interval = 50000.

    #The following definition body is from the open source code at:
    #http://senselab.med.yale.edu/ModelDB/ShowModel.asp?model=151681
    #with some minor modifications
    cvode = h.CVode()
    cvode.cache_efficient(1)
    # pc.spike_compress(0,0,1)

    pc.setup_transfer()
    mindelay = pc.set_maxstep(10)
    if RANK == 0:
        print 'mindelay = %g' % mindelay
    runtime = h.startsw()
    exchtime = pc.wait_time()

    inittime = h.startsw()
    h.stdinit()
    inittime = h.startsw() - inittime
    if RANK == 0:
        print 'init time = %g' % inittime

    while h.t < tstop:
        told = h.t
        tnext = h.t + checkpoint_interval
        if tnext > tstop:
            tnext = tstop
        pc.psolve(tnext)
        if h.t == told:
            if RANK == 0:
                print 'psolve did not advance time from t=%.20g to tnext=%.20g\n' \
                    % (h.t, tnext)
            break    
        print 'working', h.t
    runtime = h.startsw() - runtime
    comptime = pc.step_time()
    splittime = pc.vtransfer_time(1)
    gaptime = pc.vtransfer_time()
    exchtime = pc.wait_time() - exchtime
    if RANK == 0:
        print 'runtime = %g' % runtime
    print comptime, exchtime, splittime, gaptime

# Code above should be translated into the file below.

h.xopen('simulation_run.hoc')

global_sum=pc.allreduce(local_sum, 1)
print global_sum
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
	
'''call_hoc_f=
func lambda_f() { local i, x1, x2, d1, d2, lam
        if (n3d() < 2) {
                return 1e5*sqrt(diam/(4*PI*$1*Ra*cm))
        }
// above was too inaccurate with large variation in 3d diameter
// so now we use all 3-d points to get a better approximate lambda
        x1 = arc3d(0)
        d1 = diam3d(0)
        lam = 0
        for i=1, n3d()-1 {
                x2 = arc3d(i)
                d2 = diam3d(i)
                lam += (x2 - x1)/sqrt(d1 + d2)
                x1 = x2   d1 = d2
        }
        //  length of the section in units of lambda
        lam *= sqrt(2) * 1e-5*sqrt(4*PI*$1*Ra*cm)

        return L/lam
}

proc geom_nseg() {
  //soma area(.5) // make sure diam reflects 3d points
  forsec all { nseg = int((L/(0.1*lambda_f(100))+.9)/2)*2 + 1  }
}
'''
#h(call_hoc_f)


