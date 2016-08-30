import mpi4py
from mpi4py import MPI

#MPI stuff we're using
COMM = MPI.COMM_WORLD
SIZE = COMM.Get_size()
RANK = COMM.Get_rank()



from neuron import h
pc = h.ParallelContext()

id = int(pc.id())
nhost = int(pc.nhost())
nhost = int(SIZE)

print 'I am %d of %d'%(id, nhost)
print 'I am %d of %d'%(RANK, nhost)

pc.runworker()
pc.done()
h.quit()
