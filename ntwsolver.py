from neuron.rxd import plugins
from neuron.rxd.rxdException import RxDException
import ctypes
import atexit
import numpy
from neuron import rxd, h
import neuron.rxd.initializer

# connect to NTW
_c_dll = ctypes.CDLL('./connectivity.so')
_long_ptr = ctypes.POINTER(ctypes.c_long)
_float_ptr = ctypes.POINTER(ctypes.c_float)



#QUESTIONS 

# how can I send specis info say d[] values to my c code?



def _write_connectivity_flattened_c(flat_items, flat_offsets, stateList, numberOfSpecies, diffusionConstantList):
    """wrapper for the C function"""
    # The 0 at the end indicates that we are NOT using hardcoded reaction dynamics
    #print '_c_dll.write_connectivity(%r,\narray(%r),\narray(%r),\narray(%r), \n%r, \narray(%r), \n%r)' % (int(len(flat_offsets) - 1), list(flat_items), list(flat_offsets), list(stateList), numberOfSpecies, list(diffusionConstantList), 0)
    _c_dll.write_connectivity(
        int(len(flat_offsets) - 1),
        flat_items.ctypes.data_as(_long_ptr),
        flat_offsets.ctypes.data_as(_long_ptr),
        stateList.ctypes.data_as(_float_ptr),
        numberOfSpecies,
        diffusionConstantList.ctypes.data_as(_long_ptr),
        0)

def _runMainLoop_c(nextTimeStep):
    """wrapper for the C function"""
    #print '_c_dll.runMainLoop(%r)' %  nextTimeStep
    _c_dll.runMainLoop(nextTimeStep.ctypes.data_as(_float_ptr))

def _get_connectivity(region):
    indices = _get_indices_dict(region)
    xs, ys, zs = region._xs, region._ys, region._zs
    return [_get_neighbor_indices(xs[i], ys[i], zs[i], indices) for i in xrange(len(xs))]

def _get_indices_dict(region):
    indices = {}
    xs, ys, zs = region._xs, region._ys, region._zs
    for i in xrange(len(xs)):
        indices[(xs[i], ys[i], zs[i])] = i
    return indices

def _get_neighbor_indices(i, j, k, indices):
    result = []
    for pt in [(i - 1, j, k), (i + 1, j, k), (i, j - 1, k), (i, j + 1, k), (i, j, k - 1), (i, j, k + 1)]:
        lookup_value = indices.get(pt)
        if lookup_value is not None:
            result.append(lookup_value)
    return result


def _flatten(list_of_lists):
    """turns a list of lists into 2 numpy arrays of integers
    
    The first return value is the data, the second is the first index of each
    list. The second array is one longer than the first to indicate the le
    """
    import itertools
    data = numpy.array(list(itertools.chain(*list_of_lists)), dtype=int)
    start_indices = [0]
    for item in list_of_lists:
        start_indices.append(len(item) + start_indices[-1])
    start_indices = numpy.array(start_indices, dtype=int)
    return data, start_indices

def _initSVs_c(dx):
    """wrapper for the C function"""
    #print '_c_dll.initSVs(%r)' % dx
    _c_dll.initSVs(dx.ctypes.data_as(_float_ptr))

# when we're done, cleanup
atexit.register(_c_dll.finalization)

class NTWSolver(plugins.SharedMemorySolverPlugin):
    def _set_mesh(self):
        # TODO: check to make sure that all of the species still exist

        #print 'hello from _set_mesh', len(self._species)
        
        if hasattr(self, '_meshsetup'):
            # TODO: this should not actually be a problem; just throw out the old solver
            raise RxDException('mesh already configured')
        
        if len(self._species):
            self._meshsetup = True
            
            # TODO: note that a species could no longer exist
            # TODO: allow more than one region
            # TODO: note this doesn't make any sense in a mosaic sim where some parts
            #       are 1D and others are 3D
            geom_region = self._species[0]._regions[0]
            if not hasattr(geom_region, '_xs'):
                geom_region._do_init()
            
            # Geometry's neighbor information that is list of neighbor lists are created . . . .
            connectivity = _get_connectivity(geom_region) 
            #data, dataIndeces getting from here . . .
            flat_items, flat_offsets = _flatten(connectivity)

            # sending geometry and specis informations to stochastic simulator
            diffusion_constant_list = numpy.array([float(s._d) for s in self._species])
            _write_connectivity_flattened_c(flat_items, flat_offsets, rxd.node._states, len(self._species), diffusion_constant_list)

            # Subvolume's dx, dy, and dz are sending to NTW from here. here dx=dy=dz
            # TODO: ultimately this is not where dx SHOULD be read from, but for now it's the
            #       only place that knows
            _initSVs_c(numpy.array([geom_region.dx]))

    def _process_reactions(self):
        #print 'hello from _process_reactions', len(self._reactions)
        pass
    
    def advance(self, dt):
        print 'advance: h.t = %g' % h.t
        time_array = numpy.array([h.t + dt])
        _runMainLoop_c(time_array)
        print 'advanced to: %g' % (h.t + dt)
        for s in self._species:
            if s is not None: s._transfer_to_legacy()
    
    def supported_sim(self):
        """TODO: something better"""
        return True

    def update_value(self, node):
        """For a sim that did not use shared memory, this would be where we would change things if a node changed"""
        pass

rxd.set_solve_type(dimension=3)
plugins.set_solver(NTWSolver())
