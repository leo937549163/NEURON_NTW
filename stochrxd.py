import neuron.rxd as rxd
from neuron.rxd.rxdException import RxDException
from neuron import h
from neuron.rxd.node import _molecule_node
from neuron.rxd.species import _get_all_species
import subprocess
import numpy

# TODO: this currently assumes a spatially uniform diffusion constant; fix this

_the_region = None

def Region(sections, nrn_region=None, geometry=None, dimension=1, dx=None):
    global _the_region
    if dimension != 3:
        raise RxDException('stochrxd only supports 3D so far')
    if _the_region is not None:
        raise RxDException('stochrxd only supports 1 region so far')
    _the_region = rxd.Region(sections, nrn_region=nrn_region, geometry=geometry, dimension=dimension, dx=dx)
    return _the_region

def Reaction(scheme, rate1, rate2=None, regions=None, custom_dynamics=False):
    # the following check is because of a limitation in the stochastic solver
    if custom_dynamics != False:
        raise RxDException('must leave custom_dynamics == False; only mass action supported right now')
    return rxd.Reaction(scheme, rate1, rate2=rate2, regions=regions, custom_dynamics=custom_dynamics)

def Species(regions=None, d=0, name=None, charge=0, initial=None):
    """unlike Species in neuron.rxd, initial here refers to a # of molecules"""
    s = rxd.Species(regions=regions, d=d, name=name, charge=charge, initial=initial)
    # TODO: this is part of the part that assumes 3D (because for 1D node list constructed only on demand)
    for node in s._nodes:
        node._data_type = _molecule_node
    return s
    
def continuerun(tstop):
    if tstop <= h.t: return
    all_species = [s() for s in _get_all_species().values() if s() is not None]
    # NOTE: the StochasticSolver is responsible for updating state values
    s = StochasticSolver(all_species, rxd.rxd._all_reactions)    
    s.advance(tstop - h.t)
    h.t = tstop
    
class StochasticSolver:
    def __init__(self, species, reactions):
        # TODO: work with weak references instead of strong references
        self._species = species
        self._reactions = reactions
        
        # write the molecule counts
        self._write_molecule_counts_to_file()
        
        print 'connecting to SLV...'
        self._p = subprocess.Popen(['./SLV'], stdout=subprocess.PIPE, stdin=subprocess.PIPE)
        
        assert(len(species) > 0)

        #
        # send the geometry
        #
        
        indices = _get_indices_dict(species[0])

        line = 'None'
        print 'waiting to be prompted for number of cubes'
        self._advance_past('--------------------------')
        # TODO: probably don't need the .keys() bit; does len(indices) suffice?
        self._write('%d\n' % len(indices.keys()))
        
        self._advance_past('dx:')
        self._write('%g\n' % _the_region._dx)
        for i, node in enumerate(species[0]._nodes):
            print 'processing geometry node %d of %d' % (i + 1, len(species[0]._nodes))
            print 'waiting to be asked for location'
            self._advance_past('location x y z')
            # TODO: this is a place where I'm assuming 3D
            self._write('%g\n%g\n%g\n' % (node._i, node._j, node._k))
            print 'waiting to be asked about neighbors'
            self._advance_past('num neighbors')
            neighbors = _get_neighbor_indices(node, indices)
            self._write('%d\n' % len(neighbors))
            for neighbor in neighbors:
                self._write('%d\n' % neighbor)
        
        #
        # send the species information
        #
        
        self._advance_past("How many species in your domain?")
        self._write('%d\n' % len(species))
        for s in species:
            self._advance_past('Enter Name of species')
            self._write('%s\n' % s._real_name)
            self._advance_past('Enter diffusion Constant for species')
            self._write('%g\n' % s._d)
        
        #
        # TODO: send the reactions
        #
        
        self._advance_past('How many reactions?')
        if len(reactions) > 0:
            raise RxDException('reactions not yet supported')
        else:
            self._write('0\n')

        self._prepare_to_send_time()        
        
            
            

    def _write_molecule_counts_to_file(self):
        all_nodes = [s._nodes for s in self._species]
        with open('data.in', 'w') as f:
            for row in zip(*all_nodes):
                row = ' '.join(str(node.molecules) for node in row)
                f.write(row + '\n')

    def advance(self, t):
        """
            .. warning::
            
                don't do two concurrent advances because could have a
                collision on the state writing phase
        """
        self._write('{0}\n'.format(t))
        self._prepare_to_send_time()
        
        # break the ability to continue, just so we can guarantee it shuts down
        # TODO: I think I'm doing this because of a bug in _read_states... fix/remove
        self._write('-1\n')
        
        # now wait (remove this and the line above)
        import time
        time.sleep(0.1)
        
        self._update_molecule_counts()
    
    def _update_molecule_counts(self):
        all_data = self._read_states()
        for s, data in zip(self._species, all_data):
            for node, value in zip(s._nodes, data):
                node.molecules = value

    def _advance_past(self, data, show_output=False):
        line = 'None'
        while data not in line:
            line = self._p.stdout.readline()
            if show_output or True:
                print 'read: ', line
            if not line:
                import sys
                sys.exit()
        print 'exited while at: ', line

    def _prepare_to_send_time(self):
        self._advance_past('Enter next time . . .(-1 for exit )', True)

    def _write(self, data):
        print 'writing: %r' % data
        #print dir(self._p)
        self._p.stdin.write(data)

    def _read_states(self):
        # data.out contains the full time course but NEURON only automatically
        # stores the current values; this replicates the behavior where we only
        # can get the current values
        with open('data.out') as f:
            print 'begin read states'
            for i in xrange(3):
                f.readline()
            result = None
            while f.readline():
                print 'starting a dataset'
                if len(self._species) != 1:
                    raise Exception('reading states currently only supports 1 species')
                result = [_get_state_vec(f, s) for s in self._species]
                f.readline()
            return result

        
def _get_state_vec(f, species):
    # TODO: this will need changed when more than one species
    return [float(f.readline()) for s in species._nodes]
#    return [float(v) for v in f.readline().split(':')[1].strip().split(' ')]



def _get_indices_dict(species):
    region = species._regions[0]
    indices = {}
    xs, ys, zs = region._xs, region._ys, region._zs
    for i in xrange(len(xs)):
        indices[(xs[i], ys[i], zs[i])] = i
    return indices
    
def _get_neighbor_indices(node, indices):
    i, j, k = node._i, node._j, node._k
    result = []
    for pt in [(i - 1, j, k), (i + 1, j, k), (i, j - 1, k), (i, j + 1, k), (i, j, k - 1), (i, j, k + 1)]:
        lookup_value = indices.get(pt)
        if lookup_value is not None:
            result.append(lookup_value)
    return result

def species_to_matrix(s):
    mesh = _the_region._mesh
    data = numpy.zeros([len(mesh._xs), len(mesh._ys), len(mesh._zs)]) + numpy.nan
    for node in s._nodes:
        data[node._i, node._j, node._k] = node.molecules
        print node._i, node._j, node._k, node.molecules
    return data

def plot(data):
    from volume_slicer import VolumeSlicer
    """
    # precondition: assumes all data >= 0 with max > 0
    max_value = geom.data[layer].max()
    data = numpy.array(geom.data[layer])
    # center at 0... make outside of the domain an extreme color
    for i, j, k in itertools.product(range(geom.nx), range(geom.ny), range(geom.nz)):
        if not geom.domain[k, j, i]:
            data[k, j, i] = -max_value
    """
    m = VolumeSlicer(data=data)
    m.configure_traits()
    
