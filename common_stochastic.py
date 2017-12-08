from collections import defaultdict
import time
import numpy as np
from gperftools_wrapped import ProfStart, ProfStop

def get_connectivity(region):
    indices = _get_indices_dict(region)
    xs, ys, zs = region._xs, region._ys, region._zs
    result = []
    for i in xrange(len(xs)):
        result.append(_get_neighbor_indices(xs[i], ys[i], zs[i], indices))
    return result

def get_connectivity_1d(n):
    result = []
    for i in xrange(n):
        tmp = []
        if i >= 1:
            tmp.append(i-1)
        if i <= n-2:
            tmp.append(i+1)
        result.append(tmp)
    return result

def write_connectivity_structured(f, connectivity):
    """write connectivity to f from the structured data (list of lists)"""
    for row in connectivity:
        f.write(' '.join(str(index) for index in row) + '\n')

def flatten(list_of_lists):
    """turns a list of lists into 2 numpy arrays of integers

    The first return value is the data, the second is the first index of each
    list. The second array is one longer than the first to indicate the le
    """
    import numpy
    import itertools
    data = numpy.array(list(itertools.chain(*list_of_lists)), dtype=int)
    start_indices = [0]
    for item in list_of_lists:
        start_indices.append(len(item) + start_indices[-1])
    start_indices = numpy.array(start_indices, dtype=int)
    return data, start_indices

def length(flat_offsets, list_id):
    """return the length of an item from a flattened list of lists"""
    return flat_offsets[list_id + 1] - flat_offsets[list_id]

def item(flat_items, flat_offsets, list_id, item_id):
    """return a selected item from a flattened list of lists"""
    return flat_items[flat_offsets[list_id] + item_id]

def write_connectivity_flattened(f, flat_items, flat_offsets):
    """write connectivity to f from flattened data"""
    for list_id in xrange(len(flat_offsets) - 1):
        f.write(' '.join(str(item(flat_items, flat_offsets, list_id, item_id))
                for item_id in xrange(length(flat_offsets, list_id))) + '\n')

#
# enable connectivity with C
#

import ctypes
"""
If you're accessing a function inside NEURON, you wouldn't create a separate so.
Just get access to NEURON via:
    import neuron
    nrn_dll = neuron.nrn_dll()
"""

_c_dll = ctypes.CDLL('./connectivity.so')
reaction_function = ctypes.CDLL("./reaction_function.so")
_long_ptr = ctypes.POINTER(ctypes.c_long)
_float_ptr = ctypes.POINTER(ctypes.c_float)

#how come you got the values of flat_tems here?
def write_connectivity_flattened_c(flat_items, flat_offsets, stateList, numberOfSpecies, diffusionConstantList, CaWaveFlag):
    """wrapper for the C function"""
    #print '_c_dll.write_connectivity(%r, array(%r), array(%r), array(%r), %r, array(%r), %r)' % (int(len(flat_offsets) - 1), list(flat_items), list(flat_offsets), list(stateList), numberOfSpecies, list(diffusionConstantList), CaWaveFlag)
    _c_dll.write_connectivity(int(len(flat_offsets) - 1), flat_items.ctypes.data_as(_long_ptr), flat_offsets.ctypes.data_as(_long_ptr), stateList.ctypes.data_as(_float_ptr), numberOfSpecies, diffusionConstantList.ctypes.data_as(_long_ptr), CaWaveFlag)

#QUESTIONS
def createReaction_c(numberOfReactant,reactantMetricValue,reactantSpeciesID,numberOfProduct,productMetricValue,productSpeciesID,reactionConstant,dx,reactionFunction):
    """wrapper for the C function"""
    #print '_c_dll.createReaction(%r,%r,%r,%r,%r,%r,%r)' % (numberOfReactant,reactantMetricValue,reactantSpeciesID,numberOfProduct,productMetricValue,productSpeciesID,reactionConstant)
    _c_dll.createReaction(numberOfReactant,reactantMetricValue.ctypes.data_as(_float_ptr),reactantSpeciesID.ctypes.data_as(_long_ptr),numberOfProduct,productMetricValue.ctypes.data_as(_float_ptr),productSpeciesID.ctypes.data_as(_long_ptr),reactionConstant.ctypes.data_as(_float_ptr),dx.ctypes.data_as(_float_ptr), reactionFunction)

# how can I send specis info say d[] values to my c code?
def runMainLoop_c(nextTimeStep):
    """wrapper for the C function"""
    #print '_c_dll.runMainLoop(%r)' %  nextTimeStep
    _c_dll.runMainLoop(nextTimeStep.ctypes.data_as(_float_ptr))

def initSVs_c(dx):
    """wrapper for the C function"""
    #print '_c_dell.initSVs(%r)' % dx
    _c_dll.initSVs(dx.ctypes.data_as(_float_ptr))

def finalization_c():
    """wrapper for the C function"""
    #print '_c_dll.finalization()'
    _c_dll.finalization()

def setRandomSeed_c(SEED):
    """wrapper for the C function"""
    #print '_c_dll.setRandomSeed(%r)' % SEED
    _c_dll.setRandomSeed(SEED)

#
# Please don't call _get_indices_dict and _get_neighbor_indices directly
# They depend on details of the internal implementation and are subject to change. They work as of May 2, 2014.
#

def _get_indices_dict(region):
    indices = defaultdict(int)
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



sim_type = "1d" # "3d"

if sim_type == "1d":
    secs = [h.Section()]
    for sec in secs:
        sec.nseg = 401
        sec.L = 100.25
        sec.diam = 1
    dend_region = rxd.Region(secs)
    rxd.set_solve_type(dimension=1)
elif sim_type == "3d":
    secs = [h.Section()]
    for sec in secs:
        sec.nseg = 401
        sec.L = 100
        sec.diam = 1
    dend_region = rxd.Region(secs)
    rxd.set_solve_type(dimension=3)

#CaWaveFlag: define CaWaveFlag where 1 for CaWaveSimulation 0 for any other simulation like Ca buffer model or something else.
CaWaveFlag=0 # Very important to define it. . . . 0 means it is not Ca wave model

#Defining the Number of species in the simulation
numberOfSpecies=1 ###########################################

Ca = rxd.Species(dend_region, name='Ca', d=caD) # milli molar
# Buf = rxd.Species(dend_region, name='Buf', d=bufD)
# CaBuf = rxd.Species(dend_region, name='CaBuf', d=CaBufD)

from collections import defaultdict
x_sorted = sorted(set([node.x3d for node in Ca.nodes]))
x_indice = [[] for _ in xrange(len(x_sorted))]
for i, node in enumerate(Ca.nodes):
    x_indice[x_sorted.index(node.x3d)].append(i)

x_selected = np.median(list(set([node.x3d for node in Ca.nodes])))
# print [i if Ca.nodes[i].x3d == x_median else 0 for i in xrange(len(Ca.nodes))]

Ca.initial = lambda node: (initCa if node.x3d <= x_selected else 0.1/1000)
# Buf.initial = lambda node: (initBuf if node.x3d == x_selected else 0)
# CaBuf.initial = lambda node: (initCaBuf if node.x3d == x_selected else 0)

Ca_record = [[initCa if x <= x_selected else 0.1/1000 for x in x_sorted]]
# Buf_record = [[initBuf if x == x_selected else 0 for x in x_sorted]]
# CaBuf_record = [[initCaBuf if x == x_selected else 0 for x in x_sorted]]

#creating a diffusion constant List to send to NTW...
# diffusionConstantList=np.array([Ca._d, Buf._d, CaBuf._d])
diffusionConstantList=np.array([Ca._d])
print diffusionConstantList

#must initialize before den_region._dx call
h.finitialize()



# Geometry's neighbor information that is list of neighbor lists are created . . . .
# connectivity = get_connectivity(dend_region)
if sim_type == "1d":
    connectivity = get_connectivity_1d(len(Ca.nodes))
elif sim_type == "3d":
    connectivity = get_connectivity(dend_region)

# print [node.x3d for node in Ca.nodes]
# print "connectivity", connectivity

#data, dataIndeces getting from here . . .
flat_items, flat_offsets = flatten(connectivity)

# sending geometry and specis informations to stochastic simulator
write_connectivity_flattened_c(flat_items, flat_offsets, rxd.node._states, numberOfSpecies, diffusionConstantList, CaWaveFlag)
dx_v = dend_region._dx
dx=np.array([dx_v])

wholeVolume=dx_v*dx_v*dx_v*len(Ca.nodes)

setRandomSeed_c(int(round(time.time()*1000)))
# setRandomSeed_c(123)

# from mpl_toolkits.mplot3d import Axes3D
# import matplotlib
# matplotlib.use('agg')
# import matplotlib.pyplot as plt

# x = [node.x3d for node in Ca.nodes]
# y = [node.y3d for node in Ca.nodes]
# z = [node.z3d for node in Ca.nodes]

# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')

# ax.scatter(y, x, z)

# ax.set_xlabel('Y Label')
# ax.set_ylabel('X Label')
# ax.set_zlabel('Z Label')

# plt.show()
# plt.savefig("1.png")
#
# Creating and sending reactions to NTW
#

#Reaction-0: Ca + Buf ----> CaBuf
numberOfReactant=1
reactantMetricValue=np.array([0.0])
reactantSpeciesID=np.array([0])
numberOfProduct=1
productMetricValue=np.array([1.0])
productSpeciesID=np.array([0])
reactionConstant=np.array([kf])# kf = 1.0
#sending reaction 1 to
createReaction_c(numberOfReactant,reactantMetricValue,reactantSpeciesID,numberOfProduct,productMetricValue,productSpeciesID,reactionConstant,dx,reaction_function.reaction_2)


#Reaction-1: CaBuf -----> Ca + Buf
# numberOfReactant=1
# reactantMetricValue=np.array([1.0])
# reactantSpeciesID=np.array([2])
# numberOfProduct=2
# productMetricValue=np.array([1.0,1.0])
# productSpeciesID=np.array([0,1])
# reactionConstant=np.array([kb])# kb = 1.0
# #sending reaction 2 to NTW
# createReaction_c(numberOfReactant,reactantMetricValue,reactantSpeciesID,numberOfProduct,productMetricValue,productSpeciesID,reactionConstant,dx,reaction_function.reaction_1)

#Subvolume's dx, dy, and dz are sending to NTW from here. here dx=dy=dz
#dx_v = dend_region._dx
#dx=np.array([dx_v])
initSVs_c(dx)#sending to NTW here. . .

CONVERSION_FACTOR = 602214.129
CONVERSION_FACTOR_VOLUME = CONVERSION_FACTOR / 64.0

#Some data structures for statistics
time_steps=[]

simulation_time=0.0
step=30
stop_time=30000.0
time_steps.append(simulation_time)

def get_value_list(nodes):
    # c = np.array([node.value/CONVERSION_FACTOR_VOLUME for node in nodes])
    # return [np.average(c[ind]) for ind in x_indice]
    c = np.array([node.value for node in nodes])
    return [np.average(c[ind]) for ind in x_indice]

#print Ca.nodes[0].value
while simulation_time < stop_time :
    simulation_time=simulation_time+step
    timeArray=np.array([simulation_time])


    # ProfStart("output.prof")
    runMainLoop_c(timeArray)
    # ProfStop()

    time_steps.append(simulation_time)
    Ca_record.append(get_value_list(Ca.nodes))

    if simulation_time % 10 == 0:
        print simulation_time
    # Buf_record.append(get_value_list(Buf.nodes))
    # CaBuf_record.append(get_value_list(CaBuf.nodes))

def save_to_txt(time, var, var_name, concentration="low"):
    with open("./Data/Stochastic/States_"+concentration+"_"+var_name+".txt",'w') as f:
        # First Col: Time
        # All the other Col: Data
        for i in xrange(len(time)):
            f.write(str(time[i])+'\t'+'\t'.join(str(x) for x in list(var[i]))+"\n")
