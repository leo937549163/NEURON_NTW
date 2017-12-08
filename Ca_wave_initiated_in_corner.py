"""
Mohammad, McGill University, Canada
July 7, 2014
"""

def get_connectivity(region):
    indices = _get_indices_dict(region)
    xs, ys, zs = region._xs, region._ys, region._zs
    result = []
    for i in xrange(len(xs)):
        result.append(_get_neighbor_indices(xs[i], ys[i], zs[i], indices))
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
_long_ptr = ctypes.POINTER(ctypes.c_long)
_float_ptr = ctypes.POINTER(ctypes.c_float)

#how come you got the values of flat_tems here?
def write_connectivity_flattened_c(flat_items, flat_offsets, stateList, numberOfSpecies, diffusionConstantList,CaWaveFlag):
    """wrapper for the C function"""
    _c_dll.write_connectivity(int(len(flat_offsets) - 1), flat_items.ctypes.data_as(_long_ptr), flat_offsets.ctypes.data_as(_long_ptr), stateList.ctypes.data_as(_float_ptr),numberOfSpecies, diffusionConstantList.ctypes.data_as(_long_ptr))

#QUESTIONS 
def createReaction_c(numberOfReactant,reactantMetricValue,reactantSpeciesID,numberOfProduct,productMetricValue,productSpeciesID,reactionConstant,dx):
    """wrapper for the C function"""
    _c_dll.createReaction(numberOfReactant,reactantMetricValue.ctypes.data_as(_float_ptr),reactantSpeciesID.ctypes.data_as(_long_ptr),numberOfProduct,productMetricValue.ctypes.data_as(_float_ptr),productSpeciesID.ctypes.data_as(_long_ptr),reactionConstant.ctypes.data_as(_float_ptr),dx.ctypes.data_as(_float_ptr))

# how can I send specis info say d[] values to my c code?
def runMainLoop_c(nextTimeStep):
    """wrapper for the C function"""
    _c_dll.runMainLoop(nextTimeStep.ctypes.data_as(_float_ptr))

def initSVs_c(dx):
    """wrapper for the C function"""
    _c_dll.initSVs(dx.ctypes.data_as(_float_ptr))

def finalization_c():
    """wrapper for the C function"""
    _c_dll.finalization()

#
# Please don't call _get_indices_dict and _get_neighbor_indices directly
# They depend on details of the internal implementation and are subject to
# change. They work as of May 2, 2014.
#

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

#
# 		main 		main 		main		main		main		main
#
if __name__ == '__main__':
    from neuron import h, rxd, gui
    from matplotlib import pyplot
    import volume_slicer as vs
    import numpy as np
    
    #CaWaveFlag: define CaWaveFlag where 1 for CaWaveSimulation 0 for any other simulation like Ca buffer model or something else.
    CaWaveFlag=1 # Very important to define it. . . . 0 means it is not Ca wave model
     
    #Cylinder of 5x5 micron is created here 
    dend1 = h.Section()
    dend1.L = 5
    dend1.diam = 5
  
    dend_region = rxd.Region([dend1])
    rxd.set_solve_type(dimension=3)

    Ca_C=rxd.Species(dend_region,name='Ca_C',d=0.075, initial=4)#0.75
    Ca_E=rxd.Species(dend_region,name='Ca_E',d=0.0, initial=200)#200
    IP3=rxd.Species(dend_region,name='IP3',d=2.75, initial=0)#2.75
    inhibit=rxd.Species(dend_region,name='inhibit',d=0.0, initial=1)
    not_inhibit=rxd.Species(dend_region,name='not_inhibit',d=0.0, initial=0)
    not_active=rxd.Species(dend_region,name='not_active',d=0.0, initial=1)
    active=rxd.Species(dend_region,name='active',d=0.0, initial=0)
    
    diffusionConstantList=np.array([Ca_C._d,Ca_E._d,IP3._d,inhibit._d,not_inhibit._d,not_active._d,active._d])    

    h.finitialize()
    
    connectivity = get_connectivity(dend_region)
    
    # reaction rates
    k_rc1 = 1.0
    k_rc2 = 1.0
    k_rc3 = 1.0
    k_rc4 = 1.0
    k_rc5 = 10.0#10
    k_rc6 = 0.4 #.5
    k_rc7 = 0.70#1.5
    

    IP3_initial=4150.0	
    
    '''
    # IP3 is initiated in the centre of 5x5 micron^3 cylindre
    IP3.nodes[3749].value=IP3_initial
    IP3.nodes[3750].value=IP3_initial
    IP3.nodes[3751].value=IP3_initial
    IP3.nodes[3752].value=IP3_initial
    IP3.nodes[3753].value=IP3_initial
    '''
    
    # IP3 is initiated in the corner of 5x5 micron^3 cylindre
    IP3.nodes[0].value=IP3_initial
    IP3.nodes[1].value=IP3_initial
    IP3.nodes[2].value=IP3_initial
    IP3.nodes[3].value=IP3_initial
    IP3.nodes[4].value=IP3_initial   
	
    
    numberOfSpecies=7 #number of species   
    
    #data, dataIndeces getting from here . . .
    flat_items, flat_offsets = flatten(connectivity)

    # sending geometry and specis informations to stochastic simulator
    write_connectivity_flattened_c(flat_items, flat_offsets, rxd.node._states, numberOfSpecies, diffusionConstantList,CaWaveFlag)
    
    dx_v = dend_region._dx  
    dx=np.array([dx_v])
    
    
    #dynamics for reaction 1: 
    numberOfReactant=1
    reactantMetricValue=np.array([1.0])
    reactantSpeciesID=np.array([3])
    numberOfProduct=1
    productMetricValue=np.array([1.0])
    productSpeciesID=np.array([4])
    reactionConstant=np.array([k_rc1])
    #sending reaction 1
    createReaction_c(numberOfReactant,reactantMetricValue,reactantSpeciesID,numberOfProduct,productMetricValue,productSpeciesID,reactionConstant,dx)

    #dynamics for reaction 2: 
    numberOfReactant=1
    reactantMetricValue=np.array([1.0])
    reactantSpeciesID=np.array([4])
    numberOfProduct=1
    productMetricValue=np.array([1.0])
    productSpeciesID=np.array([3])
    reactionConstant=np.array([k_rc2])
    #sending reaction 2
    createReaction_c(numberOfReactant,reactantMetricValue,reactantSpeciesID,numberOfProduct,productMetricValue,productSpeciesID,reactionConstant,dx)
    
    #dynamics for reaction 3: 
    numberOfReactant=1
    reactantMetricValue=np.array([1.0])
    reactantSpeciesID=np.array([5])
    numberOfProduct=1
    productMetricValue=np.array([1.0])
    productSpeciesID=np.array([6])
    reactionConstant=np.array([k_rc3])
    #sending reaction 3
    createReaction_c(numberOfReactant,reactantMetricValue,reactantSpeciesID,numberOfProduct,productMetricValue,productSpeciesID,reactionConstant,dx)

    #dynamics for reaction 4: 
    numberOfReactant=1
    reactantMetricValue=np.array([1.0])
    reactantSpeciesID=np.array([6])
    numberOfProduct=1
    productMetricValue=np.array([1.0])
    productSpeciesID=np.array([5])
    reactionConstant=np.array([k_rc4])
    #sending reaction 4
    createReaction_c(numberOfReactant,reactantMetricValue,reactantSpeciesID,numberOfProduct,productMetricValue,productSpeciesID,reactionConstant,dx)
  
    #dynamics for reaction 5: Releasing
    numberOfReactant=1
    reactantMetricValue=np.array([2.0])
    reactantSpeciesID=np.array([1])
    numberOfProduct=1
    productMetricValue=np.array([2.0])
    productSpeciesID=np.array([0])
    reactionConstant=np.array([k_rc5])
    #sending reaction 5
    createReaction_c(numberOfReactant,reactantMetricValue,reactantSpeciesID,numberOfProduct,productMetricValue,productSpeciesID,reactionConstant,dx) 
    
    #dynamics for reaction 6: Pumping
    numberOfReactant=1
    reactantMetricValue=np.array([1.0])
    reactantSpeciesID=np.array([0])
    numberOfProduct=1
    productMetricValue=np.array([0.0])#Asol kaj ekhane hoche . . . .
    productSpeciesID=np.array([1])   
    reactionConstant=np.array([k_rc6])
    #sending reaction 6
    createReaction_c(numberOfReactant,reactantMetricValue,reactantSpeciesID,numberOfProduct,productMetricValue,productSpeciesID,reactionConstant,dx)
    
    #dynamics for reaction 7: Leaking
    numberOfReactant=1
    reactantMetricValue=np.array([1.0])
    reactantSpeciesID=np.array([1])
    numberOfProduct=1
    productMetricValue=np.array([1.0])
    productSpeciesID=np.array([0])
    reactionConstant=np.array([k_rc7])
    #sending reaction 7 to NTW
    createReaction_c(numberOfReactant,reactantMetricValue,reactantSpeciesID,numberOfProduct,productMetricValue,productSpeciesID,reactionConstant,dx)
       
     
    #dx=np.array([dx_v])# this value is dx = dy = dz = 0.25
    initSVs_c(dx)  #that is c code . .  
    
    time_steps=[]
    Ca_C_record=[]
    Ca_E_record=[]
    
    simulation_time=0.0
    flag=1
    i=0
    
    filter_n=5
    sample_counter=0
    timeArray=np.array([0.1])
    runMainLoop_c(timeArray)
    
    #for 3d figure ! 
    print 'Simulaiton Advanced :', simulation_time, 'ms i.e. initial states in which IP3 is initiated in the corner of the cylindre'   
    mesh = dend_region._mesh
    data = np.zeros([len(mesh._xs), len(mesh._ys), len(mesh._zs)]) + np.nan
    for node in Ca_C.nodes:
        data[node._i, node._j, node._k] = node.value    
    m=vs.VolumeSlicer(data=data)
    m.configure_traits()    
       
    step=1.0;
    #print Ca_C.nodes[0].value
    while simulation_time < 7.0:
	    simulation_time=simulation_time+step
	    timeArray=np.array([simulation_time])
	    runMainLoop_c(timeArray)
	    time_steps.append(simulation_time)
	    print 'Simulaiton Advanced :', simulation_time, 'ms'
	    #for 3d figure ! 
	    mesh = dend_region._mesh
	    data = np.zeros([len(mesh._xs), len(mesh._ys), len(mesh._zs)]) + np.nan
	    for node in Ca_C.nodes:
		  data[node._i, node._j, node._k] = node.value
	    m=vs.VolumeSlicer(data=data)
	    m.configure_traits()
	     
    finalization_c()# Ending of the simulation with statistics !
