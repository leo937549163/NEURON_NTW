from neuron import h, rxd, gui

#reaction constants
kf = 1.0 #reaction constant for reaction r0: Ca + Buf = CaBuf
kb = 0.01 #reaction constant for reaction r1: CaBuf = Ca + Buf 

# initial concentrations in milli molar
# initCa = 1.91294084 / 1000
initCa = 2.0 / 1000
# initBuf = 4.0 / 1000
# initCaBuf = 0.0001 / 1000

# Diffusion constants (micro_m2/ms)
# caD = 0.001 # to have a good plot by time
caD = 0.1
# bufD = 0.0
# CaBufD = 0.0