"""
Mohammad, McGill University, Canada
July 7, 2014
Rewritten: 11-11-2015
"""

execfile("common_low_neuron.py")
execfile("common_stochastic.py")

save_to_txt(time_steps, Ca_record, "Ca", "low")
# save_to_txt(time_steps, Buf_record, "Buf", "low")
# save_to_txt(time_steps, CaBuf_record, "CaBuf", "low")

# Ending of the simulation with statistics!
finalization_c()
