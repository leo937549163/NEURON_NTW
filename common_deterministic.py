import numpy as np

secs = [h.Section()]
# for sec in secs:
#     sec.nseg = 41
#     sec.L = 10.25
#     sec.diam = 1
for sec in secs:
    sec.nseg = 401
    sec.L = 100.25
    sec.diam = 1

dend_region = rxd.Region(secs)
rxd.set_solve_type(dimension=1)

Ca = rxd.Species(dend_region, name='Ca', d=caD) # milli molar
# Buf = rxd.Species(dend_region, name='Buf', d=bufD)
# CaBuf = rxd.Species(dend_region, name='CaBuf', d=CaBufD)

# print sorted([node.x3d for node in Ca.nodes])
x_selected = np.median(list(set([node.x3d for node in Ca.nodes]))) # np.median

Ca.initial = lambda node: (initCa if node.x3d <= x_selected else 0.1/1000)
# Buf.initial = lambda node: (initBuf if node.x3d==x_selected else 0)
# CaBuf.initial = lambda node:(initCaBuf if node.x3d==x_selected else 0)

# buffering = rxd.Reaction(Ca + Buf <> CaBuf, kf, kb)

# All in milli molar
# High   2.0/1000
# Low    0.1/1000
# Middle 0.6/1000
# h_ca, l_ca, m_ca = 1.91294084/1000, 0.106274491/1000, 0.637646946/1000
h_ca, l_ca, m_ca = 2.0/1000, 0.1/1000, 0.6/1000
r_123_123 = rxd.Rate(Ca, (l_ca-Ca)*(h_ca-Ca)*(m_ca-Ca)*1000)

# print [node.x3d for node in Ca.nodes]

def get_record_list(nodes):
    record = []
    for node in nodes:
        record.append(h.Vector())
        record[-1].record(node._ref_concentration)
    return record

h.init()

record_Ca = get_record_list(Ca.nodes)
# record_Buf = get_record_list(Buf.nodes)
# record_CaBuf = get_record_list(CaBuf.nodes)

times = h.Vector()
times.record(h._ref_t)

h.init()

# simulate 300 ms
h.dt = 1.0
# h.CVode().fixed_step()
# h.CVode().active(True)
# h.CVode().atol(1e-8) # absolute tolerance

# h.continuerun(300)
h.continuerun(30000)

def save_to_txt(var, var_name, concentration="low"):
    with open('./Data/Deterministic/States_'+concentration+'_'+var_name+'.txt', 'w') as f:
        # First Col: Time
        # All the other Col: Data
        for data in zip(times, *var):
            f.write("\t".join(str(x) for x in list(data))+"\n")
