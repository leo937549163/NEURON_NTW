cdef extern from "gperftools/profiler.h":
	int ProfilerStart(char* fname)
	void ProfilerStop()
	void ProfilerFlush()

def ProfStart(fname):
	return ProfilerStart(fname)

def ProfStop():
	ProfilerStop()

def ProfFlush():
	ProfilerFlush()