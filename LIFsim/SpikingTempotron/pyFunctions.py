import numpy as np
import pickle as __pkl
#from _ST import __SV
from . import _ST


def __array_w(self):
    return np.array(self.w())

def __SpikingTempotron_getVoltageTrace_1(self,il,mu):
	il.setCurrentPattern(mu)
	self.setPattern(il.getP2CurrentPattern())
	self.restart()
	trace = self._getVoltageTrace();
	s=len(trace)
	t=np.zeros(s)
	V=np.zeros(s)
	for ind,p in enumerate(trace):
		t[ind]=p.time
		V[ind]=p.V
	return t,V

def __SpikingTempotron_getVoltageTrace(self,pat):
	self.setPattern(pat)
	self.restart()
	trace = self._getVoltageTrace();
	s=len(trace)
	t=np.zeros(s)
	V=np.zeros(s)
	for ind,p in enumerate(trace):
		t[ind]=p.time
		V[ind]=p.V
	return t,V


def times2SpikeTrain(pat_times,PSP_Kernel=None,time_block_size = None):
    st = _ST.SpikeTrain()
    for i,spikes in enumerate(pat_times):
        for t_i in spikes:
            st._cpp_add_spike(i,t_i)
	#st.extend([_ST.Spike2(i,spike_time) for spike_time in spikes])
    st.sort()
    if time_block_size != None:
        st.timeBlockSize = time_block_size
        st.applyTimeBlock()
    if PSP_Kernel != None :
        st.calcExponents(PSP_Kernel)
    return st

#@property
#def __FE_ES_Proj(self):
#    return np.array(self._Proj)
#    #return np.reshape(self._Proj, (self.nSynapse+1,self.nSynapses+1) , order='F')

#@__FE_ES_Proj.setter
#def __FE_ES_Proj(self,P):
#    '''self._Proj = _ST.__vecvecdouble()
#    l = []
#    for row in P:
#        v = _ST.__vecdub()
#        v.extend(row)
#        l.append(v)
#    self._Proj.extend(l)'''
#    self.Proj = _ST.__vecdub()
#    #print P.flatten(order='F').shape
#    self._Proj.extend(P.flatten('F'))
#    #self._Proj.extend(P.transpose()[np.tril_indices_from(P)])
#    #print len(self._Proj)

@property
def __tempotron__w(self):
	return np.array([w for w in self._w])
@__tempotron__w.setter
def __tempotron__w(self,w):
	self._w = w

@property
def __HH__w(self):
	'''returns a copy of the weights'''
	return np.array(self._w)
@__HH__w.setter
def __HH__w(self,w):
	if len(w) != len(self._w):
		raise RuntimeError("Trying to set weights with wrong number of synapses")
	self._w = _ST.__vecdub()
	self._w.extend(w)
