# This file embbeds the c++ file that handles the ADMM updates in python

import cone_SDP_programming as CnP
import Interface as ADMM


MySDP = SnP.DNNSDP(Copt, Aeq, Ain, bin) # get the input of the SDP problem
MyConic = MySDP.toConic() # change the input to the standard conic problem
Solution = ADMM() # update the problem with the ADMM solver
