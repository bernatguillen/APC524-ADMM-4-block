# This file embbeds the c++ file that handles the ADMM updates in python

import input as INP
import SDP as SDP
import Interface as ADMM


MySDP = DNNSDP(Copt, Aeq, Ain, bin) # get the input of the SDP problem
MyConic = MySDP.toConic() # change the input to the standard conic problem
Solution = ADMM() # update the problem with the ADMM solver
