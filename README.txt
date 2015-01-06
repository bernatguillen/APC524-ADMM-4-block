===========
admm4block
===========

admm4block provides methods for solving conic programming problems and doubly non-negative SDP problems. Typical usage
often looks like this::

    #!/usr/bin/env python

    from admm4block import conic
    from admm4block import dnnsdp
    
    #generate Copt, Aeq, beq, Ain, bin
    mySDP = dnnsdp.DNNSDP(Copt, Aeq, beq, Ain, bin)
    mySDP.Solve(sigma, tau, tol, nsteps)

(Note the double-colon and 4-space indent formatting above.)

Paragraphs are separated by blank lines. *Italics*, **bold**,
and ``monospace`` look like this.


conic
=========

conic includes the class ConicProgrammingProblem, with the following methods:

* Solve

* InitialConditions


dnnsdp
======

dnnsdp includes the class DNNSDP, with the following methods:

* toConic

* Solve
