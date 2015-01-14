===========
admm4block
===========

admm4block provides methods for solving conic programming problems and SDP problems. Typical usage
often looks like this:

    #!/usr/bin/env python

    import admm4block
    
    #generate Copt, Aeq (a list of matrices), beq
    mySDP = admm4block.DNNSDP(Copt, Aeq, beq)
    #decide sigma, tau, tol and nsteps
    mySDP.Solve(sigma, tau, tol, nsteps)

when importing admm4block, the class ConicProgrammingProblem, SDP and DNNSDP are already imported so the user can use::

    admm4block.DNNSDP

instead of::
    
    admm4block.sdp.DNNSDP

conic
=========

conic includes the class ConicProgrammingProblem, with the following methods:

* Solve: Solves the Conic Programming Problem using the ADMM-3-block method explained in [1]. There is no "Step" public function.

* InitialConditions: Creates feasible initial conditions (x0 = Pseudoinverse of Aeq * beq, s0 = the projection on K* of x0, z0 = the projection on Kp* of x0). If any of the initial conditions is preset it is left that way.

* __init__ : Requires a nxm matrix Aeq, an array beq of dimension m, an array Copt of dimension n, and two functions K and Kp that indicate the projections onto the cones K and Kp (the duals are calculated using Moreau's decomposition theorem). In the future will admit inequalities too, for now slack variables must be added by the user.

sdp
======

sdp includes the classes SDP and DNNSDP, with the following methods:

* toConic: Converts the SDP or DNNSDP into a Conic problem and returns an object of the class ConicProgrammingProblem.

* Solve: Solves the SDP or DNNSDP problem and returns [X,s,z,y,res,message] where X is the primal variable, and s,z,y are dual variables.

* __init__ : Requires a list of nxn matrices Aeq, a list beq of the same length, and an nxn matrix Copt (symmetric). The projections onto the cones are automatically generated (thus be careful if you have a condition such as AX <= I)

Build/Install
=================

The installation is pretty straightforward using python setup.py. For best results I recommend using pip with the following commands::

    python setup.py sdist
    cd dist
    pip install admm4block

Another way is simply running::

    python setup.py install

Apparently the first way is better for ensuring easy uninstall in the future with pip uninstall admm4block.

Tests
================

The tests can be run also using setup.py::

    python setup.py test

Or importing them within the Python console.

=================
Known issues (v0.3dev)
=================

* Checking the tolerance in every iteration is not only too costly but also in some way affects the solutions (I have not been able to fix it)
* There is no adaptative control for the parameter sigma.
* The error starts decreasing and then increases (all of it coming from <X,S> = 0 !) after ~ 500 steps. It's possible that this is due to numerical error of the eigh function because it happens when n is larger than 50. The error never goes below 1e-2 for big enough matrices
* There is no multiprocessor implementation of any kind.
* Input and output is bad, for now the user has to create the list of matrices instead of giving normal conditions such as X[i][i] = 1.
* The scripts are supposed to be created but they don't get saved, or I don't know how to run them. Will put them as callable modules (__main__) in the future.
* No SOCP or LP implemented.
