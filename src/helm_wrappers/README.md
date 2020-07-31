This folder contains several routines for evaluating various layer
potentials for solving Helmholtz boundary value problems.

Recall that the general nomenclature for routine types is the following

<routine type prefix>_<pde name>_<input representation>_<data type on
output>_<routine type suffix>

For Helmholtz boundary value problems <pde name>=helm

Input representation
======================

Input representation = comb
-----------------------------
Combined field representation

u = \alpha S_{k} + \beta D_{k}

Input representation = rpcomb
-------------------------------
Right preconditioned combined field representation

u = S_{k} + i \alpha D_{k} [S_{ik}]


Output datatype 
========================

output datatype = dir
--------------------------

Returns potential

output datatype = neu
--------------------------

Returns normal derivative of data (the target points must be restricted
to the boundary)

output datatype = grad
----------------------------

Returns gradient of the potential


Currently only "comb_dir" combination, and "rpcomb_dir/neu" combinations
are supported

