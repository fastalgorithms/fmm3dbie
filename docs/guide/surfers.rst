
.. role:: matlab(code)
   :language: matlab   

Surface Discretization with Surfers
====================================

In fmm3dBIE, a smooth, regular surface is discretized by dividing
it into pieces, called "patches", which are then represented
by polynomial interpolants using appropriate discretization nodes. This
information is stored in a :matlab:`surfer` object.

The :matlab:`surfer` class documentation gives a survey of the available methods:

.. include:: ../../matlab/@surfer/surfer.m
   :literal:
   :code: matlab
   :start-after: classdef surfer
   :end-before: % author

.. note::

   To obtain the documentation of a class method which has
   overloaded the name of a MATLAB built-in, use the syntax
   :matlab:`help class_name/method_name`. For example:

   .. code:: matlab

      help surfer/area
