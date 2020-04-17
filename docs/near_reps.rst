Near quadrature correction representation
==========================================

The near quadrature correction is stored in a row-sparse compressed
format as a list between target and patches and stores the following
quantities

    - nnz - number of non-zero target-patch interactions
    - col-ind(nnz) - list of patches which interact with the targets.
      May contain repetitions, see example below
    - row-ptr(ntarg+1) - row_ptr(i) is the starting location in the
      col-ind list where the list of patches in the near field of target
      i start. If row_ptr(i) <= j < row_ptr(i+1), then target(i) and
      patch col_ind(j) are in the near-field of each other
    - nquad - number of non-zero entries in the near field quadrature
      array:

      .. math::

           nquad = \sum_{i=1}^{nnz} m_{\textrm{col-ind(j)}}

    - iquad(nnz) - iquad(i) is the location in the quadrature array
      where the matrix entries corresponding to the interaction of
      target i, and patch col_ind(j) start in wnear array
    - wnear(nquad) - near field quadrature correction array

.. _nearreps-exmp:

Example
-------
Consider the following matrix with 3 targets (rows) and 5 patches
(columns), where $\times$ denotes a combination of patch and target which are
handled through near quadrature correction, and $-$ are handled through
oversampled far quadrature

.. math::
   
   \begin{bmatrix}
   \times & - & - & - & \times \\
   - & \times & - & \times & - \\
   - & - & \times & \times & \times 
   \end{bmatrix}

Then, for this example

- row-ptr = [1,3,5,8]
- col-ind = [1,5,2,4,3,4,5]
