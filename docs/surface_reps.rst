Surface representation
=======================

Recall that the surface $\Gamma$ be defined via a collection of patches
$\Gamma_{j} = 1,2,\ldots N_{\textrm{patches}}$, where each patch 
$\Gamma_{j}$ is parametrized by a non-degenerate chart 
$\boldsymbol{x}^{j}: B \to \Gamma_{j}$, where $B$ is standard base
element. 
Each patch $\Gamma_{j}$ is sampled at a collection of discretization nodes
$(u_{i},v_{i}) \in B$, $i=1,2,\ldots m_{j}$. Let $N$ denote the total
number of discretization points, then we store the following quantities
to represent the surfaces

    - N: number of discretization points
    - $N_{\textrm{patches}}$: number of patches
    - norders($N_{\textrm{patches}}$): order of discretization of each 
      patch
    - srcvals (12,N): surface samples of $\boldsymbol{x}^{j},
      \partial_{u} \boldsymbol{x}^{j}, \partial_{v} \boldsymbol{x}^{j},$
      and $\boldsymbol{n}^{j}$, where
    
    .. math::
   
        \boldsymbol{n}^{j} = \frac{\partial_{u} \boldsymbol{x}^{j} \times 
        \partial_{v} \boldsymbol{x}^{j}}{|\partial_{u} \boldsymbol{x}^{j}
        \times \partial_{v} \boldsymbol{x}^{j}|}

    - srccoefs (9,N): Orthogonal polynomial expansions of 
      $\boldsymbol{x}^{j}, \partial_{u} \boldsymbol{x}^{j}$, 
      and $\partial_{v} \boldsymbol{x}_{j}$
    - iptype ($N_{\textrm{patches}}$): patch type
    - ixyzs ($N_{\textrm{patches}} + 1$): location in srcvals, and
      srccoefs where information for patches begin. Also implicitly
      stores $m_{j}$ = ixyzs(j+1)-ixyzs(j)

Supported base elements and discretization nodes
-------------------------------------------------

- iptype = 1: 

  .. math::

    T_{0} = \{ (u,v): u>0,v>0, u+v<1 \},

  discretized using Vioreanu Rokhlin nodes (up to order 20), 
  and the basis functions are Koornwinder polynomial expansions. 
  For $norder=p$, there are $(p+1) \cdot (p+2)/2$ discretization nodes.


- iptype = 11: 

  .. math::

    Q_{0} = \{ (u,v) \in (-1,1)^2  \},

  discretized with tensor product Gauss-Legendre nodes, and 
  the basis functions are tensor product Legendre polynomials. 
  For $norder=p$, there are $(p+1)^2$ discretization nodes.


- iptype = 12: 

  .. math::

    Q_{0} = \{ (u,v) \in (-1,1)^2  \},

  discretized with tensor product Chebyshev nodes, and 
  the basis functions are tensor product Chebyshev polynomials. 
  For $norder=p$, there are $(p+1)^2$ discretization nodes.


Supported input formats 
--------------------------

Here are the list of input formats that are currently supported

- :ref:`.go3`

.. _.go3:

.go3
*****
The .go3 file format is a storage format where each patch is iptype=1,
and discretized using the same order Vioreanu-Rokhlin nodes.:: 

    norder-1
    N_{patches}
    srcvals(1,1)
    srcvals(1,2)
        .
        .
        .
    srcvals(1,N)    
    srcvals(2,1)
        .
        .
        .
    srcvals(12,N)

See ``geometries/sphere_192_o03.go3`` for an example.



