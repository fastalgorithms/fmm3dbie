Locally corrected quadratures
------------------------------

Let the surface $\Gamma$ be defined via a collection of patches
$\Gamma_{j} = 1,2,\ldots N_{\textrm{patches}}$, where each patch 
$\Gamma_{j}$ is parametrized by a non-degenerate chart 
$\boldsymbol{x}^{j}: B \to \Gamma_{j}$, where $B$ is standard base
element. Given a density $\sigma$ defined on $\Gamma$, consider the 
evaluation of the layer potential $\mathcal{S}[\sigma](x)$ at a
collection of targets $x=\boldsymbol{t}_{j}$, $j=1,2\ldots n_{t}$. These
targets could be anywhere in the volume including on the surface
$\Gamma$. 

.. math::

   \mathcal{S}[\sigma](x_{i}) = 
   \int_{\Gamma}  G(x_{i},y) \sigma(y) da(y)
   = \sum_{j=1}^{N_{\textrm{patches}}} \int_{\Gamma_{j}} G_{k}(x_{i},y)
   \sigma(y) da(y) 

If a patch $\Gamma_{j}$ is close to (as compred to the size of the patch) 
a target $x_{i}$, then the integrand
is nearly singular and the integral becomes difficult to evaluate
accurately as compared to when $x_{i}$ is far from the patch
$\Gamma_{j}$. Locally corrected quadrature methods precompute the
quadrature for all near interactions between patches and targets, and
use appropriately oversampled quadratures for the rest of the
interactions.



Near-far split
==================

Let $c_{j}$ denote the centroid of a patch given by

.. math::
   
   c_{j} = \int_{\Gamma_{j}} y da(y) \, ,

and let $R_{j}$ denote the smallest radius $R$ such that a sphere of
radius R centered at $c_{j}$ completely contains $\Gamma_{j}$, i.e.

.. math::

   R_{j} = \min_{R} \{ R : \Gamma_{j} \subset B_{R}(c_{j}) \} \, .

Then given $\eta>0$, the $\eta$-scaled near field of the patch
$\Gamma_{j}$ is given by

.. math::

   N_{\eta}(\Gamma_{j}) = \{ x : d(c_{j},x) \leq \eta R_{j} \} \, .

Given $N_{\eta}(\Gamma_{j})$, let $T_{\eta}(x_{i})$ denote the dual
list -- the collection of patches $\Gamma_{j}$ for which $x_{i}$ is in
its $\eta$-scaled near field,

.. math::
   
   T_{\eta}(x_{i}) = \{ \Gamma_{j} : x_{i} \in N_{\eta}(\Gamma_{j} \} =
   \{ \Gamma_{j} : d(x_{i},c_{j}) \leq \eta R_{j} \}

The integral for $\mathcal{S}[\sigma](x)$ can be split into two parts

.. math::

   \mathcal{S}[\sigma](x) &= \sum_{\Gamma_{j} \in T_{\eta}(x)}
   \int_{\Gamma_{j}} G(x,y)\sigma(y) da(y) + 
   \sum_{\Gamma_{\ell} \not \in T_{\eta}(x)} \int_{\Gamma_{j}}
   G(x,y)\sigma(y) da(y) \\
   &= \mathcal{S}_{\textrm{near}}[\sigma](x) +
   \mathcal{S}_{\textrm{far}}[\sigma](x)

