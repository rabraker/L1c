.. _lin_solve:

Linear Solvers
==============


Many of the optimizations require us to solve the a large system of equations to acquire a descent method. For example, in Newton descent, we solve

.. math::

   H d_x = -g

for :math:`d_x` where :math:`H` is the Hessian and :math:`g` the gradient. For large scale problems, :math:`H` will not fit memory, though we often have a method to compute the product :math:`Hx` efficiently and without storing the entire matrix.

Thus, we resign ourselves to computing an approximate solution. There are many methods to do this. Implemented so far in L1c are conjugate gradient and pre-conditioned conjugate gradient solvers.

Note that basic usage of L1c does not require one to use these functions. However, they seem like they could prove useful to users and so they are exposed.


Linear Solver APIs
------------------

.. doxygengroup:: lin_solve
   :content-only:




