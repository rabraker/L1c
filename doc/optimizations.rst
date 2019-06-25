Optimizations
=============



Nesta
------
This is an implementation of the l1-specialized version of Nesterov's algorithm described in :cite:`becker_nesta_2011`. This method is to be prefered for large scale problems.

.. doxygengroup:: nesta
   :content-only:
   :members:


TV de-noising via Bregman Splitting
-----------------------------------
This algorithm is based on the Goldstien paper :cite:`goldstein_splitbregman_2009`.

.. doxygengroup:: bregman
   :content-only:
   :members:

Quadratically Constrained Log-barrier
--------------------------------------
This is a log-barrier newton method. It is based on the algorithm described in the l1-magic documentation :cite:`l1magic`. This algorithm is not recommended for large scale problems.

.. doxygengroup:: l1qc_lb
   :content-only:
   :members:

