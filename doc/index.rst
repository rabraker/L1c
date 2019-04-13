.. L1c documentation master file, created by
   sphinx-quickstart on Tue Mar 26 19:06:17 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

=====================
Documentation for L1C
=====================

.. toctree::
   :maxdepth: 3
   :hidden:

   building
   api
   usage
   interfaces
   zreferences

.. ------------- Links Used ---------------------------
.. _l1-magic: https://statweb.stanford.edu/\~candes/l1magic

.. _NESTA: http://statweb.stanford.edu/\~candes/nesta
.. _SparseLab: https://sparselab.stanford.edu
.. _SplitBregman: https://www.ece.rice.edu/~tag7/Tom_Goldstein/Split_Bregman.html
.. _Bregman Splitting: https://www.ece.rice.edu/~tag7/Tom_Goldstein/Split_Bregman.html




The goal of this project is to provide a c library that is useful for the application of compressed sensing and related L1-regularized problems.
To date, lots of research has been produced that develops methods to efficiently solve various formulations of the compressed sensing problem.
To their great credit, many of those reasearchers have released code along with their publications.
See, for example, l1-magic_ , NESTA_ , SparseLab_, and SplitBregman_.
So, a good question is, why re-do it in c?

* Much of the published code is written in Matlab. This means if you want to use the code in a Python or R or Julia or a c/c++ project, you must translate the Matlab code to your language. By contrast, c is something like the langua franca of computing. It is relatively straigtforward to interface c code with most languages, including Matlab.

* CS optimizations are computationally intensive. Thus, for those who wish to *apply* CS to real world problems, it is beneficial to write optimized implementations of those optimizations and c fits that requirement (Fortran might be better, but I have forgotten it).

* c is fun (probably the main reason).



Features
--------

Present capabilites and features include

* Solving the basis pursuit denoising problem

  .. math::

     \min_{x} ||x||_{1}  \quad \text{s.t.} \quad ||Ax -b||_{2} < \epsilon,



  using a log-barrier algorithm, which is based on the implementation in l1-magic_,
  .. <<sec:l1qc_mod, with a few modifications>>. 

* Solving the Anistropic Total Variation (TV) denoising using  `Bregman Splitting`_. Given a noisy :math:`n` by :math:`m` image :math:`f`, this solves the problem

.. math::

   \min_{u} ||\nabla_{x}||_{1} + ||\nabla_{y}||_{1} + ||u-f||_{2}. 

* 1D and 2D discrete cosine transforms (using FFTW).

* Python and Matlab bindings for the high-level optimization routines.



Performance
===========
So far, using `l1C` gives me a speed increase of between 2 and 7 times faster compared to the original matlab code, depending on the problem and computer I run it on.

If you compile with FFTW+OpenBlas, it is important that both libraries are compiled with openmp. I don't quite understand what happens, but if this is not the case, I see only single processor being used and performance suffers dramatically. 

If you have a CPU with hyperthreading, it is important to export the environmental variable

`export OMP_NUM_THREADS=N`

where N is the number of *physical* cores. Essentially, if you have HT, this is half the the number of processors you see in a resource monitor, which shows the number of *logical* cores. The code currently can not detect this, and for number crunching applications like this one, HT is detrimental.

Setting `OMP_BIND_PROC=true` seems to cost us about 1 second.


