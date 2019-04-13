.. _transforms:

Linear Transforms
==================


Compressed sensing revolves around finding sparse solutions to the system

.. math::

   y = Ax

where :math:`y\in\mathbb{R}^n` are the observations and :math:`x\in\mathbb{R}^m` is the true signal.
In general, :math:`n\lt\lt m`. The transform :math:`A` can take various forms (see, e.g., :cite:`candes_stable_2006`)


* The elements of :math:`A` are iid guassian.
  This type of transform can be implemented in L1c using the general matrix transforms. See :cpp:func:`l1c_setup_matrix_transforms`.

* Fourier Ensemble. L1c does not yet implement this.

* General orthogonal measurement ensemble (a generalization of the Fourier ensemble).
  Here, we take measurements in the time or spatial domain according to a measurement matrix :math:`\Phi`

  .. math::

     y = \Phi x


  and assume they will be (approximately) sparse in a different basis, e.g., a wavelet or discrete cosine transform (DCT) basis. Let :math:`M` map a vector in the sparsity basis back to the time (resp., spatial) domain. Then we want to solve


  .. math::
     
     y & = \Phi x\\
       & = \Phi M z

  Thus far, L1c implements a special version of this where :math:`\Phi=E` is the identity with rows removed and :math:`M` is the IDCT or 2D IDCT. There are many other possibilities here.


If you wish to add additional transforms, see `dct1.c` or `matrix_transforms.c` and the struct :cpp:type:`l1c_AxFuns`.


Transform API
-------------

.. doxygengroup:: transforms
   :content-only:


