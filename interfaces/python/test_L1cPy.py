"""
Unittests for L1cPy. To run this, execute
>>> python3 test_L1cPy.py -v
"""

import unittest
import json
import numpy as np
from numpy.testing import assert_almost_equal
from numpy.linalg import norm

import _l1cPy_module as l1cPy


class TestBregman(unittest.TestCase):
    """These are tests for the .
    """

    def setUp(self):
        self.fpath = 'example_img_data.json'

        with open(self.fpath) as json_data:
            d = json.load(json_data)

        img_vec = np.array(d['x_orig'], ndmin=1)
        n = int(np.sqrt(len(img_vec)))
        m = n

        np.random.seed(0)
        img_vec_noisy = img_vec + np.random.rand(n*m)*0.25
        self.Img_noisy = np.reshape(img_vec_noisy, (n, m))
        self.Img_orig = np.reshape(img_vec, (n, m))
        self.n = n
        self.m = m


    def test_breg_dims(self):
        """Check that we get the right dimensions back out.
        Espcially, with a non-square input.
        """

        Img_sub = self.Img_noisy[1:, :]
        n, m = Img_sub.shape
        Img_clean = l1cPy.breg_anistropic_TV(Img_sub,
                                             max_iter=100, max_jac_iter=1,
                                             tol=0.001, mu=5)
        self.assertEqual(n, Img_clean.shape[0])
        self.assertEqual(m, Img_clean.shape[1])


    def test_breg_errors(self):
        """Check that we get errors when we should."""

        wrong_data=np.random.rand(3,3)
        with self.assertRaises(IndexError):
            l1cPy.breg_anistropic_TV(5)
            # , max_iter=100, max_jac_iter=1,
                                                 # tol=0.001, mu=5)
        with self.assertRaises(TypeError):
            l1cPy.breg_anistropic_TV(self.Img_noisy, max_iter=wrong_data)

        with self.assertRaises(TypeError):
            l1cPy.breg_anistropic_TV(self.Img_noisy, wrong_data)



    def test_bregman_anisTV(self):
        """
        Check the bregman TV denoising produces a reasonable output.
        Not obvious to me what else we can check.
        """
        mu = 5
        Img_clean = l1cPy.breg_anistropic_TV(self.Img_noisy,
                                             max_iter=10000, max_jac_iter=10,
                                             tol=0.00001, mu=mu)

        Tvy_clean = norm(np.diff(Img_clean, axis=0).flatten(), ord=1)
        Tvy_dirty = norm(np.diff(self.Img_noisy, axis=0).flatten(), ord=1)

        Tvx_clean = norm(np.diff(Img_clean, axis=1).flatten(), ord=1)
        Tvx_dirty = norm(np.diff(self.Img_noisy, axis=1).flatten(), ord=1)

        self.assertLess(Tvy_clean, Tvy_dirty)
        self.assertLess(Tvx_clean, Tvx_dirty)

        assert_almost_equal(Img_clean, self.Img_orig, decimal=0.5)


if __name__ == '__main__':
    unittest.main()
