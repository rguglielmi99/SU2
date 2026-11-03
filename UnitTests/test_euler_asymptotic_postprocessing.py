#!/usr/bin/env python3
"""Scientific sanity checks independent of any measured convergence outcome."""
import unittest
import numpy as np
from export_euler_asymptotic_study import conservative, primitive, masks, norms


class Diagnostics(unittest.TestCase):
    def test_richardson_cancels_only_leading_temporal_error(self):
        exact = np.array([1., .5, .3, 1.])
        def field(dt):
            return exact+dt*np.array([.1, -.3, .2, .4])+dt**2*np.array([.2, .1, -.1, .3])
        coarse = 2*field(.05)-field(.1)
        fine = 2*field(.025)-field(.05)
        np.testing.assert_allclose(coarse-exact, 4*(fine-exact), rtol=1e-11, atol=1e-15)

    def test_conservative_extrapolation_precedes_conversion(self):
        w1, w2 = np.array([[1.1, .4, .3, 1.2]]), np.array([[1.05, .5, .35, 1.1]])
        correct = primitive(2*conservative(w2)-conservative(w1))
        incorrect = 2*w2-w1
        self.assertGreater(np.max(np.abs(correct-incorrect)), 1e-4)
        np.testing.assert_allclose(primitive(conservative(w1)), w1, atol=1e-14)

    def test_region_energy_is_additive_despite_nonuniform_volumes(self):
        layer = np.array([0, 1, 2, 3, 4])
        region = masks(np.column_stack((np.linspace(0, 1, 5), np.linspace(0, 1, 5))), layer)
        count = sum(region[k].astype(int) for k in ('boundary', 'ring1', 'ring2', 'deep3'))
        np.testing.assert_array_equal(count, np.ones(5))
        error = np.arange(20).reshape(5, 4)*.01
        v = np.array([.02, .03, .04, .05, .06])
        full, energy = norms(error, v)
        self.assertFalse(np.allclose(full['rms'], full['l2']))
        pieces = sum(norms(error[region[k]], v[region[k]])[1] for k in ('boundary', 'ring1', 'ring2', 'deep3'))
        np.testing.assert_allclose(pieces, energy, atol=1e-15)


if __name__ == '__main__':
    unittest.main()
