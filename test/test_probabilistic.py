import unittest

import numpy as np

import KGBN


class PBNLoadAndUpdateTest(unittest.TestCase):
    def test_load_pbn_from_string(self):
        np.random.seed(1)
        network = KGBN.load_pbn_from_string(
            """
            A = A, 1.0
            B = A, 0.7
            B = !A, 0.3
            """,
            initial_state={"A": 1, "B": 0},
        )

        self.assertEqual(network.N, 2)
        self.assertEqual(network.nf.tolist(), [1, 2])
        self.assertTrue(np.allclose(network.cij[1, :2], [0.7, 0.3]))

        trajectory = network.update(iterations=3)
        self.assertEqual(trajectory.shape, (4, 2))

    def test_load_network_auto_detects_pbn(self):
        network = KGBN.load_network(
            """
            A = A, 1
            B = A, 0.5
            B = !A, 0.5
            """,
            initial_state=[1, 0],
        )

        self.assertEqual(network.nf.tolist(), [1, 2])


if __name__ == "__main__":
    unittest.main()
