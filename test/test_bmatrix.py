import unittest

import KGBN


class BNLoadAndUpdateTest(unittest.TestCase):
    def test_load_bn_from_string_and_update(self):
        network = KGBN.load_network_from_string(
            """
            A = A
            B = A
            C = !B
            """,
            initial_state={"A": 1, "B": 0, "C": 0},
        )

        self.assertEqual(network.N, 3)
        self.assertEqual(network.nodeDict, {"A": 0, "B": 1, "C": 2})

        results = network.update(iterations=2)
        self.assertEqual(results[1, :].tolist(), [1, 1, 0])
        self.assertEqual(results[2, :].tolist(), [1, 1, 0])

    def test_load_network_auto_detects_bn(self):
        network = KGBN.load_network(
            """
            A = A
            B = A
            """,
            initial_state=[1, 0],
        )

        self.assertEqual(network.N, 2)
        self.assertEqual(network.update(iterations=1)[1, :].tolist(), [1, 1])


if __name__ == "__main__":
    unittest.main()
