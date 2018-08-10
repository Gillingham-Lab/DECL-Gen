import unittest
from random import randint

from DECLGen.evaluation import AlignmentResult

class ResultTestCase(unittest.TestCase):
    attributes = [
        "reads_useful",
        "valid_pairs",
        "invalid_pairs",
        "low_quality_skips",
        "both_low_quality_skips",
        "r1_low_quality_skips",
        "r2_low_quality_skips"
    ]

    def test_if_freshly_created_result_is_set_to_0(self):
        result = AlignmentResult()

        for a in self.attributes:
            with self.subTest(i=a):
                self.assertEqual(0, result[a])

    def test_if_assigning_number_to_result_saves_and_returns_that_number(self):
        result = AlignmentResult()

        for a in self.attributes:
            with self.subTest(i=a):
                i = randint(0, 2e10)
                result[a] = i

                self.assertEqual(i, result[a])

    def test_if_addition_of_results_produces_correct_results(self):
        result1 = AlignmentResult()
        result2 = AlignmentResult()
        solution = {}

        for a in self.attributes:
            i1 = randint(0, 2e10)
            i2 = randint(0, 2e10)

            result1[a] = i1
            result2[a] = i2
            solution[a] = i1 + i2

        combined1 = result1 + result2
        combined2 = result2 + result1

        for a in self.attributes:
            with self.subTest(i=a):
                self.assertEqual(solution[a], combined1[a])
                self.assertEqual(solution[a], combined2[a])
