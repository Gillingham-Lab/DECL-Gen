import unittest

from DECLGen import Terminal

class TerminalTestCase(unittest.TestCase):
    def test_terminal_colours(self):
        self.assertIsNotNone(Terminal.bold)
        self.assertIsNotNone(Terminal.underline)
        self.assertIsNotNone(Terminal.red)
        self.assertIsNotNone(Terminal.normal)
        self.assertIsNotNone(Terminal.green)

    def test_terminal_shape(self):
        t = Terminal()

        self.assertIsInstance(t.width, int)
        self.assertGreaterEqual(t.width, 0)

        self.assertIsInstance(t.height, int)
        self.assertGreaterEqual(t.height, 0)

        self.assertTupleEqual((t.width, t.height), t.shape)