import math
import time
from blessings import Terminal


class ProgressBar:
    term = None
    location = None
    width = None

    def __init__(self, terminal: Terminal, location = None):
        self.term = terminal
        if location is None:
            location = (0, terminal.height - 1)
        self.location = location
        self.width = terminal.width

    def write(self, string):
        print(string, end='')

    def start(self):
        self.bar_length = self.width - 10
        self.write("|>{}|{:>3}%\r".format(" "*self.bar_length, "0"))

    def update(self, percentage):
        bar_length = math.floor(percentage * self.bar_length)
        white_length = self.bar_length - bar_length
        self.write("|{}>{}|{: >3}%\r".format("="*bar_length, " "*white_length, math.floor(percentage*100)))

    def finish(self):
        self.write("|{}| Done\n".format("="*(self.bar_length+1)))



