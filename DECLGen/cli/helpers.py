import math
import time


class ProgressBar:
    term = None
    location = None
    width = None
    desc: str = None

    def __init__(self, terminal, location = None, desc=None):
        self.term = terminal

        if location is None:
            location = (0, terminal.height - 1)

        self.location = location
        self.width = terminal.width
        self.desc = desc

    def write(self, string):
        print(string, end='')

    def start(self):
        self.bar_length = self.width - 10

        if self.desc is None:
            self.write("|>{}|{:>3}%\r".format(" "*self.bar_length, "0"))
        else:
            self.bar_length -= (len(self.desc) + 2)
            self.write("{}: |>{}|{:>3}%\r".format(self.desc, " "*self.bar_length, "0"))


    def update(self, percentage):
        bar_length = math.floor(percentage * self.bar_length)
        white_length = self.bar_length - bar_length

        if self.desc is None:
            self.write("|{}>{}|{: >3}%\r".format("="*bar_length, " "*white_length, math.floor(percentage*100)))
        else:
            self.write("{}: |{}>{}|{: >3}%\r".format(self.desc, "="*bar_length, " "*white_length, math.floor(percentage*100)))


    def finish(self):
        if self.desc is None:
            self.write("|{}| Done\n".format("="*(self.bar_length+1)))
        else:
            self.write("{}: |{}| Done\n".format(self.desc, "="*(self.bar_length+1)))

    def __enter__(self):
        self.start()

    def __exit__(self,  type, value, tb):
        self.finish()

