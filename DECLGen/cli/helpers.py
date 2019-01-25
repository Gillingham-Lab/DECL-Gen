import math
import time
from DECLGen.terminal import Terminal


class ProgressBar:
    """
    Shows a progress bar to indicate loading. This class can both be used directly and within a context operator.
    """
    term: Terminal = None
    location = None
    width = None
    desc: str = None

    def __init__(self, terminal, location = None, desc=None):
        """

        :param terminal: A reference to the active DECLGen.terminal.Terminal. This is used to determine terminal dimensions.
        :param location: Currently not really supported. Dummy parameter.
        :param desc: Adds an optional description in front of the progress bar.
        """
        self.term = terminal

        if location is None:
            location = (0, terminal.height - 1)

        self.location = location
        self.width = terminal.width
        self.desc = desc

    def _write(self, string):
        """
        Writes the progress bar. Internal use only.
        :param string:
        :return:
        """
        print(string, end='')

    def start(self) -> None:
        """
        Starts the progress bar
        :return:
        """
        self.bar_length = self.width - 10

        if self.desc is None:
            self._write("▓▒{} {:>3}%\r".format(" " * self.bar_length, "0"))
        else:
            self.bar_length -= (len(self.desc) + 2)
            self._write("{}: ▓▒{} {:>3}%\r".format(self.desc, "░" * self.bar_length, "0"))


    def update(self, percentage) -> None:
        """
        Updates the progress bar to a certain percentage. Cannot ever be set above 100%.
        :param percentage: The new percentage the progress bar should display
        :return: None
        """
        percentage = min(1, percentage)
        bar_length = math.floor(percentage * self.bar_length)
        white_length = self.bar_length - bar_length

        if self.desc is None:
            self._write("▓{}▒{} {: >3}%\r".format("▓" * bar_length, "░" * white_length, math.floor(percentage * 100)))
        else:
            self._write("{}: ▓{}▒{} {: >3}%\r".format(self.desc, "▓" * bar_length, "░" * white_length, math.floor(percentage * 100)))


    def finish(self) -> None:
        """
        Finished the progress bar by displaying a full bar and "Done" instead of the percentage.
        :return:
        """
        if self.desc is None:
            self._write("▓{}  Done\n".format("▓" * (self.bar_length + 1)))
        else:
            self._write("{}: ▓{}  Done\n".format(self.desc, "▓" * (self.bar_length + 1)))

    def __enter__(self):
        self.start()
        return self

    def __exit__(self,  type, value, tb):
        self.finish()

