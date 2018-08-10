import shutil
from colorama import init, deinit, Fore, Back, Style
from typing import Tuple

class Terminal:
    bold = Fore.CYAN
    underline = Fore.CYAN
    red = Fore.RED
    normal = Style.RESET_ALL
    green = Fore.GREEN

    def __init__(self):
        init()

    def __del__(self):
        deinit()

    @property
    def height(self) -> int:
        return shutil.get_terminal_size()[1]

    @property
    def width(self) -> int:
        return shutil.get_terminal_size()[0]

    @property
    def shape(self) -> Tuple[int, int]:
        return self.width, self.height