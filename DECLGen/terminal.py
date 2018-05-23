import shutil
from colorama import init, Fore, Back, Style

class Terminal:
    bold = Fore.CYAN
    underline = Fore.CYAN
    red = Fore.RED
    normal = Style.RESET_ALL

    def __init__(self):
        init(autoreset=False, convert=True)

    @property
    def height(self):
        return shutil.get_terminal_size()[1]

    @property
    def width(self):
        return shutil.get_terminal_size()[0]