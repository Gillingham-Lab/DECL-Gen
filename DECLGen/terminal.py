import shutil
from colorama import init, deinit, Fore, Back, Style
from typing import Tuple, List, Sequence, Any, Callable
from .exceptions import DECLException

class Terminal:
    bold = Fore.CYAN
    underline = Fore.CYAN
    red = Fore.RED
    yellow = Fore.YELLOW
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

    def highlight(self, msg):
        return "{t.bold}{msg}{t.normal}".format(msg=msg, t=self)

    def val_good(self, msg):
        return "{t.green}{msg}{t.normal}".format(msg=msg, t=self)

    def val_warning(self, msg):
        return "{t.yellow}{msg}{t.normal}".format(msg=msg, t=self)

    def val_bad(self, msg):
        return "{t.red}{msg}{t.normal}".format(msg=msg, t=self)

    def table(self, shape: Sequence[int], first_column: bool = False, first_row: bool = False):
        """ Returns a TerminalTable object to easily display a table. """
        return TerminalTable(
            self,
            shape,
            first_column=first_column,
            first_row=first_row
        )

    def dl(
        self,
        shape_key: int,
        shape_val: int,
        highlight_key: bool = False,
        list_item: str = "-"
    ):
        """ Returns a TerminalDefinitionList object to quickly make a definition list. """
        return TerminalDefinitionList(
            self,
            shape_key,
            shape_val,
            highlight_key=highlight_key,
            list_item=list_item
        )


class TerminalTable:
    t: Terminal = None
    shape: Tuple[int] = None
    first_column: bool = None
    first_row: bool = None
    rows: List = []

    def __init__(self,
         t: Terminal,
         shape: Tuple[int],
         first_column: bool = False,
         first_row: bool = False
    ):
        self.t = t
        self.shape = shape
        self.first_column = first_column
        self.first_row = first_row
        self.rows = []

    def add_row(self, *data):
        if len(data) != len(self.shape):
            raise DECLException("TerminalTable row must have same number of columns as given in shape.")

        self.rows.append(data)

    def display(self):
        row_normal = ["{{:<{:d}}}".format(x) for x in self.shape]
        row_first = [x for x in row_normal]
        row_first[0] = "{t.bold}" + row_first[0] + "{t.normal}"

        row_normal = " ".join(row_normal)
        row_first = " ".join(row_first)

        for x in range(len(self.rows)):
            if self.first_row is True and x == 0:
                tmp = "{t.underline}" + row_normal + "{t.normal}"
            elif self.first_column is True:
                tmp = row_first
            else:
                tmp = row_normal

            print(tmp.format(*self.rows[x], t=self.t))

class TerminalDefinitionList:
    t: Terminal = None
    shape_key: int = None
    shape_val: int = None
    highlight_key: bool = None
    list_item: str = None
    t: List[Tuple[str, Any, Callable]] = None

    def __init__(self,
         t: Terminal,
         shape_key: int, shape_val: int,
         highlight_key: bool = False,
         list_item: str = "-"
     ):
        self.t = t
        self.shape_key = shape_key
        self.shape_val = shape_val
        self.highlight_key = highlight_key
        self.rows = []
        self.list_item = list_item

    def add_row(self, key: str, val: Any, format = None):
        if format is None:
            format = lambda x: x
        self.rows.append((key, val, format))

    def display(self):
        row_normal = [
            "{{key:<{:d}}}".format(self.shape_key),
            "{{val:>{:d}}}".format(self.shape_val),
        ]
        row_highlight = [x for x in row_normal]
        row_highlight[0] = "{t.bold}" + row_highlight[0] + "{t.normal}"


        if self.highlight_key:
            tmp = row_highlight
        else:
            tmp = row_normal

        for row in self.rows:
            print(" {} ".format(self.list_item), end="")
            print(tmp[0].format(key=row[0], t=self.t), end="")
            print(row[2](tmp[1].format(val=row[1], t=self.t)))