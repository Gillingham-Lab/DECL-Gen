import shutil
from colorama import init, deinit, Fore, Back, Style
from typing import Tuple, List, Sequence, Any, Callable
from .exceptions import DECLException


class Terminal:
    """
    Terminal helper class to standardize output in DECL gen.
    """
    bold = Fore.CYAN
    underline = Fore.CYAN
    red = Fore.RED
    yellow = Fore.YELLOW
    normal = Style.RESET_ALL
    green = Fore.GREEN

    def __init__(self):
        # Initializes colorama.
        init()

    def __del__(self):
        # Uninitializes colorama.
        deinit()

    @property
    def height(self) -> int:
        """
        Property to return the height of the terminal.
        :return: Height of the terminal in characters.
        """
        return shutil.get_terminal_size()[1]

    @property
    def width(self) -> int:
        """
        Property to return the width of the terminal.
        :return: Width of the terminal in characters.
        """
        return shutil.get_terminal_size()[0]

    @property
    def shape(self) -> Tuple[int, int]:
        """
        Property to return the shape of the terminal.
        :return: Tuple (width, height) in characters.
        """
        return self.width, self.height

    def highlight(self, msg):
        """
        Highlights a message.
        :param msg: Input message.
        :return: Output message.
        """
        return "{t.bold}{msg}{t.normal}".format(msg=msg, t=self)

    def val_good(self, msg):
        """
        Formats a message to be good.
        :param msg: Input message.
        :return: Output message.
        """
        return "{t.green}{msg}{t.normal}".format(msg=msg, t=self)

    def val_warning(self, msg):
        """
        Formats a message as a warning.
        :param msg: Input message.
        :return: Output message.
        """
        return "{t.yellow}{msg}{t.normal}".format(msg=msg, t=self)

    def val_bad(self, msg):
        """
        Formats a message as an error.
        :param msg: Input message.
        :return: Output message.
        """
        return "{t.red}{msg}{t.normal}".format(msg=msg, t=self)

    def decide(self, msg: str, question: str, answers: Sequence[str]):
        """
        Asks for a decision between different options. Re-Asks if the user does not give an answer.
        :param msg: A message or explanation to give.
        :param question: The question to ask.
        :param answers: Possible answers.
        :return: The given answer.
        """
        print("{t.yellow}{msg}{t.normal}".format(msg=msg, t=self))

        answers_string = "/".join(answers)

        while True:
            answer = input("{question} ({answers}) ".format(question=question, answers=answers_string)).strip()

            if answer in answers_string:
                return answer

            if answer[0] in answers_string:
                return answer


    def table(self, shape: Sequence[int], first_column: bool = False, first_row: bool = False):
        """
        Generates a TerminalTable object with a specific shape and options.
        :param shape: A sequence of how logn each column should be.
        :param first_column: True of first column should be highlighted.
        :param first_row: True of first row should be highlighted.
        :return: TerminalTable object.
        """
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
        list_item: str = "-",
        values_justify_right = False,
    ):
        """
        Returns a TerminalDefinitionList object that allows to quickly make a definition list.
        :param shape_key: Width of the key in characters.
        :param shape_val: Width of the value in characters.
        :param highlight_key: True if the key should be highlighted.
        :param list_item: A character or character combination to symbolize each list entry.
        :param values_justify_right: True if values should be justified on the right.
        :return: TerminalDefinitionList
        """
        return TerminalDefinitionList(
            self,
            shape_key,
            shape_val,
            highlight_key=highlight_key,
            list_item=list_item,
            values_justify_right=values_justify_right,
        )


class TerminalTable:
    """
    An object for generating a table.
    """
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
        """
        Constructor.
        :param t: Reference for the terminal instance.
        :param shape: List of integers describing the width of each column.
        :param first_column: True if the first column should be highlighted.
        :param first_row: False if the first column should be highlighted.
        """
        self.t = t
        self.shape = shape
        self.first_column = first_column
        self.first_row = first_row
        self.rows = []

    def add_row(self, *data):
        """
        Adds a row to the table. Must have the same number of columns as defined in shape.
        :param data:
        :return:
        """
        if len(data) != len(self.shape):
            raise DECLException("TerminalTable row must have same number of columns as given in shape.")

        self.rows.append(data)

    def display(self):
        """
        Renders the table and prints it.
        :return: None
        """
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
    """
    A definition list to display key: value relations as a list.
    """
    t: Terminal = None
    shape_key: int = None
    shape_val: int = None
    highlight_key: bool = None
    list_item: str = None
    rows: List[Tuple[str, Any, Callable]] = None
    values_justify_right: bool = None

    def __init__(
            self,
            t: Terminal,
            shape_key: int, shape_val: int,
            highlight_key: bool = False,
            list_item: str = "-",
            values_justify_right = False,
     ):
        """
        Constructor.
        :param t:  Reference to the terminal instance.
        :param shape_key: Width of the key in characters.
        :param shape_val: Width of the value in characters.
        :param highlight_key: True if the key should be highlighted.
        :param list_item: A character or character combination to symbolize each list entry.
        :param values_justify_right: True if values should be justified on the right.
        """
        self.t = t
        self.shape_key = shape_key
        self.shape_val = shape_val
        self.highlight_key = highlight_key
        self.rows = []
        self.list_item = list_item
        self.values_justify_right = values_justify_right

    def add_row(self, key: str, val: Any, format = None):
        """
        Adds a row to the definition list.
        :param key: Row title
        :param val: Row value
        :param format: Formatting string used in .format()
        :return:
        """
        if format is None:
            format = lambda x: x
        self.rows.append((key, val, format))

    def display(self):
        """
        Renders the definition list and prints the output.
        :return:
        """
        row_normal = [
            "{{key:<{:d}}}".format(self.shape_key),
            "{{val:{:}{:d}}}".format(
                ">" if self.values_justify_right is True else "<",
                self.shape_val
            ),
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