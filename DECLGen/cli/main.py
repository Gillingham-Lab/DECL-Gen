import argh
from DECLGen.cli.commands import about, commands


def main():
    # Make commands available for argh to generate cli interface
    p = argh.ArghParser()
    p.add_commands(sorted(commands, key=lambda x: x.__name__))

    p.set_default_command(about)

    argh.dispatch(p)
