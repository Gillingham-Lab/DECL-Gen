import argh
from DECLGen.cli.tool_commands import commands
from DECLGen.cli import about


def main():
    # Make gen_commands available for argh to generate cli interface
    p = argh.ArghParser()
    p.add_commands([about] + sorted(commands, key=lambda x: x.__name__))

    p.set_default_command(about)

    argh.dispatch(p)
