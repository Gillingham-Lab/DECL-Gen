import argh
from DECLGen.cli.gen_commands import about, commands


def main():
    # Make gen_commands available for argh to generate cli interface
    p = argh.ArghParser()
    p.add_commands(commands)

    p.set_default_command(about)

    argh.dispatch(p)
