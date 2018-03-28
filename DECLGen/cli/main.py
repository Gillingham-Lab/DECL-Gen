#!/usr/bin/env python
import argh

def main():
    # Prepare commands
    commands = [

    ]

    # Make commands available for argh to generate cli interface
    p = argh.ArghParser()
    p.add_commands(sorted(commands, key=lambda x: x.__name__))

    p.set_default_command(about)