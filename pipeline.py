#!/usr/bin/env python

import sys
import os

def subprocess_cmd(command):
    # run a terminal command with stdout piping enabled
    import subprocess as sp
    process = sp.Popen(command,stdout=sp.PIPE, shell=True)
    proc_stdout = process.communicate()[0].strip()
    print(proc_stdout)

def my_debugger(vars):
    # starts interactive Python terminal at location in script
    # call with pl.my_debugger(globals().copy()) anywhere in your script
    # or call my_debugger(locals().copy()) from anywhere within this package
    import readline # optional, will allow Up/Down/History in the console
    import code
    # vars = globals().copy() # in python "global" variables are actually module-level
    vars.update(locals())
    shell = code.InteractiveConsole(vars)
    shell.interact()

# debug = True
# input_files = ["foo.txt", "bar.txt"]
#
# environment_setup = """
# echo "module load goes here"
# """
#
# task_command = """
# {}
# echo '{}'
# """.format(environment_setup, ' '.join(input_files))
#
# if debug == True:
#     task_command = '''
# set -x
# {}
#     '''.format(task_command)
#
# print("Command is:",task_command)
#
# subprocess_cmd(task_command)

class Task():
    '''
    A base class for pipeline tasks
    '''
    def __init__(self):
        self.environment = 'echo "module load goes here"'

        self.hpc_params='submit_job_to_cluster -b y -j y -pe threaded 4'

        self.command = '{}'.format('\n'.join([self.environment, self.hpc_params]))
    def __repr__(self):
        return self.command

class PrintFoo(Task):
    def __init__(self):
        Task.__init__(self)
        self.command = '{} echo "foooo"'.format(self.command)

x = PrintFoo()
my_debugger(globals().copy())
