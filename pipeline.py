#!/usr/bin/env python

import sys
import os

def subprocess_cmd(command):
    '''
    run a terminal command with stdout piping enabled
    '''
    import subprocess as sp
    process = sp.Popen(command,stdout=sp.PIPE, shell=True)
    proc_stdout = process.communicate()[0].strip()
    print(proc_stdout)

def my_debugger(vars):
    '''
    starts interactive Python terminal at location in script
    call with pl.my_debugger(globals().copy()) anywhere in your script
    or call my_debugger(locals().copy()) from anywhere within this package
    '''
    import readline # optional, will allow Up/Down/History in the console
    import code
    # vars = globals().copy() # in python "global" variables are actually module-level
    vars.update(locals())
    shell = code.InteractiveConsole(vars)
    shell.interact()

def mkdirs(path, return_path=False):
    '''
    make a directory, and all parent dir's in the path
    '''
    import sys
    import os
    import errno

    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise
    if return_path:
        return path


class Task():
    '''
    A base class for pipeline tasks
    '''
    def __init__(self):
        '''
        Initialize the object with default attributes
        '''
        # attributes that build the system command
        self.environment = '' # "module load goes here"
        self.hpc_params = '' # submit_job_to_cluster -b y -j y -pe threaded 4
        self.pre_processing = ''
        self.sys_command = ''
        self.post_processing = ''
        # base locations
        self.task_name = ''
        self.input = ''
        self.output_dir_base = 'output'
        self.output_dir = os.path.join(self.output_dir_base,self.task_name)
    def build_command(self):
        '''
        Build the system command to pass to subprocess & run in the shell
        '''
        command = '{}'.format('\n'.join([
        self.environment,
        self.hpc_params,
        self.pre_processing,
        self.sys_command,
        self.post_processing
        ]))
        return(command)
    def run(self):
        '''
        Run the system command in the shell
        '''
        command = self.build_command()
        subprocess_cmd(command)
    def __repr__(self):
        return self.build_command()

class printFoo(Task):
    def __init__(self):
        Task.__init__(self)
        self.sys_command = '''
        echo foo
        '''
class runFastQC(Task):
    def __init__(self):
        Task.__init__(self)
        self.environment = '''
        module load fastqc/0.11.4
        '''
        self.sys_command ='''
        fastqc --version
        '''

def run_pipeline(task_list):
    for task in task_list:
        task.run()


foo = printFoo()
fastqc = runFastQC()
task_list = [
foo,
fastqc
]

my_debugger(globals().copy())

run_pipeline(task_list)
