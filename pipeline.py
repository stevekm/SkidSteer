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
    def __init__(self, name):
        '''
        Initialize the object with default attributes
        These can be overwritten in your Task subclasses
        '''
        # attributes that build the system command
        self.environment = '' # "module load goes here"
        self.hpc_params = '' # submit_job_to_cluster -b y -j y -pe threaded 4
        self.pre_processing = ''
        self.sys_command = ''
        self.post_processing = ''
        # base location & settings
        self.task_name = name
        self.input_file = None
        self.output_dir_base = 'output'
        self.output_dir = os.path.join(self.output_dir_base,self.task_name)
    def build_command(self):
        '''
        Default method to build the system command to be passed to the shell
        This can be overwritten in your Task subclasses
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
        Task.__init__(self, "printfoo")
        self.sys_command = '''
echo foo
        '''
class runFastQC(Task):
    def __init__(self):
        Task.__init__(self, 'fastqc')
        self.environment = 'module load fastqc/0.11.4'
        self.sys_command ='fastqc'
        self.threads = 4
        self.params = '--nogroup'
    def build_command(self):
        '''
        Method for building the FastQC command
        ex: fastqc --threads "$THREADS" --nogroup --outdir "$tmpFastQCdir1" "$INPUTFILE"
        '''
        command = '{}\n{}'.format(
        self.environment,
        self.sys_command
        )
        command = '{} --threads {} {} --outdir {}'.format(command, self.threads, self.params, self.output_dir)
        command = '{} {}'.format(command, ' '.join([file for file in self.input_file]))
        return(command)

def get_sample_files(sampleID, dir, suffix):
    '''
    Find files for a given sample in the supplied directory
    '''
    # do things
    print()

def run_pipeline(task_list):
    '''
    Run all tasks supplied to the pipeline
    '''
    import time
    for task in task_list:
        print("--------------------")
        print("Now running task: {}".format(task.task_name))
        print("Command is:\n{}".format(task))
        task.run()
        time.sleep(3) # delays for 3 seconds

# create task objects
foo = printFoo()

fastqc = runFastQC()
fastqc.input_file = ["input/test_R1.fastq"]

# add task objects to task list
task_list = [
foo,
fastqc
]

# run the pipeline on the task objects
run_pipeline(task_list)
# my_debugger(globals().copy())
