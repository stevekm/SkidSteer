#!/usr/bin/env python

from pipeline_functions import *
import sys
import os

class PipelineTask(object):
    '''
    A base class for all pipeline tasks.
    This holds attributes and methods that all pipeline tasks share
    '''
    def __init__(self, name, sampleID = None, input_task = None, input_dir = None, input_suffix = None, output_suffix = None, *args, **kwargs):
        self.name = name
        self.sampleID = sampleID
        # Set input files and dirs
        if input_suffix != None: # no input suffix = no input files
            if input_task != None:
                self.set_input_from_task(input_task, input_suffix)
            else:
                self.set_inputs(input_dir, input_suffix)
        # set output files and dirs
        self.output_dir_base = 'output'
        self.output_dir = os.path.join(self.output_dir_base, self.name)
        self.output_suffix = output_suffix
    def set_input_from_task(self, task, suffix):
        '''
        Set the input dir files from another task's output dir
        '''
        if issubclass(type(task), PipelineTask):
            self.input_dir = task.output_dir
            self.input_suffix = suffix
        else:
            print('ERROR: Task is not a sublcass of type PipelineTask')
    def set_inputs(self, input_dir, suffix):
        '''
        Set the Input files and input dir for the object
        '''
        self.input_dir = input_dir
        self.input_suffix = suffix
        # my_debugger(locals().copy())
        if self.input_dir != None:
            self.input_files = find_sample_files(self.input_dir, self.sampleID, self.input_suffix)
        else:
            self.input_files = []


class ShellTask(PipelineTask):
    '''
    A base class for pipeline tasks that will be run in the system shell
    '''
    def __init__(self, name, *args, **kwargs):
        '''
        Initialize the object with default attributes
        These can be overwritten in your Task subclasses
        '''
        PipelineTask.__init__(self, name, *args, **kwargs)
        # attributes that build the system command
        self.environment = ''
        self.hpc_params = ''
        self.pre_processing = ''
        self.sys_command = ''
        self.post_processing = ''
    def build_command(self):
        '''
        Default method to build the system command to be passed to the shell
        This can be overwritten in your Task subclasses
        If you don't plan to use the command you can set
        command = False
        just make sure you end with
        return(command)
        or the pipeline will break
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
        Don't run in subprocess if command = False
        '''
        command = self.build_command()
        if command != False:
            print("--------------------")
            print("Now running task: {}".format(self.name))
            print("Command is:\n{}".format(self))
            subprocess_cmd(command)
    def __repr__(self):
        return self.build_command()

class PythonTask(PipelineTask):
    '''
    A base class for tasks which will be run using custom internal Python code
    '''
    def __init__(self, name, *args, **kwargs):
        '''
        Initialize the object with default attributes
        These can be overwritten in your Task subclasses
        '''
        PipelineTask.__init__(self, name, *args, **kwargs)
    def build_command(self):
        '''
        Default method to setup the Python code to execute
        This can be overwritten in your Task subclasses
        '''
        print("No commands have been set for task {}".format(self.name))
    def run(self):
        '''
        Run the custom Python code that has been set
        '''
        print("--------------------")
        print("Now running task: {}".format(self.name))
        self.build_command()
    def __repr__(self):
        return self.name
