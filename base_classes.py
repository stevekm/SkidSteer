#!/usr/bin/env python

class ShellTask():
    '''
    A base class for pipeline tasks that will be run in the system shell
    '''
    def __init__(self, name, sampleID = None):
        '''
        Initialize the object with default attributes
        These can be overwritten in your Task subclasses
        '''
        import sys
        import os
        self.sampleID = sampleID
        # attributes that build the system command
        self.environment = ''
        self.hpc_params = ''
        self.pre_processing = ''
        self.sys_command = ''
        self.post_processing = ''
        # base location & settings
        self.task_name = name
        self.input_files = None
        self.input_dir = None
        self.input_suffix = None
        self.inputs = (None, None) # ("dir", ".suffix")
        self.input_from_task = None
        self.output_dir_base = 'output'
        self.output_dir = os.path.join(self.output_dir_base,self.task_name)
    def set_inputs(self):
        '''
        Set the Input files and input dir for the object
        '''
        # do a thing

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
            print("Now running task: {}".format(self.task_name))
            print("Command is:\n{}".format(self))
            subprocess_cmd(command)
    def __repr__(self):
        return self.build_command()

class PythonTask():
    '''
    A base class for tasks which will be run using custom internal Python code
    '''
    def __init__(self, name, sampleID = None):
        '''
        Initialize the object with default attributes
        These can be overwritten in your Task subclasses
        '''
        import sys
        import os
        self.task_name = name
        self.input_files = None
        self.input_from_task = None
        self.output_dir_base = 'output'
        self.output_dir = os.path.join(self.output_dir_base,self.task_name)
    def build_command(self):
        '''
        Default method to setup the Python code to execute
        This can be overwritten in your Task subclasses
        '''
        print("No commands have been set for task {}".format(self.task_name))
    def run(self):
        '''
        Run the custom Python code that has been set
        '''
        print("--------------------")
        print("Now running task: {}".format(self.task_name))
        self.build_command()
    def __repr__(self):
        return self.task_name
