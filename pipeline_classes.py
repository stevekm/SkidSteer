#!/usr/bin/env python

from base_classes import *
from pipeline_functions import *

class printFoo(ShellTask):
    '''
    A class that invokes 'echo' from the shell
    '''
    def __init__(self):
        ShellTask.__init__(self, "printfoo")
        self.sys_command = '''
echo foo
        '''

class convertFastq2Txt(ShellTask):
    '''
    Convert a .fastq file to a .txt file
    '''
    def __init__(self):
        ShellTask.__init__(self, "txt2fastq")
    def build_command(self):
        '''
        Do the conversion
        '''
        command = 'echo "In this step I am going to convert a file"'
        return(command)


class runFastQC(ShellTask):
    '''
    A task that runs FastQC
    '''
    def __init__(self):
        ShellTask.__init__(self, 'fastqc')
        # super(runFastQC, self).__init__(**kwargs)
        self.environment = 'module load fastqc/0.11.4'
        self.sys_command ='fastqc'
        self.threads = 4
        self.params = '--nogroup'
        self.input_suffix = ".fastq"
    def build_command(self):
        '''
        Method for building the FastQC command
        ex: fastqc --threads "$THREADS" --nogroup --outdir /path/to/output file.fastq.gz
        '''
        command = '{}\n{}'.format(
        self.environment,
        self.sys_command
        )
        command = '{} --threads {} {} --outdir {}'.format(command, self.threads, self.params, self.output_dir)
        command = '{} {}'.format(command, ' '.join([file for file in self.input_files]))
        return(command)

class doPythonCode(PythonTask):
    '''
    An example of a class that executes some custom Python code
    '''
    def __init__(self):
        PythonTask.__init__(self, 'python-task')
    def build_command(self):
        print("This is some Python code")


class convertFile(PythonTask):
    '''
    Write lines from input file to output file in a reversed order
    '''
    def __init__(self):
        PythonTask.__init__(self, 'convert-file')
    def build_command(self):
        for file in self.input_file:
            with open(file) as input_file, open('newtest', 'w') as new:
                for line in old:
                    if line.rsplit('|', 1)[-1].strip() == 'number3':
                        new.write('this is replacement|number7\n')
                    else:
                        new.write(line)
