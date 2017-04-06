#!/usr/bin/env python

def subprocess_cmd(command):
    '''
    Run a terminal command with stdout piping enabled
    '''
    import subprocess as sp
    process = sp.Popen(command,stdout=sp.PIPE, shell=True)
    proc_stdout = process.communicate()[0].strip()
    print(proc_stdout)

def my_debugger(vars):
    '''
    Starts interactive Python terminal at location in script
    call with my_debugger(globals().copy()) anywhere in your script
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
    Make a directory, and all parent dir's in the path
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


def find_sample_files(dir, sampleID, file_suffix):
    import os
    file_list = []
    for file in os.listdir("input"):
        if file.endswith(file_suffix):
            if file.startswith(sampleID):
                file_list.append(os.path.join(dir, file))
    return(file_list)


def run_pipeline(task_list, sampleID):
    '''
    Run all tasks supplied to the pipeline
    '''
    import time
    for task in task_list:
        # make the outdir if present
        if task.output_dir != None:
            mkdirs(task.output_dir)
        task.run()
        time.sleep(3) # delays for 3 seconds


def transform_file_suffix(input_file, input_suffix, output_suffix):
    '''
    Change the file extension on a file
    '''
    import os
    if input_file.endswith(input_suffix):
        filename, file_extension = os.path.splitext(input_file)
        output_filename = filename + output_suffix
        return(output_filename)

# def transform_file(input_dir, input_suffix, output_suffix, sampleID):
#     '''
#     Return the output file for a given sample input file
#     using glob
#     '''
#     import glob
#     import os
#     # get the input files
#     glob_pattern = os.path.join(input_dir, input_suffix)
#     input_file_list = []
#     for item in glob.glob(glob_pattern):
#         input_file_list.append(item)

def build_io_files(input_dir, input_suffix, output_dir, output_suffix, sampleID):
    '''
    Return a list of the output files to be created for a given set of input criteria
    '''
    # print(build_io_files(input_dir = "input", input_suffix = ".fastq", output_dir = "output", output_suffix = ".txt", sampleID = sampleID))
    # get input files from input dir
    input_files = find_sample_files(dir = input_dir, sampleID = sampleID, file_suffix = input_suffix)
    # get output file paths
    output_files = []
    for file in input_files:
        output_files.append(transform_file_suffix(input_file = file, input_suffix = input_suffix, output_suffix = output_suffix))
    io_dict = {"input": input_files, "output": output_files}
    return(io_dict)

def set_input_from_task(self, task, suffix):
    '''
    Set the input dir files from another task's output dir
    '''
    self.input_dir = task.output_dir
    self.input_suffix = suffix
    print('')
