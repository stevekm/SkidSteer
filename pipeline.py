#!/usr/bin/env python

# system packages
import sys
import os

# pipeline packages
from base_classes import *
from pipeline_classes import *
from pipeline_functions import *


myfile = "/path/to/thing.txt"
transform_file_suffix(myfile, ".txt", ".foo")

# convert_fastq_to_txt = convertFastq2Txt()
# convert_fastq_to_txt.input_dir = 'input'
# convert_fastq_to_txt.input_suffix = '*.fastq'

sample_ID = "test2"
file_suffix = ".fastq"

print(build_io_files_list(input_dir = "input", input_suffix = ".fastq", output_dir = "output", output_suffix = ".txt", sample_ID = sample_ID))

sys.exit()
# my_debugger(globals().copy())

# create task objects
foo = printFoo()

fastqc = runFastQC()
fastqc.input_files = ["input/test_R1.fastq"]

do_python = doPythonCode()

# add task objects to task list
task_list = [
foo,
fastqc,
do_python
]

# run the pipeline on the task objects
run_pipeline(task_list)
# my_debugger(globals().copy())
