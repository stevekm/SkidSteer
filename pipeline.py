#!/usr/bin/env python

# ~~~~~~ LOAD PACKAGES ~~~~~ #
# system packages
import sys
import os

# pipeline packages
from base_classes import *
from pipeline_classes import *
from pipeline_functions import *

if __name__ == "__main__":
    # ~~~~~~ CREATE TASK OBKECTS ~~~~~ #
    sampleID = "test2"
    # print(build_io_files_list(input_dir = "input", input_suffix = ".fastq", output_dir = "output", output_suffix = ".txt", sampleID = sampleID))


    # ~~~~~~ CREATE TASK OBJECTS ~~~~~ #
    foo = printFoo()

    fastqc = runFastQC(input_dir = 'input') #
    # fastqc.input_dir = 'input'

    do_python = doPythonCode()

    convert_fastq = convertFastq2Txt()

    # ~~~~~~ BUILD TASK LIST ~~~~~ #
    task_list = [
    foo,
    fastqc,
    do_python,
    convert_fastq
    ]

    # ~~~~~~ RUN THE TASKS ~~~~~ #
    # my_debugger(globals().copy())
    run_pipeline(task_list, sampleID)
