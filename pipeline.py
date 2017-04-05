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
    print(build_io_files(input_dir = "input", input_suffix = ".fastq", output_dir = "output", output_suffix = ".txt", sampleID = sampleID))
    # ~~~~~~ CREATE TASK OBJECTS ~~~~~ #
    foo = printFoo(sampleID = sampleID)

    fastqc = runFastQC(input_dir = 'input', sampleID = sampleID) #
    align_bowtie2 = alignBowtie2(input_dir = 'input', sampleID = sampleID)

    # my_debugger(globals().copy())
    do_python = doPythonCode(sampleID = sampleID)

    convert_fastq = convertFastq2Txt(sampleID = sampleID)

    # ~~~~~~ BUILD TASK LIST ~~~~~ #
    task_list = [
    foo,
    fastqc,
    align_bowtie2,
    do_python,
    convert_fastq
    ]

    # ~~~~~~ RUN THE TASKS ~~~~~ #
    run_pipeline(task_list, sampleID)
