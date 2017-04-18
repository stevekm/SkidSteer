# SkidSteer 
`SkidSteer` is a Python class-based framework for data analysis workflows, designed with bioinformatic pipelines in mind and inspired by the [Ruffus](http://www.ruffus.org.uk/) [package](https://pypi.python.org/pypi/ruffus). 

<a href="https://clipartfest.com/download/eb8d08b9430bf5d8348d61abdcdd197c5887f51f" title="Image from clipartfest.com"><img src="https://img.clipartfest.com/f46f2fc8b1f194af8219ab976525008f_-skid-loader-a-small-skid-skid-steer-mulcher-clipart-logo_249-194.jpeg" width="249" alt="skid steer mulcher clipart logo" /></a>

# Classes

## Base Classes

These classes are defined in `base_classes.py`

### `PipelineTask`

This class represents a task to be performed by the pipeline. It has attributes and methods that set the task's input & output directories and files. 

### `ShellTask`

A `PipelineTask` which ultimately invokes a command to run on the shell, using a wrapper to the Python `subprocess` package. 

### `PythonTask`

A `PipelineTask` which can hold custom Python code. 

## User Defined Classes

In order to use `SkidSteer`, the user should create subclasses of `ShellTask` and `PythonTask` for their pipeline tasks. These should be held in the `pipeline_classes.py` file. See [pipeline_classes.py](https://github.com/stevekm/SkidSteer/blob/master/pipeline_classes.py) for examples.

# Methods

Some notable object methods include:

- `build_command`: All classes should have a `build_command` method, which is used to build a shell command to run (`ShellTask`) or should contain custom Python code to run (`PythonTask`). This method should be filled in by the user.

- `run`: All classes will come with the `run` method, which invokes `build_command` and runs the generated command, if applicable. The user should not modify this method. 

# Pipeline

The complete list of tasks to be run in the pipeline are held in `pipeline.py`. Task objects should be created in the `build_pipeline_task_list` function, and added to a `task_list` in the desired order of completion. The `run_pipeline` function will evaluate and run each task in the list. 

# Design Goals

The following goals were considered in the development of the `SkidSteer` package:

- modularity of adding and managing tasks
- linear pipeline task completion
- standardized input/output file mapping (as much as possible)
- room for user to implement own methods & techniques inside of tasks
- object oriented system
- linear object class inheritance (no diamond inheritance or multiple inheritance)
- Python 2.7/3+ compatible
- **the main pipeline script (e.g. `pipeline.py`) should be easy to read and follow; sequester complex program logic in external modules (e.g. `base-classes.py`, etc.)**

The following aspects were actively avoided in the creation of `SkidSteer`:

- complex branching of inheritance or pipeline tasks
- excessively deep inheritance; currently at 2 base levels + 1 user level, try not to go any deeper

# Example Output

The following output can been seen by running the provided example [pipeline](https://github.com/stevekm/SkidSteer/blob/4c8a9100e89f48b77e866a13865f82c206ce5523/pipeline.py)

```
$ ./pipeline.py
--------------------
Running Pipeline

Sample ID is: test2


--------------------
Now running task: printfoo

Command is:




echo foo


foo
--------------------
Now running task: fastqc

Command is:
module load fastqc/0.11.4
fastqc --threads 4 --nogroup --outdir output/fastqc input/test2.fastq
/bin/sh: module: command not found
Started analysis of test2.fastq
Analysis complete for test2.fastq
--------------------
Now running task: bowtie2

Command is:
bowtie2 --local -x "~/ref/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome" -q -U input/test2.fastq --threads 1 | samtools view -@ "1" -Sb1 - | samtools sort -m 2G -@ "1" - -o  "output/bowtie2/test2.bam"
samtools index "output/bowtie2/test2.bam"
10 reads; of these:
  10 (100.00%) were unpaired; of these:
    0 (0.00%) aligned 0 times
    9 (90.00%) aligned exactly 1 time
    1 (10.00%) aligned >1 times
100.00% overall alignment rate

--------------------
Now running task: bamToBedgraph

Command is:
export LC_COLLATE=C
bedtools genomecov -ibam "output/bowtie2/test2.bam" -g "data/hg19.chrom.sizes_full.txt" -bg > "output/bamToBedgraph/test2.bedgraph"
sort -k1,1 -k2,2n output/bamToBedgraph/test2.bedgraph > output/bamToBedgraph/tmp && /bin/mv output/bamToBedgraph/tmp output/bamToBedgraph/test2.bedgraph

--------------------
Now running task: bamToBed

Command is:
export LC_COLLATE=C
bamToBed -i output/bowtie2/test2.bam > output/bamToBed/test2.bed
sort -k1,1 -k2,2n output/bamToBed/test2.bed > output/bamToBed/tmp && /bin/mv output/bamToBed/tmp output/bamToBed/test2.bed

--------------------
Now running task: bedToBigBed

Command is:
bedToBigBed "output/bamToBed/test2.bed" "data/hg19.chrom.sizes_full.txt" "output/bedToBigBed/test2.bigbed"
pass1 - making usageList (8 chroms): 1 millis
pass2 - checking and writing primary data (10 records, 6 fields): 2 millis

--------------------
Now running task: bedgraphToBigWig

Command is:
bedGraphToBigWig "output/bamToBedgraph/test2.bedgraph" "data/hg19.chrom.sizes_full.txt" "output/bedgraphToBigWig/test2.bw"

--------------------
Now running task: python-task

This is some Python code
--------------------
Now running task: txt2fastq

Command is:
echo "In this step I am going to convert a file"
In this step I am going to convert a file
```

Output can be seen here:

```
$ tree output
output
|-- bamToBed
|   `-- test2.bed
|-- bamToBedgraph
|   `-- test2.bedgraph
|-- bedToBigBed
|   `-- test2.bigbed
|-- bedgraphToBigWig
|   `-- test2.bw
|-- bowtie2
|   |-- test2.bam
|   `-- test2.bam.bai
|-- fastqc
|   |-- test2_fastqc.html
|   `-- test2_fastqc.zip
|-- printfoo
|-- python-task
`-- txt2fastq
```

# Why `SkidSteer`?
Because the best pipeline is one you built yourself.
