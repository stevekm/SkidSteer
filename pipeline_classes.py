#!/usr/bin/env python

from base_classes import *
from pipeline_functions import *

class printFoo(ShellTask):
    '''
    A class that invokes 'echo' from the shell
    '''
    def __init__(self, *args, **kwargs):
        ShellTask.__init__(self, "printfoo", *args, **kwargs)
        self.sys_command = '''
echo foo
        '''

class convertFastq2Txt(ShellTask):
    '''
    Convert a .fastq file to a .txt file
    '''
    def __init__(self, *args, **kwargs):
        ShellTask.__init__(self, "txt2fastq", *args, **kwargs)
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
    def __init__(self, *args, **kwargs):
        ShellTask.__init__(self, 'fastqc', input_suffix = ".fastq", *args, **kwargs)
        self.environment = 'module load fastqc/0.11.4'
        self.sys_command ='fastqc'
        self.threads = 4
        self.params = '--nogroup'
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

class alignBowtie2(ShellTask):
    '''
    A task that runs bowtie2 and samtools to align a fastq then convert to .bam format
    '''
    def __init__(self, *args, **kwargs):
        ShellTask.__init__(self, 'bowtie2', input_suffix = ".fastq", output_suffix = ".bam", *args, **kwargs)
        self.sys_command ='bowtie2'
        self.threads = 1
        self.params = '--local -x "~/ref/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome"'
        self.mem = "2G"
    def build_command(self):
        '''
        Method for building the bowtie2 command
        ex: bowtie2 --threads "$THREADS" --local -x "$tmpGenome" -q -U "$INPUTFILE" | samtools view -@ "$THREADS" -Sb1 - | samtools sort -m 10G -@ "$THREADS" - "$OUTFILE"
        '''
        # set the base command and params
        command = '{}'.format(self.sys_command)
        command = '{} {}'.format(command, self.params)
        # list to hold command for each input file
        command_list = []
        for input_fastq in self.input_files:
            output_bam = transform_file(input_file = input_fastq, input_suffix = self.input_suffix, output_dir = self.output_dir, output_suffix = self.output_suffix)
            file_command = '{} -q -U {} --threads {} | samtools view -@ "{}" -Sb1 - | samtools sort -m {} -@ "{}" - -o  "{}"\nsamtools index "{}"'.format(command, input_fastq, self.threads, self.threads, self.mem, self.threads, output_bam, output_bam)
            command_list.append(file_command)
        # make a single command string for all the input files
        command = '\n'.join(command_list)
        return(command)


class doPythonCode(PythonTask):
    '''
    An example of a class that executes some custom Python code
    '''
    def __init__(self, *args, **kwargs):
        PythonTask.__init__(self, 'python-task', *args, **kwargs)
    def build_command(self):
        print("This is some Python code")


class convertFile(PythonTask):
    '''
    Write lines from input file to output file in a reversed order
    '''
    def __init__(self, *args, **kwargs):
        PythonTask.__init__(self, 'convert-file', *args, **kwargs)
    def build_command(self):
        for file in self.input_file:
            with open(file) as input_file, open('newtest', 'w') as new:
                for line in old:
                    if line.rsplit('|', 1)[-1].strip() == 'number3':
                        new.write('this is replacement|number7\n')
                    else:
                        new.write(line)

class bamToBed(ShellTask):
    '''
    Convert a .bam file to .bed format with bedtools
    '''
    def __init__(self, *args, **kwargs):
        ShellTask.__init__(self, 'bamToBed', input_suffix = ".bam", output_suffix = ".bed", *args, **kwargs)
        self.sys_command ='bamToBed'
        self.pre_command = 'export LC_COLLATE=C'
        import os
        self.tmp_file = os.path.join(self.output_dir, "tmp")
    def build_command(self):
        '''
        Method for building the bamToBed command
        ex: bamToBed -i test.bam > test.bed
        '''
        # set the base command and params
        command = '{}'.format(self.sys_command)
        # list to hold command for each input file
        command_list = []
        for input_bam in self.input_files:
            output_bed = transform_file(input_file = input_bam, input_suffix = self.input_suffix, output_dir = self.output_dir, output_suffix = self.output_suffix)
            file_command = '{} -i {} > {}'.format(command, input_bam, output_bed)
            file_command = '{}\nsort -k1,1 -k2,2n {} > {} && /bin/mv {} {}'.format(file_command, output_bed, self.tmp_file, self.tmp_file, output_bed)
            command_list.append(file_command)
        # make a single command string for all the input files
        command = '\n'.join(command_list)
        # add the pre_processing step to the command
        command = '{}\n{}'.format(self.pre_command, command)
        return(command)

class bedToBigBed(ShellTask):
    '''
    Convert a BED file to BigBED format
    '''
    def __init__(self, *args, **kwargs):
        ShellTask.__init__(self, 'bedToBigBed', input_suffix = ".bed", output_suffix = ".bigbed", *args, **kwargs)
        self.sys_command ='bedToBigBed'
        self.chrom_sizes = "data/hg19.chrom.sizes_full.txt"
    def build_command(self):
        '''
        Method for building the bedToBigBed command
        ex: bedToBigBed test.bed "$chrom_sizes" test.bigBed
        '''
        # set the base command and params
        command = '{}'.format(self.sys_command)
        # list to hold command for each input file
        command_list = []
        for input_bed in self.input_files:
            output_bed = transform_file(input_file = input_bed, input_suffix = self.input_suffix, output_dir = self.output_dir, output_suffix = self.output_suffix)
            file_command = '{} "{}" "{}" "{}"'.format(command, input_bed, self.chrom_sizes, output_bed)
            command_list.append(file_command)
        # make a single command string for all the input files
        command = '\n'.join(command_list)
        return(command)


class bamToBedgraph(ShellTask):
    '''
    Convert a BAM file to Bedgraph format
    '''
    def __init__(self, *args, **kwargs):
        ShellTask.__init__(self, 'bamToBedgraph', input_suffix = ".bam", output_suffix = ".bedgraph", *args, **kwargs)
        self.sys_command ='bedtools genomecov'
        self.chrom_sizes = "data/hg19.chrom.sizes_full.txt"
        self.pre_command = 'export LC_COLLATE=C'
        import os
        self.tmp_file = os.path.join(self.output_dir, "tmp")
    def build_command(self):
        '''
        Method for building the bam to bedgraph command
        ex: bedtools genomecov -ibam test.bam -g "$chrom_sizes" -bg > test.bedgraph
        '''
        # set the base command and params
        command = '{}'.format(self.sys_command)
        # list to hold command for each input file
        command_list = []
        for input_bam in self.input_files:
            output_bedgraph = transform_file(input_file = input_bam, input_suffix = self.input_suffix, output_dir = self.output_dir, output_suffix = self.output_suffix)
            file_command = '{} -ibam "{}" -g "{}" -bg > "{}"'.format(command, input_bam, self.chrom_sizes, output_bedgraph)
            file_command = '{}\nsort -k1,1 -k2,2n {} > {} && /bin/mv {} {}'.format(file_command, output_bedgraph, self.tmp_file, self.tmp_file, output_bedgraph)
            command_list.append(file_command)
        # make a single command string for all the input files
        command = '\n'.join(command_list)
        # add the pre_processing step to the command
        command = '{}\n{}'.format(self.pre_command, command)
        return(command)


class bedgraphToBigWig(ShellTask):
    '''
    Convert a Bedgraph file to BigWig format
    '''
    def __init__(self, *args, **kwargs):
        ShellTask.__init__(self, 'bedgraphToBigWig', input_suffix = ".bedgraph", output_suffix = ".bw", *args, **kwargs)
        self.sys_command ='bedGraphToBigWig'
        self.chrom_sizes = "data/hg19.chrom.sizes_full.txt"

    def build_command(self):
        '''
        Method for building the bam to bedgraph command
        ex: bedGraphToBigWig test.bedgraph "$chrom_sizes" test.bigwig
        '''
        # set the base command and params
        command = '{}'.format(self.sys_command)
        # list to hold command for each input file
        command_list = []
        for input_bg in self.input_files:
            output_bw = transform_file(input_file = input_bg, input_suffix = self.input_suffix, output_dir = self.output_dir, output_suffix = self.output_suffix)
            file_command = '{} "{}" "{}" "{}"'.format(command, input_bg, self.chrom_sizes, output_bw)
            command_list.append(file_command)
        # make a single command string for all the input files
        command = '\n'.join(command_list)
        return(command)
