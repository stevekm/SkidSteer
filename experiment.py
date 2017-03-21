#!/usr/bin/env python
# python 2.7


import sys
import os

def my_debugger(vars):
    # starts interactive Python terminal at location in script
    # call with pl.my_debugger(globals().copy()) anywhere in your script
    # or call my_debugger(locals().copy()) from anywhere within this package
    import readline # optional, will allow Up/Down/History in the console
    import code
    # vars = globals().copy() # in python "global" variables are actually module-level
    vars.update(locals())
    shell = code.InteractiveConsole(vars)
    shell.interact()

def load_tabular(path, *args, **kwargs):
    # load a tabular file into a dataframe with pandas
    import pandas as pd
    df = pd.read_table(path, *args, **kwargs)
    return df


class experimentFile:
    # a file associated with the experiment
    def __init__(self, filetype, path):
        self.filetype = filetype
        self.path = path

class queryTable(experimentFile):
    # a 'query table' made from a VCF file
    def __init__(self, path):
        experimentFile.__init__(self, filetype = "query_table", path = path)
        self.path = path
    def load(self):
        x = load_tabular(self.path, sep='\t', header=0, na_values=['.'])
        return x

class sampleSheet(experimentFile):
    # a samplesheet for an experiment
    def __init__(self, path):
        experimentFile.__init__(self, filetype = "sample_sheet", path = path)
        self.path = path
    def load(self):
        x = load_tabular(self.path, header=1)
        return x

class experimentSample:
    # a sample in the experiment
    def __init__(self, name):
        self.name = name
        self.files = {}
    def addFile(self, x):
        if isinstance(x, experimentFile): self.files[x.path] = os.path.abspath(x.path)
        if not isinstance(x, experimentFile): self.files[x] = experimentFile('unspecified', os.path.abspath(x))
        # self.files.append(x)

class experiment:
    # contains samples
    def __init__(self, name, samplesheet):
        self.name = name
        self.samples = {}
        self.samplesheet = sampleSheet(samplesheet).load()
        self.files = []
        self.dir = os.path.dirname(samplesheet)
    def addSample(self, sample):
        if isinstance(sample, experimentSample): self.samples[sample.name] = sample
        if not isinstance(sample, experimentSample): self.samples[sample] = experimentSample(sample)
    def listSamples(self):
        for key in self.samples.keys():
            print key

def parse_sample_dir(mypath):
    # get the samples from a dir; one sample per subdir in the parent dir
    # my_debugger(locals().copy())
    samples_list = []
    dir_list = [d for d in os.listdir(mypath) if os.path.isdir(os.path.join(mypath, d))]
    # experimentSample(f)
    for d in dir_list:
        mydirpath = os.path.join(mypath,d)
        file_list = [f for f in os.listdir(mydirpath) if os.path.isfile(os.path.join(mydirpath, f))]
        mysample = experimentSample(d)
        for f in file_list: mysample.addFile(f)
        samples_list.append(mysample)
    return samples_list



x = queryTable("16-05/IonXpress_003/IonXpress_003_query.tsv")

s = experimentSample("16")

s.addFile(x)

a = experiment(name = "MyWonderfulExperiment", samplesheet = "16-05/sample_barcode_IDs.tsv")
a.addSample(s)
a.addSample("foo")

my_debugger(globals().copy())
a.samples['16'].files['16-05/IonXpress_003/IonXpress_003_query.tsv'].load()
parse_sample_dir('16-05')
my_debugger(globals().copy())
