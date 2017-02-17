#!/usr/bin/env python
# python 2.7


#

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

def load_tabular(path, sep='\t', header=0, na_values=['.'], *args, **kwargs):
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
        x = load_tabular(self.path)
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
        if isinstance(x, experimentFile): self.files[x.path] = x
        if not isinstance(x, experimentFile): self.files[x] = experimentFile('unspecified', x)
        # self.files.append(x)

class experiment:
    # contains samples
    def __init__(self, name, samplesheet):
        self.name = name
        self.samples = {}
        self.samplesheet = sampleSheet(samplesheet).load()
        self.files = []
    def addSample(self, sample):
        if isinstance(sample, experimentSample): self.samples[sample.name] = sample
        if not isinstance(sample, experimentSample): self.samples[sample] = experimentSample(sample)
    def listSamples(self):
        for key in self.samples.keys():
            print key



x = queryTable("16-05/IonXpress_003/IonXpress_003_query.tsv")

s = experimentSample("16")
s.addFile(x)

a = experiment("MyWonderfulExperiment", "16-05/sample_barcode_IDs.tsv")
a.addSample(s)
a.addSample("foo")

a.samples['16'].files['16-05/IonXpress_003/IonXpress_003_query.tsv'].load()

my_debugger(globals().copy())
