import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn; seaborn.set_style('whitegrid')
import itertools
import numpy

import os
import pkgutil


if not pkgutil.find_loader('tensorflow'):
    os.environ['KERAS_BACKEND'] = 'theano'

from avocado import Avocado

celltypes = ['sample' + str(i) for i in np.arange(1, 164)]
assays = ['H3K27ac', 'H3K27me3', 'H3K9me3','H3K36me3','H3K4me1','H3K4me3']

test_set_sample = ['sample' + str(i) for i in np.arange(128, 164)]
test_set_assays = ['H3K36me3','H3K4me1','H3K4me3']

data = {}
for celltype, assay in itertools.product(celltypes, assays):
    if (celltype in test_set_sample) and (assay in test_set_assays):
        continue
        
    filename = '/data/npz/{}_{}.pilot.arcsinh.npz'.format(celltype, assay)
    data[(celltype, assay)] = numpy.load(filename)['arr_0']

model = Avocado(celltypes, 
                assays, 
                n_layers=1, 
                n_nodes=128, 
                n_assay_factors=24, 
                n_celltype_factors=32,
                n_25bp_factors=10, 
                n_250bp_factors=20, 
                n_5kbp_factors=30, 
                n_genomic_positions=273119,
                batch_size=10000)

model.fit(data, n_epochs=10, epoch_size=28)
model.save("/model/test-model")

for celltype, assay in itertools.product(test_set_sample, test_set_assays):
    y_hat  = model.predict(celltype,assay)
    np.savez('/result/{}_{}.predict.npz'.format(celltype,assay), y_hat)
