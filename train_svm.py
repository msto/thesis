#!/usr/bin/env python


import numpy as np
import pandas as pd
from sklearn import tree
from sklearn import cross_validation as cv
from sklearn import svm
from sklearn import preprocessing as pp
from sklearn import grid_search as gs


def svm_train(bedpe, svtype):
    col_names = ['chrA', 'startA', 'endA', 'chrB', 'startB', 'endB', 'name',
                 'score', 'strandA', 'strandB', 'num_reads', 'mapqA', 'mapqB',
                 'uniqA', 'uniqB', 'global', 'localA', 'localB', 'poly_count',
                 'samples', 'call', 'branch']

    sv = pd.read_table(bedpe, names=col_names, header=None)
    sv['spanA'] = np.absolute(sv['endA'] - sv['startA'])
    sv['spanB'] = np.absolute(sv['endB'] - sv['startB'])
    sv['local_diff'] = np.absolute(sv['localA'] - sv['localB'])
    sv['mapq_diff'] = np.absolute(sv['mapqA'] - sv['mapqB'])

    # Compute size of varation by taking difference of innermost read positions
    if svtype == 'del' or svtype == 'dup':
        # End is innermost point on forward strand; start is innermost
        # on reverse
        A_mask = sv['strandA'].str.contains('\+')
        B_mask = sv['strandB'].str.contains('\+')
        innerA = A_mask * sv['endA'] + (1 - A_mask) * sv['startA']
        innerB = B_mask * sv['endB'] + (1 - B_mask) * sv['startB']

        # if svtype == 'inv':
        #     A_mask = sv['strandA'].str.contains('\.')
        #     B_mask = sv['strandB'].str.contains('\.')
        #     minA = np.amin(sv[['startA', 'endA']], axis=1)
        #     maxB = np.amax(sv[['startB', 'endB']], axis=1)

        #     innerA = A_mask * minA + (1 - A_mask) * innerA
        #     innerB = B_mask * maxB + (1 - B_mask) * innerB

        sv['size'] = np.absolute(innerB - innerA)

        features = ['num_reads', 'mapqA', 'mapqB', 'uniqA', 'uniqB', 'global',
                    'localA', 'localB', 'poly_count', 'local_diff',
                    'mapq_diff', 'size']
    else:
        features = ['num_reads', 'mapqA', 'mapqB', 'uniqA', 'uniqB', 'global',
                    'localA', 'localB', 'poly_count', 'local_diff',
                    'mapq_diff']

    # Feature vectors
    X = sv[features]
    X = X.values
    X_scaled = pp.scale(X)

    # Labels
    Y = sv['call']
    Y = Y.replace(to_replace='Valid', value=1)
    Y = Y.replace(to_replace='Invalid', value=0)
    Y = Y.values

