#!/usr/bin/env python


import argparse
import numpy as np
import pandas as pd


def split_train(bedpe):
    """
    Randomly select 15% of valid and 15% of invalid data to use as test set

    Parameters
    ----------
    bedpe : str
        Filepath to bedpe to sample test set from

    Returns
    -------
    test : pandas.DataFrame
        Randomly selected 15% of input
    train : pandas.DataFrame
        The remaining 85% of input
    """

    col_names = ['chrA', 'startA', 'endA', 'chrB', 'startB', 'endB', 'name',
                 'score', 'strandA', 'strandB', 'size', 'mapqA', 'mapqB',
                 'uniqA', 'uniqB', 'global', 'localA', 'localB', 'poly_count',
                 'samples', 'call', 'branch']

    sv = pd.read_table(bedpe, names=col_names, header=None)

    valids = sv[sv['call'].str.match('Valid')]
    invalids = sv[sv['call'].str.match('Invalid')]

    # Randomize
    valid_rows = np.random.permutation(valids.index.values)
    invalid_rows = np.random.permutation(invalids.index.values)

    # Select 15% of each set to save for testing
    valid_test = valids.ix[valid_rows[:int(0.15 * len(valids.index))]]
    valid_train = valids.ix[valid_rows[int(0.15 * len(valids.index)):]]
    invalid_test = invalids.ix[invalid_rows[:int(0.15 * len(invalids.index))]]
    invalid_train = invalids.ix[invalid_rows[int(0.15 * len(invalids.index)):]]

    test = pd.concat([valid_test, invalid_test])
    train = pd.concat([valid_train, invalid_train])

    return test, train


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('bedpe')
    parser.add_argument('svtype')

    args = parser.parse_args()

    test, train = split_train(args.bedpe)

    test.to_csv('Eco.%s.test.bedpe' % args.svtype,
                sep='\t', header=False, index=False)
    train.to_csv('Eco.%s.train.bedpe' % args.svtype,
                 sep='\t', header=False, index=False)

if __name__ == '__main__':
    main()
