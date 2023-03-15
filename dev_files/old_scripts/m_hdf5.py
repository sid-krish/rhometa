#!/usr/bin/env python
import glob
import logging
import os
import warnings
from collections import namedtuple
from collections.abc import Iterable
from functools import partial

import pandas as pd
from ldpop import rhos_from_string

warnings.filterwarnings('ignore', category=pd.io.pytables.PerformanceWarning)

logger = logging.getLogger(__name__)


class NotFoundException(Exception):
    def __init__(self, message):
        super().__init__(message)


class DuplicateRunsException(Exception):
    def __init__(self):
        super().__init__('Catalog contains duplicate run entries')


class ForeignKeyException(Exception):
    def __init__(self, message):
        super().__init__(message)


class TableStore(object):
    _CATALOG_KEY = 'catalog'
    _CATALOG_COLUMNS = ['theta', 'rhos', 'sample_size', 'pop_sizes', 'times', 'is_approx']
    SPEC_T = namedtuple('Row', ['theta', 'rhos', 'sample_size', 'pop_sizes', 'times', 'is_approx'])

    @staticmethod
    def _open(file_name, mode):
        return pd.HDFStore(file_name, mode=mode, complevel=9)

    @staticmethod
    def _run(theta, rhos, sample_size=1, pop_sizes=[1.0, ], times=[], is_approx=True):
        assert isinstance(rhos, Iterable), 'rhos must be iterable (eg. a collection)'
        assert isinstance(pop_sizes, Iterable), 'pop_sizes must be iterable (eg. a collection)'
        assert isinstance(times, Iterable), 'times must be iterable (eg. a collection)'
        return TableStore.SPEC_T(theta, tuple(rhos), sample_size, tuple(pop_sizes), tuple(times), is_approx)

    @staticmethod
    def _make_fk(_run_idx, _depth):
        try:
            _run_idx = int(_run_idx)
            _depth = int(_depth)
        except ValueError:
            raise ForeignKeyException('Either run index or depth was not castable to "int"')
        return 'run_{}/depth_{}'.format(_run_idx, _depth)

    def __init__(self, file_name, mode='w+'):
        if not file_name.endswith('.h5'):
            file_name = '{}.h5'.format(file_name)
        self.file_name = file_name
        self.mode = mode

    def __enter__(self):
        self.h5store = self._open(self.file_name, self.mode)
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()

    def close(self):
        if self.h5store is not None:
            self.h5store.close()

    def init_store(self):
        if TableStore._CATALOG_KEY not in self.h5store:
            logger.info('Initialising new h5 store as file: {}'.format(self.file_name))
            _empty_df = pd.DataFrame(columns=TableStore._CATALOG_COLUMNS)
            self.h5store[TableStore._CATALOG_KEY] = _empty_df
        else:
            logger.warning('TableStore instance already initialized')

    def _find_run(self, run, must_exist=False):
        df = self.get_run_catalog().query(
            'theta==@run.theta and '
            'rhos==@run.rhos and '
            'sample_size==@run.sample_size and '
            'pop_sizes==@run.pop_sizes and '
            'times==@run.times and '
            'is_approx==@run.is_approx')
        if must_exist and len(df) != 1:
            raise NotFoundException('Run {} did not exist in catalog'.format(run))
        if len(df) >= 2:
            raise DuplicateRunsException()
        if len(df) == 0:
            return None
        return df.index[0], df.iloc[0]

    def _run_exists(self, run):
        try:
            self._find_run(run, True)
        except NotFoundException:
            return False

    def get_run_catalog(self):
        return self.h5store[TableStore._CATALOG_KEY]

    def add_tableset(self, downsample_path, theta, rhos, sample_size=1,
                     pop_sizes=[1.0, ], times=[], is_approx=True, file_type='pickle'):

        _catalog = self.get_run_catalog()
        _new_run = TableStore._run(theta, rhos, sample_size, pop_sizes, times, is_approx)
        if self._run_exists(_new_run):
            raise IOError('Catalog already contains a run with the specification: {}'.format(_new_run))

        # replace catalog with update
        logger.info('Adding new run specification to catalog: {}'.format(_new_run))
        _catalog = pd.concat([_catalog, pd.DataFrame([_new_run])], ignore_index=True, sort=False)
        self.h5store[TableStore._CATALOG_KEY] = _catalog

        # get the index of the run
        _run_idx, _run = self._find_run(_new_run, True)

        if file_type == 'pickle':
            file_ext = '.pkl'
            _reader = pd.read_pickle
        elif file_type == 'csv':
            file_ext = '.csv'
            _reader = partial(pd.read_csv, index_col=0)
        else:
            raise IOError('unknown table type: {}'.format(file_type))

        input_tables = [(int(os.path.splitext(fn)[0].split('_')[-1]), fn) for fn in
                        glob.glob('{}/*{}'.format(downsample_path, file_ext))]
        input_tables = sorted(input_tables, key=lambda x: x[0])
        assert len(input_tables) > 0, 'no input tables were found on path: {}'.format(downsample_path)
        depth_list = []
        for _depth, _fn in input_tables:
            logger.info('Reading depth: {}, from file: {}'.format(_depth, _fn))
            _df = _reader(_fn)
            _df.index.set_names('00 01 10 11', inplace=True)
            depth_list.append(_depth)
            rho_range = tuple(_df.columns.astype(float).to_series().agg([min, max]))
            assert rho_range == (min(rhos), max(rhos)), 'rho range in table disagrees with arguments: {}'.format(_fn)
            self.add_table(_run_idx, _df, _depth)

    def add_table(self, _run_idx, lk_table, depth):
        assert isinstance(lk_table, pd.DataFrame)
        lk_table['foreign_key'] = _run_idx
        _key = TableStore._make_fk(_run_idx, depth)
        assert _key not in self.h5store, 'Table for depth {} already exists'.format(depth)
        self.h5store[_key] = lk_table

    def fetch_table(self, run, depth):
        _run_idx, _run = self._find_run(run, True)
        _key = TableStore._make_fk(_run_idx, depth)
        if _key not in self.h5store:
            raise NotFoundException('Lk table not found for run {} and depth {}'.format(run, depth))
        return self.h5store[_key]


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser('HDF5 Lookup Table Store')
    parser.add_argument('--verbose', default=False, action='store_true', help='Verbose logging')
    parser.add_argument('--table-format', default='csv', help='Downsampled lookup table format',
                        choices=['csv', 'pickle'])
    parser.add_argument('--init', default=False, action='store_true', help='Initialise a new h5 store')
    parser.add_argument('-t', '--theta', type=float, help='Theta used in table generation')
    parser.add_argument('-r', '--rho-range', help='Range of rho values used in table generation')
    parser.add_argument('-n', '--sample-size', default=1, type=int, help='Sample size used in table generation')
    parser.add_argument('lk_table_path', help='Path to the set of downsampled lookup tables')
    parser.add_argument('h5_store_file', help='Output file name for H5 store')
    args = parser.parse_args()

    _mode = 'w' if args.init else 'a'

    logging.captureWarnings(True)
    logger = logging.getLogger('main')

    # root log listens to everything
    root = logging.getLogger('')
    root.setLevel(logging.DEBUG)
    logging.getLogger('numba').setLevel(logging.INFO)

    # log message format
    formatter = logging.Formatter(fmt='%(levelname)-8s | %(asctime)s | %(name)7s | %(message)s')

    # Runtime console listens to INFO by default
    ch = logging.StreamHandler()
    if args.verbose:
        ch.setLevel(logging.DEBUG)
    else:
        ch.setLevel(logging.INFO)
    ch.setFormatter(formatter)
    root.addHandler(ch)

    with TableStore(args.h5_store_file, mode=_mode) as table_store:
        if args.init:
            table_store.init_store()
        rho_range = rhos_from_string(args.rho_range)
        table_store.add_tableset(args.lk_table_path, args.theta, rho_range, args.sample_size,
                                 [1.0, ], [], True, file_type=args.table_format)
