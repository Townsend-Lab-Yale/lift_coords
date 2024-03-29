#!/usr/bin/env python
import os
import shlex
import logging
import subprocess
from datetime import datetime

import pandas as pd

from . import paths
from .admin import copy_initial_data


__author__ = 'Stephen Gaffney'
_logger = logging.getLogger(__name__)
VALID_BUILDS = frozenset(['grch37', 'grch38', 'hg19', 'hg38'])
_CHAIN_HG19_37 = 'hg19_to_GRCh37.chain'
_CHAIN_HG19_HG38 = 'hg19_to_hg38.chain.gz'
_CHAIN_HG38_HG19 = 'hg38_to_hg19.chain'
_CHAIN_HG38_38 = 'hg38_to_GRCh38.chain'
_CHAIN_37_HG38 = 'GRCh37_to_hg38.chain'
_CHAIN_37_HG19 = 'GRCh37_to_hg19.chain'
_CHAIN_37_38 = 'GRCh37_to_GRCh38.chain'
_CHAIN_38_37 = 'GRCh38_to_GRCh37.chain'
_CHAIN_38_HG38 = 'GRCh38_to_hg38.chain'

_CHAIN_LISTS = {
    ('grch37', 'grch38'): [_CHAIN_37_38],
    ('grch37', 'hg19'): [_CHAIN_37_HG19],
    ('grch37', 'hg38'): [_CHAIN_37_HG38],
    ('grch38', 'grch37'): [_CHAIN_38_37],
    ('grch38', 'hg19'): [_CHAIN_38_37, _CHAIN_37_HG19],
    ('grch38', 'hg38'): [_CHAIN_38_HG38],
    ('hg19', 'grch37'): [_CHAIN_HG19_37],
    ('hg19', 'grch38'): [_CHAIN_HG19_HG38, _CHAIN_HG38_38],
    ('hg19', 'hg38'): [_CHAIN_HG19_HG38],
    ('hg38', 'grch37'): [_CHAIN_HG38_HG19, _CHAIN_HG19_37],
    ('hg38', 'grch38'): [_CHAIN_HG38_38],
    ('hg38', 'hg19'): [_CHAIN_HG38_HG19],
}


def lift_over(df, build_in: str, build_out: str, keep_orig=False,
              keep_intermediate=False):
    """Lift table of site info from one build to another.

    Args:
        df (pd.DataFrame): table that includes chrom + position coordinates
        build_in: Name of input build. See VALID_BUILDS for options.
        build_out: Name of output build.
        keep_orig (bool): whether to keep initial coordinate columns as {col}_orig.
        keep_intermediate (bool): whether to keep intermediate liftOver bed files.
    Returns:
        new_table (pd.DataFrame), inds_unlifted
    """
    if build_in not in VALID_BUILDS or build_out not in VALID_BUILDS:
        raise BuildArgumentException(f"Build must be one of {VALID_BUILDS}.")
    new_table, unlifted = _lift_with_chains(
        df, keep_orig=keep_orig,
        chain_list=_CHAIN_LISTS[(build_in, build_out)], new_build_name='GRCh37',
        keep_intermediate=keep_intermediate)
    return new_table, unlifted


def _prettify_build_str(build: str) -> str:
    """Get build string suited to MAF Build column.

    >>> _prettify_build_str('grch37')
    'GRCh37'

    >>> _prettify_build_str('hg19')
    'hg19'
    """
    if build.startswith('grch'):
        return build.replace('grch', 'GRCh')
    return build


def _lift_with_chains(df, keep_orig=False, chain_list=None,
                      new_build_name=None, keep_intermediate=False):
    """Lift in two steps, using two chain files (chain1, and chain2).

    Args:
        df (pd.DataFrame): data table that includes chrom + position coordinates
        keep_orig (bool): whether to keep initial coordinate columns as {col}_orig
        chain_list (list): chain file basenames, e.g. [_CHAIN_HG38_HG19,]
        new_build_name (str): OPTIONAL. name to go in new Build column
            (if input table has a build column)
    """
    if chain_list is None:
        raise TypeError("At least one chain file is required.")
    build_col = _match_column_str('build', df)
    chr_col = _match_column_str('chrom', df)
    start_col = _match_column_str('start', df)
    end_col = _match_column_str('end', df)
    end_col = end_col if end_col else start_col
    use_cols = [chr_col, start_col, end_col]

    df_orig = df
    df = df[use_cols].copy()

    # note: index must be unique
    df['ind'] = df.index.values

    bed_path = _make_bed(chr_col, start_col, end_col, index_col='ind',
                         data=df, out_dir=paths.WORK_DIR)

    current_bed = bed_path
    out_paths, unlift_paths = [], []
    for ind, chain in enumerate(chain_list):
        chain_path = os.path.join(paths.CHAIN_DIR, chain)
        out_path, unlift_path = _lift_bed(current_bed, chain_path=chain_path,
                                          out_dir=paths.WORK_DIR)
        out_paths.append(out_path)
        unlift_paths.append(unlift_path)
        current_bed = out_path

    n = pd.read_csv(current_bed, header=None, sep='\t',
                    names=[chr_col, start_col, end_col, 'ind'])
    n.Start_Position += 1
    n.set_index('ind', inplace=True)

    if not keep_intermediate:
        _logger.info('Cleaning up intermediate files.')
        for path in out_paths + unlift_paths + [bed_path]:
            os.remove(path)

    df2 = df_orig.join(n, how='left', lsuffix='_orig')

    unlifted = df_orig[df2.Start_Position.isnull()]
    _logger.info(f'{len(unlifted)} rows failed liftover.')

    df2.dropna(axis=0, how='any', subset=[start_col, end_col], inplace=True)
    df2[start_col] = df2[start_col].astype(int)
    df2[end_col] = df2[end_col].astype(int)

    if build_col:
        if keep_orig:
            df2.rename(columns={build_col: build_col + '_orig'}, inplace=True)
            if new_build_name:
                df2[build_col] = new_build_name
        else:
            if new_build_name:
                df2[build_col] = new_build_name
            df2.drop([i + '_orig' for i in [chr_col, start_col, end_col]], axis=1,
                     inplace=True)

    return df2, unlifted


def _lift_bed(bed_in_path, chain_path=None, out_dir=None):
    """Lifts bed file using provided chain."""

    if not chain_path:
        raise Exception("chain_path required.")
    if not out_dir:
        out_dir = os.getcwd()

    # n_pos = len(chrom)
    # if n_pos != len(start_pos) != len(end_pos):
    #     raise Exception("Parameters must have same length.")
    now_str = datetime.isoformat(datetime.now())
    out_path = os.path.join(out_dir, f'lift_{now_str}_out')
    unlifted_path = os.path.join(out_dir, f'lift_{now_str}_unlifted')

    # liftOver oldFile map.chain newFile unMapped
    cmd = ['liftOver', bed_in_path, chain_path, out_path, unlifted_path]
    _logger.info(f"Running {' '.join(shlex.quote(arg) for arg in cmd)}")

    with open(os.devnull, "r") as fnullin:
        with subprocess.Popen(cmd, stdin=fnullin, stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE,
                              universal_newlines=True) as p:
            for line in p.stdout:
                msg = line.rstrip('\n')
                _logger.info(f"LIFTOVER: {msg}")
            for line in p.stderr:
                msg = line.rstrip('\n')
                _logger.error(f"LIFTOVER: {msg}")
    return out_path, unlifted_path


def _make_bed(chrom='Chromosome', start_pos='Start_Position',
              end_pos='End_Position', index_col='ind', data=None,
              out_dir=None):
    """Create bed file from dataframe."""
    if not out_dir:
        out_dir = os.getcwd()

    now_str = datetime.isoformat(datetime.now())
    bedfilepath = os.path.join(out_dir, 'lift_' + now_str + '_in')

    bed_df = data[[chrom, start_pos, end_pos, index_col]].copy()
    bed_df[start_pos] = bed_df[start_pos] - 1
    bed_df.to_csv(bedfilepath, sep='\t', header=False, index=False)

    return bedfilepath


def _match_column_str(val=None, df=None):
    """Get first column which, in lowercase, contains input text."""

    matches = [i for i in df.columns if val in i.lower()]
    if matches:
        return matches[0]
    else:
        return None


class BuildArgumentException(Exception):
    pass

copy_initial_data()
