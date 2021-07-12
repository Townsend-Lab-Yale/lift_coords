#!/usr/bin/env python
import os
import logging
import subprocess
from datetime import datetime

import pandas as pd

from . import paths


__author__ = 'Stephen Gaffney'

CHAIN_HG38_19 = 'hg38_to_hg19.chain'
CHAIN_37_HG38 = 'GRCh37_to_hg38.chain'
CHAIN_37_19 = 'GRCh37ToHg19.chain'
CHAIN_19_37 = 'hg19_to_GRCh37.chain'
_logger = logging.getLogger(__name__)


def lift_with_chains(df, keep_orig=False, chain_list=None,
                     new_build_name=None):
    """Lift in two steps, using two chain files (chain1, and chain2).

    Args:
        df (pd.DataFrame): data table that includes chrom + position coordinates
        keep_orig (bool): whether to keep initial coordinate columns as {col}_orig
        chain_list (list): chain file basenames, e.g. [CHAIN_HG38_19,]
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

    df2 = df_orig.join(n, how='left', lsuffix='_orig')

    inds_unlifted = df2[df2.Start_Position.isnull()].index.values
    n_unlifted = len(inds_unlifted)
    _logger.info('%s rows failed liftover.', n_unlifted)

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
    # CLEAN UP TEMP FILES
    # os.remove(bed_path)
    # for path in [unlift_path1, unlift_path2, out_path, out2_path]:
    #     os.remove(os.path.join(out_dir, path))

    return df2, inds_unlifted


def lift_38_to_grch37(df, keep_orig=False):
    """Convert hg38 to grch37."""

    df2, inds_unlifted = lift_with_chains(
        df, keep_orig=keep_orig,
        chain_list=[CHAIN_HG38_19, CHAIN_19_37], new_build_name='GRCh37')
    return df2, inds_unlifted


def lift_37_to_hg19(df, keep_orig=False):
    """Convert hg38 to grch37."""

    df2, inds_unlifted = lift_with_chains(
        df, keep_orig=keep_orig,
        chain_list=[CHAIN_37_19], new_build_name='hg19')
    return df2, inds_unlifted


def _lift_bed(bed_in_path, chain_path=None, out_dir=None):
    """Lifts ."""

    if not chain_path:
        raise Exception("chain_path required.")
    if not out_dir:
        out_dir = os.getcwd()

    # n_pos = len(chrom)
    # if n_pos != len(start_pos) != len(end_pos):
    #     raise Exception("Parameters must have same length.")
    now_str = datetime.isoformat(datetime.now())
    out_path = os.path.join(out_dir, 'lift_' + now_str + '_out')
    unlifted_path = os.path.join(out_dir, 'lift_' + now_str + '_unlifted')

    # liftOver oldFile map.chain newFile unMapped
    cmd = ['liftOver', bed_in_path, chain_path, out_path, unlifted_path]
    with open(os.devnull, "r") as fnullin:
        with open(os.devnull, "w") as outfile:
            subprocess.check_call(cmd, stdout=outfile, stderr=subprocess.STDOUT,
                                  stdin=fnullin)
            # subprocess.check_call(cmd,stdout=outfile,stdin=fnullin)
            # proc = subprocess.Popen(cmd,stdout=outfile,stdin=fnullin)
            # (a,b)=proc.communicate()

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