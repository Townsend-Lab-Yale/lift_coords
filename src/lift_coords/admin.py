import logging
import os
import shutil
import pathlib
from importlib import resources

from . import paths

_logger = logging.getLogger(__name__)


def copy_initial_data():
    os.makedirs(paths.DATA_ROOT, exist_ok=True)
    added_data = []
    for file_name in [
        'GRCh37_to_GRCh38.chain.gz',
        'GRCh37_to_hg38.chain.gz',
        'GRCh38_to_GRCh37.chain.gz',
        'hg19_to_GRCh37.chain.gz',
        'hg19_to_hg38.chain.gz',
        'hg38_to_GRCh38.chain.gz',
        'hg38_to_hg19.chain.gz',
    ]:
        new_path = pathlib.Path(paths.CHAIN_DIR).joinpath(file_name)
        if not new_path.exists():
            added_data.append(file_name)
            with resources.path('lift_coords.refs.data', file_name) as path:
                shutil.copy(path, new_path)
    if added_data:
        _logger.info(f"Copied reference data to {paths.DATA_ROOT}: {added_data}.")
