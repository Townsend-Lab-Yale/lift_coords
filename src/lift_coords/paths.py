import os
import pathlib

import appdirs


_app_dirs = appdirs.AppDirs('lift_coords')
DATA_ROOT = _app_dirs.user_data_dir
CHAIN_DIR = pathlib.Path(DATA_ROOT).joinpath('data')
WORK_DIR = pathlib.Path(DATA_ROOT).joinpath('temp')

os.makedirs(CHAIN_DIR, exist_ok=True)
os.makedirs(WORK_DIR, exist_ok=True)
