import os
import sys
from pathlib import Path


directory = os.path.curdir
root = Path('/Users/rionbr/Downloads/All_hits_processed')

folders = [path for path in root.iterdir() if path.is_dir()]

for folder in folders:
    subfolder = folder / 'processed_8bit_weighted'

    # MOVE
    for file in subfolder.iterdir():
        file.rename(folder / file.name)

    # REMOVE subfolder
    subfolder.rmdir()

    # RENAME folder
    if '8bit_weighted' in folder.name:
        folder.rename(folder.parent / folder.name.split('_')[0])