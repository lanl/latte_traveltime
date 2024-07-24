#
# © 2024. Triad National Security, LLC. All rights reserved.
#
# This program was produced under U.S. Government contract 89233218CNA000001 
# for Los Alamos National Laboratory (LANL), which is operated by 
# Triad National Security, LLC for the U.S. Department of Energy/National Nuclear 
# Security Administration. All rights in the program are reserved by 
# Triad National Security, LLC, and the U.S. Department of Energy/National 
# Nuclear Security Administration. The Government is granted for itself and 
# others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide 
# license in this material to reproduce, prepare. derivative works, 
# distribute copies to the public, perform publicly and display publicly, 
# and to permit others to do so.
#
# Author:
#    Kai Gao, kaigao@lanl.gov
#


import sys
import os
from os.path import expanduser

# Get file size in bytes
def get_filesize(file):
    return os.path.getsize(file)

# Check file existence
def check_filexist(file):
    return os.path.exists(file)

# Make directory
def make_directory(dir):
	os.makedirs(dir, exist_ok=True)
	
# Get home directory
def get_homedir():
	return expanduser("~")