#!/bin/bash

set -e

bindir=../../bin

# Create the model, the geometry and setup.json
python3 create_model_and_geometry.py

# Gradients with adjoint_receiver = nearest and linear, and the
# finite-difference reference, with four LATTE runs at a time
python3 compare_gradients.py --bindir $bindir --jobs 4

# Print the errors and plot the figures
python3 plot.py
