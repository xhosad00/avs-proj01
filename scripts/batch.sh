#!/bin/bash

rm -f batch.npz img.png

# run  calculation with batch calculator
./mandelbrot -c batch -s 1024 -i 100 batch.npz

# create image
python3 ../scripts/visualise.py batch.npz --save img.png

# Check if the process was successful
if [[ -f "batch.npz" && -f "img.png" ]]; then
  echo "Script completed successfully. Output saved to batch.npz and img.png."
else
  echo "Error: Script did not complete successfully or files were not created."
fi
