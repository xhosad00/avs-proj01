#!/bin/bash

rm -f line.npz img.png

# run  calculation with line calculator
./mandelbrot -c line -s 1024 -i 100 line.npz

# create image
python3 ../scripts/visualise.py line.npz --save img.png

# Check if the process was successful
if [[ -f "line.npz" && -f "img.png" ]]; then
  echo "Script completed successfully. Output saved to line.npz and img.png."
else
  echo "Error: Script did not complete successfully or files were not created."
fi
