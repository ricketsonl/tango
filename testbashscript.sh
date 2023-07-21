#!/bin/bash

python runDummyXGC.py &
python shestakov_files.py

rm converged.done
