#!/bin/bash

~/anaconda3/bin/python runDummyXGC.py &
~/anaconda3/bin/python shestakov_files.py

rm converged.done
