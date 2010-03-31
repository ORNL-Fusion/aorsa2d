#!/bin/bash

mpirun -n 4 xterm -geom 120x40 -e gdb -ex run ~/code/aorsa2d/xaorsa2d
