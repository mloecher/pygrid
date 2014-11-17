pygrid
======

WORK IN PROGRESS, DO NOT USE

The goal is to make a python module for gridding, with:
- Different implementations for speed/compatibility balance
    - Numpy vectors
    - Cython/openMP
    - pyCUDA
    - outside libraries linked in
- Flexible kernel and deappodization choices
- Methods to pick the best implementation/kernel for a given trajectory
    - Benchmarking and RMSE measures across the various techniques
