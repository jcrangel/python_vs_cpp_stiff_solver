# Comparison beetween ODE solvers in c++ and python
This solves the stiff system from :
J. C. Rangel-Reyes, J. C. Chimal-Eguía, and E. Castillo-Montiel, “Dendritic Immunotherapy Improvement for an Optimal Control Murine Model,” Computational and Mathematical Methods in Medicine, vol. 2017, Article ID 5291823, 9 pages, 2017. https://doi.org/10.1155/2017/5291823.


## Requirements 
1. Boost library
2. Eigen library

## c++
Compile with 
```
g++ -w -O2 -std=c++11 DC_therapy.cpp -I/home/rangel/boost_1_69_0/ -I/home/rangel/eigen
```
Plot graph from plot_points.ipynb or 
```
python  plot_points.py
```

## Python 

Run 
```
python DC_therapy.py
```

## Timing
* c++ version took : 2.12742 sec
* python version took: 5 min!
