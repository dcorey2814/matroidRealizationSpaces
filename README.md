# matroidRealizationSpaces

Code used in the paper <a href="https://arxiv.org/abs/xxx"> Singular matroid realization spaces</a> by Daniel Corey and Dante Luber. This code works with OSCAR version 0.12.1. Note: it may be essential to use this version of OSCAR. To ensure this, do the following. First, run `julia --project=.` in the terminal from the root of this project. Next, open `julia` and run the following:

```
julia> using Pkg
julia> Pkg.instantiate()
```

## Functions and documentation

All new functions are in the `src` directory. The documentation is contained in the `jupyter` notebook `functionDocumentation.ipynb`. **Note** In the file `src/matroid_realization.jl`, some of these functions and the structure `MatroidRealizationSpace`  were developed with Lukas K&uuml;hne, and will be incorporated in a later version of Oscar. 

## Smoothness of realization spaces

The code used to determine smoothness (and in some cases, irreducibility) may be found in the folders `d3n9`, `d3n10`, `d3n11`, `d4n8`, `d4n9`, where, e.g., `d3n9` contains the data for rank $d=3$ matroids on $n=9$ elements. In each folder is a `jupyter` notebook that contains all instructions for the smoothness tests. 

## Singular realization spaces

The `jupyter` notebook `matroids-3-12.ipynb` contains code to verify that the $(3,12)$--matroid in section 4.2 has a singular realization space, as well as code used in Proposition 6.7 to verify that the coordinatewise valuation of the Pl&uuml;cker coordinates of the row spans of the matrices $C_1, C_2$ both equal the corank vector of this matroid. 
