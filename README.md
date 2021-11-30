# Algorithms for computing coding and non-coding dinucleotide forces, in Julia

This repository contains the code to compute dinucleotide forces (for a definition
give a look a [this paper](https://www.pnas.org/content/111/13/5054.short)),
as well as a simple example of its use.
For any question or comment concerning the code please contact Andrea Di Gioacchino,
<andrea.dgioacchino@gmail.com>.

All of the code is written in Julia.

## Depdendencies
[FiniteDiff.jl](https://github.com/JuliaDiff/FiniteDiff.jl) must be installed
to use the script. 
Moreover an installation of [jupyter](https://jupyter.org) is needed to run the `*.ipynb` notebook.

## Repository structure:
- The `NoncodingForces_v2_1.jl` file containts the script to compute forces
in the non-coding case (without caring about codons); 
- The `example.ipynb` notebook contains a short tutorial of the
usage of the scripts.


## Acknowledgements
Previous versions of these algorithms were developed by Simona Cocco,
Benjamin D. Greenbaum, Rémi Monasson, Alexander Solovyov and Petr Šulc.
