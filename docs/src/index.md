```@meta
CurrentModule = DinucleotideForces
```

# DinucleotideForces

Documentation for [DinucleotideForces](https://github.com/adigioacchino/DinucleotideForces.jl), a Julia package to computed dinucleotide forces in nucleic-acid sequences.

## Installation
This package is not registered. Install with:

```julia
import Pkg
Pkg.add(url="https://github.com/adigioacchino/DinucleotideForces.jl")
```

## Dinucleotide forces: definition
Nucleotide/dinucleotide forces for a sequence are parameters that capture pressures to increase or decrease the usage of a given nucleotide/dinucleotide in the sequence.
So a positive force on a nucleotide corresponds to a pressure to have more occurrences of the given dinucleotide, while a negative pressure corresponds to a lack of occurrences of the dinucleotide in the considered sequence.
Forces are inferred through a [maximum entropy model](https://en.wikipedia.org/wiki/Principle_of_maximum_entropy).
For a more formal definition give a look at [this paper](https://www.pnas.org/content/111/13/5054.short).

## Citing
If this package is used for research purposes, please consider citing [this](https://doi.org/10.1093/molbev/msab036) 
and [this](https://www.pnas.org/content/111/13/5054.short) papers.

## Examples

### Basic usage
We will take sub-sequence from the Influenza H5N1 PB2 segment (strain used: A/Anhui/1/2005) and we will use it as working example to test the function of the scripts.
We start by loading the script and the sequence:
```julia
using DinucleotideForces
example_seq = "ATGGAGAGAATAAAAGAATTAAGGGATCTAATGTCACAGTCCCGCACTCGCGAGATACTAACAAAAACCACTGTGGACCATATGGCCATAATCAAGAAGTACACATCAGGAAGACAAGAGAAGAACCCTGCTCTCAGAATGAAATGGATGATGGCAATGAAATATCCAATCACAGCGGACAAGAGAATAACAGAGATGATTCCTGAAAGGAATGAACAAGGGCAGACGCTCTGGAGCAAGACAAATGATGCCGGATCGGACAGGTTGATGGTGTCTCCCTTAGCTGTAACTTGGTGGAATAGGAATGGGCCGACGACAAGTGCAGTCCATTATCCAAAGGTTTACAAAACATACTTTGAGAAGGCT"
```

```@setup load_module_sequence
using DinucleotideForces
example_seq = "ATGGAGAGAATAAAAGAATTAAGGGATCTAATGTCACAGTCCCGCACTCGCGAGATACTAACAAAAACCACTGTGGACCATATGGCCATAATCAAGAAGTACACATCAGGAAGACAAGAGAAGAACCCTGCTCTCAGAATGAAATGGATGATGGCAATGAAATATCCAATCACAGCGGACAAGAGAATAACAGAGATGATTCCTGAAAGGAATGAACAAGGGCAGACGCTCTGGAGCAAGACAAATGATGCCGGATCGGACAGGTTGATGGTGTCTCCCTTAGCTGTAACTTGGTGGAATAGGAATGGGCCGACGACAAGTGCAGTCCATTATCCAAAGGTTTACAAAACATACTTTGAGAAGGCT"
```

Let's start by computing the force on the CpG motif along the full genome: this can be done with the only exported function of the package, `DimerForce`, as shown in the following code block.
```@repl load_module_sequence
motifs = ["CG"];
DimerForce(example_seq, motifs)
```
The resulting output will be a dictionary where the force is associated to each motif specified (in this case only the "CG" motif).

When no nucleotide biases (frequencies) are specified, as above, the corresponding forces are inferred together with the forces. 
We will use the word *fields* to indicate forces on nucleotides, while *forces* will be from now on restricted to dinucleotides.
Notice that $\sum_n e^{h_n} = 1$, where $n \in \{A,C,G,T\}$. This is the result of a gauge choice, because in general the model probabilities are invariant under the transformation $ h_n \to h_n+K$, which can then be used to make the fields interpretable as the logarithm of a frequency (more on this is discussed later in section [Focus on gauge choices](@ref)).

A user-specified bias can be given, and in this case fields are not inferred:
```@repl load_module_sequence
nt_bias = [0.25, 0.25, 0.25, 0.25]; # probs for A, C, G, T
motifs = ["CG"];
DimerForce(example_seq, motifs; freqs=nt_bias)
```
The package easily allows to compute forces on two or more dinucleotides:
```@repl load_module_sequence
motifs = ["CG", "TA"];
DimerForce(example_seq, motifs)
```

```@repl load_module_sequence
nts = ["A", "C", "G", "T"];
motifs = [a*b for a in nts for b in nts];
DimerForce(example_seq, motifs)
```
Look at how the CpG force changed in the cases: this is due to the fact that different motifs interact. 

### Focus on gauge choices
When the full set of fields and forces is inferred, the system of equations solved to obtain these parameters is underdetermined. This means that a gauge choice must be made, for instance by setting some forces to zero by passing less than 16 dinucleotides as motifs.
If such a choise is not made, the package selects the gauge so that:
- the exponential of the fields sum to 1;
- the forces for dinucleotides of the form NN and those of the form NT are put to zero.

Notice that the script allows for flexible choices: if not all the dinucleotides are given, the gauge is always chosen such that the maximum possible number of the nucleotides of the form NN and NT have forces equal to zero, as one can check by running the next code blocks.
```@repl load_module_sequence
motifs = ["AC", "AG", "AT", "CA", "GA", "TA"];
DimerForce(example_seq, motifs)
```

```@repl load_module_sequence
motifs = [
    "AC","AG","AT","CA",
    "CG","CT","GA","GC",
    "GT","TA","TC","TG",
    ];
DimerForce(example_seq, motifs)
```

Notice that when dinucleotides are not given, this is equivalent to fix to 0 their forces. 
Depending on which dinucleotides are not given, this might or not be a specific gauge - in other words, it might or not result in an equivalent model.
For instance, the following three cells result in equivalent models (although the third has different parameters because it is in another gauge), while the fourth does not:
```@repl load_module_sequence
motifs = ["AC", "AG", "AT", "CA", "GA", "TA"];
DimerForce(example_seq, motifs)
```
```@repl load_module_sequence
motifs = ["AC", "AG", "CA", "GA", "TA"];
DimerForce(example_seq, motifs)
```
```@repl load_module_sequence
motifs = ["AG", "AT", "CA", "GA", "TA"];
DimerForce(example_seq, motifs)
```
```@repl load_module_sequence
motifs = ["AG", "CA", "GA", "TA"];
DimerForce(example_seq, motifs)
```
In some cases (notably when all the dinucleotides are given), fixing fields is equivalent to fix (part of) a gauge, so the model is independendent on this choice (although the values of the inferred parameters depend on that).

### Inferring forces from multiple sequences
The package is able to infer fields and forces from a set of sequences:
```@repl load_module_sequence
motifs = ["CG"];
seqs = [example_seq, example_seq];
DimerForce(seqs, motifs)
```
The sequences do not need to have the same length. However, the values of the inferred parameters weakly depend on the sequence length, that has to be provided
```@repl load_module_sequence
motifs = ["CG"];
seqs = [example_seq, example_seq[1:100]];
DimerForce(seqs, motifs, 2500)
```
or it will arbitarily set by the package to the standard value of 5000
```@repl load_module_sequence
motifs = ["CG"];
seqs = [example_seq, example_seq[1:100]];
DimerForce(seqs, motifs)
```
Notice, in the last two examples, how the changes due to the different lengths used during the inference of the parameters only marginally affect their values.

### Computing energies and loglikelihoods
Given a sequence, the package allows to compute easily its energy assigned by a model
```@repl load_module_sequence
model = DimerForce(example_seq, ["CG"]);
DinucleotideForces.ComputeEnergy(example_seq, model)
```
and its log-likelihood
```@repl load_module_sequence
model = DimerForce(example_seq, ["CG"]) # hide
DinucleotideForces.ComputeLoglikelihood(example_seq, model)
```
To make sense of this large negative log-likelihood, consider that in the uniform case we have $4^L$ sequences, so the log-likelihood of each of them would be:
```@repl load_module_sequence
- length(example_seq) * log(4)
```

Notice that the function `ComputeLoglikelihood` computes the partition function each time it is runned. This can slow down the computation when the function is called many times with the same model and to prevent this the log of the partition function can be passed directly:
```@repl load_module_sequence
model = DimerForce(example_seq, ["CG", "AC"]);
ks = keys(model);
k1 = [k for k in ks if length(k)==1];
k2 = [k for k in ks if length(k)==2];
fields = [model[k] for k in k1];
forces = [model[k] for k in k2];
logZ = DinucleotideForces.eval_log_Z(fields, forces, k2, length(example_seq));
DinucleotideForces.ComputeLoglikelihood(example_seq, model; logZ=logZ)
```

### Sampling sequences
The framework used by this package introduced a probability distribution over the sequences given a set of fields and/or forces.
Sequences can be sampled from this distribution exactly (that is, without using any Markov-chain based algorithm such as Metropolis) using the simple topology of the graph of interactions between nucleotides (a tree) and the fact that the partition function can be computed exactly.
The funcion `SampleSequence` of the package implements this idea, and it is very easy to use, its only arguments being the model that has to be used to sample from and the length of the desired sequence:
```@repl load_module_sequence
nts = ["A", "C", "G", "T"];
motifs = [a*b for a in nts for b in nts];
model = DimerForce(example_seq, motifs);
DinucleotideForces.SampleSequence(model, 1000)
```

## Known issues
There are issues with the `add_pseudocount` argument of `DimerForce`. 
In principle, using `add_pseudocount=true` should add a single pseudocount to the number of each motif inferred (nucleotide or di-nucleotide).
However, from the output of the following two examples it looks that something is not working properly
```@repl load_module_sequence
nts = ["A", "C", "G", "T"];
motifs = [a*b for a in nts for b in nts];
DimerForce(example_seq, motifs)
DimerForce(example_seq, motifs, add_pseudocount=true)
```

This is even worse if a shorter sequence is used
```@repl load_module_sequence
nts = ["A", "C", "G", "T"]; # hide
motifs = [a*b for a in nts for b in nts]; # hide
ex2 = example_seq[1:100]
DimerForce(ex2, motifs)
DimerForce(ex2, motifs, add_pseudocount=true)
```

## API

```@index
```

```@autodocs
Modules = [DinucleotideForces]
```
