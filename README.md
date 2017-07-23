# PASTRI
------

Version 0.1

PASTRI is an algorithm that infers tumor phylogenies from one or more bulk DNA sequencing samples.

If you use this software in your research, please cite

``
Satas, G., & Raphael, B. J. (2017). Tumor phylogeny inference using tree-constrained importance sampling. Bioinformatics, 33(14), i152-i160.
``

# Table of Contents

1. License
2. Dependencies
3. Running PASTRI  
    3.1 Usage  
    3.2 File Format  
    3.3 Input Files  
    3.4 Output Files  
4. Basic Examples  
5. Development 
6. Support
    

## 1 License
----

```
Copyright 2017 Brown University, Providence, RI.

                         All Rights Reserved

Permission to use, copy, modify, and distribute this software and its  
documentation for any purpose other than its incorporation into a  
commercial product is hereby granted without fee, provided that the  
above copyright notice appear in all copies and that both that  
copyright notice and this permission notice appear in supporting  
documentation, and that the name of Brown University not be used in  
advertising or publicity pertaining to distribution of the software  
without specific, written prior permission.  

BROWN UNIVERSITY DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS SOFTWARE,  
INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR ANY  
PARTICULAR PURPOSE.  IN NO EVENT SHALL BROWN UNIVERSITY BE LIABLE FOR  
ANY SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES  
WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN    
ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF  
OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.  
```

## 2 Dependencies
---

We recommend the [Anaconda distribution](https://www.continuum.io/downloads) of Python 2 which contains all required 
dependencies. PASTRI was tested with Anaconda 4.4.0.

Dependencies:

1. Python 2.x
2. numpy
3. scipy


## 3 Running PASTRI
-----

### 3.1 Usage
```
python src/RunPASTRI.py [-h] [-n NUM_ITERS] [-o OUTPUT_PREFIX]
                    path/to/data_file path/to/proposal_file

positional arguments:
  data_file
  proposal_file

optional arguments:
  -h, --help            show help message and exit
  -n NUM_ITERS, --num_iters NUM_ITERS
  -o OUTPUT_PREFIX, --output_prefix OUTPUT_PREFIX
```

### 3.2 File Format

All files are organized as a series of lists or matrices, separated by a blank line
The format features a header giving the name or a description of the component, 
followed by the shape of the matrix, and then the matrix in tab separated format.

```
> [Name]
(# of columns, # of rows)
x_1,1    x_1,2  ...   x_1,c
|        |            | 
x_r,1    x_r,2  ...   x_r,c
```

### 3.3 Input Files

1. Allele Counts File

See `example/example.input` for an example. 
Matrix A is the variant read count matrix. Each row is an SNV, and each column is a sample.
Matrix D is the total (variant + reference) read count matrix. Each row is an SNV and each column is a sample. 

2. Proposal Distribution File

See `example/example.proposal` for an example.
Matrices Alpha and Beta correspond to parameters for a beta distribution, where each row corresponds to a
cluster of SNVs and each column corresponds to a sample. 


### 3.4 Output Files

1. Tree Posterior

See `example/example.trees` for example, after running basic example in section 3. 
Each matrix correponds to a tree topology. The name is formatted as:

```
 > rank:id:Log-likelihood
```

The provided matrix is in perfect phylogeny format. 


## 4 Basic Example
-----

An example input file is provided in `example/` directory. This example uses tree with 5 samples, 20 mutations, and 8 clusters. 
Clone PASTRI repository to your local machine. In the repository run 

```
python src/RunPASTRI.py example/example.input example/example.proposal -o example/
```

This will run PASTRI on a basic example with 20 mutations and 5 samples, with 8 clusters of mutations. 
PASTRI will execute 1000 iterations and then report the posterior distributions over trees in a `example/example.trees`. 


## 5 Development
-----

Follow PASTRI development on our [Trello board](https://trello.com/b/4AKPd5GN/pastri-development) to see in progress and upcoming features. 

## 6 Support
-----

For support, please open an issue on the GithHub page or email gryte_satas (at) brown (dot) edu.