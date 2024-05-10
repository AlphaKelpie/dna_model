# dna_model

Reproduce the code used in Tomohiro Yanao et al. paper to create a DNA structure model.

[![Standard](https://img.shields.io/badge/c%2B%2B-20-blue.svg)](https://en.wikipedia.org/wiki/C%2B%2B#Standardization)
[![Standard](https://img.shields.io/badge/Python-3.10-red.svg)](https://www.python.org)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](https://opensource.org/licenses/MIT)

## Abstract

In a paper of the 2015, Tomohiro Yanao, Sosuke Sano and Kenichi Yoshikawa have improved their previously model of the DNA double helix structure.
They have implemented this model to analyse DNA wrapping, crossover, and braiding, that are of fundamental interest on genome packaging, gene regulation, and enzyme recognition.
The study has explored elastic mechanisms for the selection of these DNA phenomena based on a coarse-grained model.
The DNA model consists of two elastic chains that mutually intertwine in a right-handed manner forming a double-stranded helix with the distinction between major and minor grooves.
Although individual potential energy functions of the DNA model have no asymmetry in terms of left and right twist, the model as a whole exhibits an asymmetric propensity to writhe in the left direction upon bending due to the right-handed helical geometry.
Monte Carlo simulations of this model suggest that DNA has a propensity to prefer left-handed wrapping around a spherical core particle and also around a uniform rod due to the asymmetric elastic coupling between bending and writhing.
This result indicates an elastic origin of the uniform left-handed wrapping of DNA in nucleosomes and also has implications on the wrapping of double-stranded DNA around rod- like molecules.
Monte Carlo simulations of the DNA model also suggest that two juxtaposed DNA molecules can braid each other spontaneously under moderate attractive interactions preferring left-handed braiding due to the asymmetric coupling between bending and writhing.
This result suggests the importance of asymmetric elasticity in the selection of chirality in braiding of a pair of DNA molecules.
I have recreated the simulations to check their results.
The main aspects of the study have been achieved: the DNA model prefers left-handed writhing and braiding and the distributions of the parameters reflect the expected ones.
Some differences happen, maybe due to statistical issues and to a complex landscape potential.

To see my report and all the all quotations, see ```report.pdf```

## Code

### Simulations

In ```src``` directory there is the code.

The program to make simulations is in ```c++ 20``` and different cpp files performe different simulations:

- ```main.cpp``` is the base.
  It creates a DNA molecule, calculates energy and saves coordinates;
- ```bend_writhe.cpp``` simulates different constant bending and dihedral angles (section 3.1);
- ```histone.cpp``` simulates wrapping around an histone (section 3.2);
- ```rod_1.cpp``` simulates wrapping around a rod with local interaction (first part of section 3.3);
- ```rod_2.cpp``` simulates wrapping around a rod with global interaction (second part of section 3.3);
- ```brand.cpp``` simulates braiding and crossover of two DNA molecules (section 3.4 and 3.5).

The header files regard:

- ```functions.hpp```: all the functions used in all the simulations;
- ```base.hpp```: class of a base-pairs. It represents the backbone structure and the two sugar-phosphate chains;
- ```coordinates.hpp```: structure of the coordinates used for backbone structure and the two sugar-phosphate chains;
- ```outputs.hpp```: structure used for energy, wrapping and chirality;

```tqdm.hpp``` and ```èigen-3.4.0``` are imported libraries, so thank you for [tqdm](https://github.com/mraggi/tqdm-cpp) and [eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page).

To compile the simulations:

```bash
g++ -std=c++2a -I ./eigen-3.4.0 <file_name>
```

### Perform simulations

I have written some python scripts to perform the repeated simulation parallelized.
Their names are the same of the cpp file, but for ```brand.cpp``` there are two: ```brand_F.py``` for section 3.4 and ```brand_D.py``` for section 3.5.
You need ```subprocess```, ```datetime```, and ```multiprocessing```.
```send_email.py``` is a module that I have written: it send you an email when the code starts and when the code ends.
You must add your gmail account with a "app password" (see [italian](https://support.google.com/mail/answer/185833?hl=it-419) and [english](https://support.google.com/accounts/answer/185833?hl=en) Google support pages), and the recipient.

### Figures

All the python scripts that end with ```*_fig.py``` creates figures from the outputs of the simulations.
Check that all the folders exist and are correct.

## Citation

If you find the code in this repository useful and you use it for any purpose, please consider citing it

```BibTeX
@misc{dna_model,
  author = {Rondini, Tommaso},
  title = {DNA coarse-grain elastic model to reproduce characteristics of higher-order structures},
  year = {2024},
  publisher = {GitHub},
  howpublished = {\url{https://github.com/AlphaKelpie/dna_model}}
}
```
