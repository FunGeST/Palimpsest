<p align="center">
  <a href="https://github.com/FunGeST/Palimpsest">
    <img height="150" src="https://github.com/FunGeST/Palimpsest/blob/master/Files/Palimpsest.jpg">
  </a>
  <h1 align="center">Palimpsest</h1>
</p>

Cancer genomes are altered by various mutational processes and, like palimpsests,
bear the signatures of these successive processes. The Palimpsest R package provides a complete workflow for the characterization and visualization of mutational signatures and their evolution along tumor development. The package covers a wide range of functions for extracting both base substitution and structural variant signatures, inferring the clonality of each alteration and analyzing the evolution of mutational processes between early clonal and late subclonal events. Palimpsest also estimates the probability of each mutation being due to each process to predict the mechanisms at the
origin of driver events. Palimpsest is an easy-to-use toolset for reconstructing the natural history of a tumor using whole exome or whole genome sequencing data.

Installation
========
Install from the GitHub repository using devtools:

    install.packages("devtools")
    library(devtools)
    devtools::install_github("FunGeST/Palimpsest")

Dependencies
========
Recommended to install R package "bedr" in order to perform structural variant signature analysis. The bedr API gives access to "BEDTools" as well as offers additional utilities for genomic regions processing. To gain the functionality of bedr package you will need to have the [BEDTools](http://bedtools.readthedocs.io/en/latest/content/installation.html) program installed and in your default PATH.

Input files
========
The input to Palimpsest are two main files -- one for **mutational catalogue** data and
one for **copy number information** and a thrid supplementary file with **sample annotation** for gender and tumor purity. 

**The input files should have the following columns (the header is required, but the order of the columns can change). For proper format, please refer the example files provided with the package.**

`1]. mut_data`: __simple somatic substitutions data__

* `Sample`: Sample identifier. Any alphanumeric string.
* `Type`: Mutation type [SNV/Indel].
* `CHROM`: Chromosome. Between chr1 and chr22 or the chrX or chrY ('chr' prefix required).
* `POS`: Mutation position. A positive integer.
* `REF`: Reference base(s): Each base must be one of A,C,G,T (upper case). Multiple bases are permitted. The value in the POS field refers to the position of the first base in the string.
* `ALT`: Alternate base(s): Each base must be one of A,C,G,T (upper case). Multiple bases are permitted. The value in the POS field refers to the position of the first base in the string.
* `Tumor_Varcount`: Number of variant bases at the position in the tumor sample.
* `Tumor_Depth`: Tumor sample sequencing depth at the position.
* `Normal_Depth`: Normal sample sequencing depth at the position.
* `Gene_Name`: OPTIONAL column for representing mutated gene name.
* `Driver`: OPTIONAL column indicating the driver events to be annotated in oncogenic timeline plots.

`2]. cna_data`: __copy number alteration data__

* `Sample`: Sample identifier. Any alphanumeric string.
* `CHROM`: Chromosome. Between chr1 and chr22 or the chrX or chrY ('chr' prefix required).
* `POS_START`: Start position of segmented chromosome.
* `POS_END`: End position of segmented chromosome.
* `LogR`: LogR information.
* `Nmin`: Minor allele copy number.
* `Nmaj`: Major allele copy number.
* `Ploidy`: Tumor ploidy.

`3]. annot_data`: __sample annotation data__

* `Sample`: Sample identifier. Any alphanumeric string.
* `Gender`: Gender information for patient [M/F].
* `Purity`: Tumor purity estimate (Represented as fraction; ranging between 0.01 - 1).

**Optional:**

**__Structural variant signature analysis__**

`4]. sv_data`: __structural variant data__

* `Sample`: Sample identifier. Any alphanumeric string.
* `Type`: Type of structural variant: [INV/DEL//DUP/BND](https://samtools.github.io/hts-specs/VCFv4.1.pdf).
* `CHROM_1`: Chromosome of the first breakpoint. Between chr1 and chr22 or the chrX or chrY ('chr' prefix required).
* `POS_1`: Position of the first breakpoint. A positive integer.
* `CHROM_2`: Chromosome of the second breakpoint. Between chr1 and chr22 or the chrX or chrY ('chr' prefix required).
* `POS_2`: Position of the second breakpoint. A positive integer.
* `Tumor_Varcount`: Column for variant allele count information.
* `Tumor_Depth`: Column for tumor sequencing depth information.
* `Normal_Depth`: Column for normal sequencing depth information.
* `Driver`: OPTIONAL column indicating the driver events to be annotated in oncogenic timeline plots.


Running Palimpsest
================

* The RUNNING_PALIMPSEST_EXAMPLE folder contains example data and an R script for a typical Palimpsest analysis. Please try!</br>
* [*Introduction to Palimpsest*](http://nbviewer.jupyter.org/github/FunGeST/Palimpsest/blob/master/Files/vignette_palimpsest.pdf) provides a comprehensive example of the Palimpsest workflow with detailed  explanations of each function.</br> 
* Please refer to the following paper for extensive description of the statistical methods used in the package: [Mutational signatures reveal the dynamic interplay of risk factors and cellular processes during liver tumorigenesis](https://www.nature.com/articles/s41467-017-01358-x), with authors Letouzé, E#, Shinde, J.#, Renault, V., Couchy, G., Blanc, J.-F., Tubacher, E., Bayard, Q., Bacq, D., Meyer, V., Semhoun, J., Bioulac-Sage, P., Prévôt, S., Azoulay, D., Paradis, V., Imbeaud, S., Deleuze, J.-F. and Zucman-Rossi, J. 

<p align="center">
  <a href="https://github.com/FunGeST/Palimpsest">
    <img height="550" src="https://github.com/FunGeST/Palimpsest/blob/master/Files/RUNNING_PALIMPSEST.png">
  </a>
</p>

Fig.(A) Workflow illustrating a typical analysis with Palimpsest. Taking as input somatic mutations, copy-number alterations (CNAs) and structural variants, the package classifies variants as clonal and subclonal, extracts mutational and structural variant
signatures separately in early clonal and late subclonal events, and estimates the probability of each alteration being due to each process. The timing of chromosome duplications is also estimated from the ratio of duplicated/non-duplicated mutations to reconstruct the complete natural history of the tumor. (B) Example of output representing, for one tumor, the number of clonal and subclonal mutations, their distribution per mutation signature, the driver alterations (colored according to the most likely causal mutational process) and CNA timing.

### The package is still under development. Will be updated soon.


License
========

Copyright (C) 2017 Jayendra Shinde

Palimpsest is a free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
