PCAWG7-data-preparation: prepare mutational spectra files from
VCF-like "simple files"

Copyright (C) 2018 Steven G. Rozen and Mi Ni Huang

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see
    <https://www.gnu.org/licenses/>

Contact: steverozen@gmail.com 

0. See Section 4 below for testing 

1. System requirements and software dependencies

1.1 Tested on Debian GNU/Linux 7.11

Depends on:
 Debian GNU/Linux 7.11
 Python 2.7.3
 Python package pandas 0.8.0

1.2 Please note that running on Windows machine may need to install 
    a specific version of “pynrrd” package separately. See below for 
    problem description and installation instructions:
    https://stackoverflow.com/questions/32261480/allen-brain-institute-mouse-connectivity-api-and-mouse-connectivity-cache-exam/32279226#32279226

1.3 No required non-standard hardware

2. Installation guide

No installation required - runs from source code
Install time - Nil

3. General instructions

There are three scripts that make catalogs of mutational spectra from
input ".simple" files, which are stripped-down VCF-line files.
Locations of mutations in simple files are 1-based.  Dinucleotides
("dinucs") are indicated by two lines of single nucleotide
substitutions in a .simple file for substitutions. It is necessary to
preserve the entire directory tree in the .zip file for this work, as
some inputs come from the annotation directory.

Running these scripts creates a directory called cache under the current
working directory. This is for temporary files and should be deleted
after running. It has been useful for debugging for the developers

Use --help to see documentation for each script

<path to python> make_spectra_snvs.py (--genome|--exome) --fasta <path_to_fasta> --output <outputname> <simple_file>

<path to python> make_spectra_indels.py (--genome|--exome) --fasta <path_to_fasta> -output <outputname> <simple_file>

<path to python> make_spectra_dinucs.py --output <outputname> <simple_file>

The scripts only support hg18, hg19, grch37, and grch38 human
reference genomes.

Please see demo_run/test.chr22.simple.

The output file format can also be found in the demo_run folder.

4. Demo/testing

Assuming that you have a bash shell and are starting in the main
directory in this .zip file

cd demo_run

<path to python> ../make_spectra_snvs.py --exome --fasta demo.chr22.fasta --output testout test.chr22.simple
diff -w out-96 testout-96
diff -w out-192 testout-192
diff -w out-1536 testout-1536

<path to python> ../make_spectra_indels.py --exome --fasta demo.chr22.fasta --output testout test.chr22.simple
diff -w out-indel testout-indel

<path to python> ../make_spectra_dinucs.py --output testout test.chr22.simple
diff -w out-dinuc testout-dinuc

Each script should take < 5 seconds to execute.

## Clean
rm -f testout*
rm -rf cache/

