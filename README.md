[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat-square)](http://bioconda.github.io/recipes/fithic/README.html)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/fithic/badges/version.svg)](https://anaconda.org/bioconda/fithic)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/fithic/badges/downloads.svg)](https://anaconda.org/bioconda/fithic)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/fithic/badges/license.svg)](https://anaconda.org/bioconda/fithic)

Fit-Hi-C was initially developed by Ferhat Ay, Timothy Bailey, and William Noble January 19th, 2014. It is currently maintained and updated by Ferhat Ay (ferhatay@lji<span></span>&#46;org) and Arya Kaul (akaul@lji<span></span>&#46;org) at the [Ay Lab](http://www.lji.org/faculty-research/labs/ay/#overview) in the La Jolla Institute for Allergy and Immunology.

Please use the [Google Group](https://groups.google.com/forum/#!forum/fithic) for discussions/bug reports/analysis questions. Sending an email to fithic@googlegroups<span></span>&#46;com will also post directly to the Group.

## Installation
Fit-Hi-C may be installed through one of three ways. 

 1. Bioconda
 2. Github
 3. Pip

Out of all of the following, we recommend installing through bioconda to automatically install all dependencies.

### Bioconda Installation

If this is your first time using the conda distribution system, we recommend using Miniconda as your preferred conda distribution system. This is chosen because it is the most lightweight out of all of Anaconda's distribution system; however, you're welcome to use any one you would like. More information on each may be found [here](http://bioconda.github.io/faqs.html#conda-anaconda-minconda). 

Once your conda distribution platform has been set up, you need to add the `bioconda` channel to access the bioinformatics recipes hosted there. If you're unfamiliar with `bioconda`, I highly recommend you check out the wonderful work they've done ( [link here!](http://bioconda.github.io/index.html)). To set up the `bioconda` channel, run the following:
```
conda config --add channels defaults
conda config --add channels conda-forge
conda config --add channels bioconda
```
Afterwards, simply run:
```
conda install fithic
```
This command will automatically install a command-line executable version of Fit-Hi-C along with all of its dependencies. Please run:

```
fithic -V
```
and ensure that the version number matches the version number in the `Anaconda Cloud` badge at the top of this README.

### Github Installation
Install `git` and then run:
```
git clone https://github.com/ay-lab/fithic.git
```

You will need the following dependencies installed to run Fit-Hi-C:

 - Python 3.+
 - Numpy 1.14.+
 - Scipy 1.1.+
 - Scikit-learn 0.19.+ 
 - SortedContainers 2.0.+
 - Matplotlib 2.2.+

This will create a direct clone of Fit-Hi-C in your working directory with the name fithic. You may now run Fit-Hi-Cv.2 by calling the `fithic.py` file in the `fithic/` directory:

```
python fithic.py --ARGUMENTS
```

Cloning into the repository **does not** automatically install the command line version of Fit-Hi-C, and if you desire that functionality follow the other installation instructions. 

### PyPi Installation
Ensure that you have `Python3` successfully installed on your computer.
Then run:
```
pip install fithic
```
After this is done, run `fithic --help` to ensure all necessary dependencies have been installed. Some users report that this automatically installs all dependencies, while some say it does not. If `fithic` does not work properly, install the necessary dependencies using `pip`. i.e.

```
pip install DEPENDENCY
```

## Testing
A good part of any software installation is being able to run tests on the correct installation of it. 

### Bioconda/PyPi Installation
Run the command:
```
svn export https://github.com/ay-lab/fithic/trunk/fithic/tests/
```
(If you receive an error, ensure that you have svn installed correctly.)  This command will generate a tests folder in your working directory. Going into that and running `./run_tests-cli.sh` will automatically run Fit-Hi-C on a variety of data and if everything was installed correctly you should see a final message that Fit-Hi-C executed correctly!

### Github Installation
Simply navigate into the repo and run `./fithic/tests/run_tests-git.sh`. If everything is working fine, you will see a final message that Fit-Hi-C executed correctly!

## Using Fit-Hi-C
Congratulations! If you have gotten to this point, then you have a working, fully installed version of Fit-Hi-C running on your computer. Good on you! But that was the easy part, now comes the difficult question. 

***How do I correctly use Fit-Hi-C to analyze my desired Hi-C dataset?***

Correctly answering this question requires navigating through several basic understandings:

 1. What *exactly* is Hi-C data?
 2. What does Fit-Hi-C *tell* me about this Hi-C data?
 3. *Why* is what Fit-Hi-C tells me important?

If you feel utterly comfortable with the answers to these three questions, then feel free to skip to the next section. If you are unclear about the answer to any of the above, then read on!

### What is Hi-C data?
While a beautiful, Latex type-set, easy-to-understand, and comprehensive document is being created, read [this](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0745-7) !

### What does Fit-Hi-C tell me about this Hi-C data?
At the beginning of this README, I stated:

> Fit-Hi-C is a tool for assigning statistical confidence estimates to
> chromosomal contact maps produced by genome architecture assays.

The phrase 'chromosomal contact maps produced by genome architecture assays' may be faithfully reduced to 'Hi-C data.' Applying that change yields:
    
> Fit-Hi-C is a tool for assigning statistical confidence estimates to Hi-C data.

Much less scary! From the above sentence, the only real phrase that could be misinterpreted is 'assigning statistical confidence estimates.' What does that mean? Well to find out you should read the paper Dr. Ay wrote (found [here](https://genome.cshlp.org/content/24/6/999))!


### Why is what Fit-Hi-C tells me important?
Fit-Hi-C tells you what contacts are *significant.* This is incredibly important because not all of the contacts seen in your Hi-C data are truly unexpected interactions. By assigning statistical confidence to each interaction, you will be able to determine which interactions are the most important and consequently, which ones warrant further investigation.  


## Running Fit-Hi-C
### Arguments

----------


#### Required Arguments
##### -f, --fragments
The -f argument is used to pass in a full path to what we deem a 'fragments file,' Each line will have 5 entries. The second and fifth fields can be
   any integer as they are not needed in most cases. The first field
   is the chromosome name or number, the third field is the coordinate
   of the midpoint of the fragment on that chromosome, the fourth
   field is the total number of observed mid-range reads (contact
   counts) that involve the specified fragment.  The fields can be
   separated by space or tab. All possible fragments need to be listed
   in this file.  One example file would look like below (excluding
   the header which is not a part of input):

| chr |  extraField | fragmentMid|marginalizedContactCount|mappable? (0/1)|
|---|---|---|---|---|
| 1 | 0 | 15000|234|1
|1|0|25000|0|0
|...|...|...|...|...

##### -i, --interactions
The interactions file contains a list of mid-range contacts between the
   fragments/windows/meta-fragments listed in the first file
   above. Each fragment will be represented by its chromosome and
   midpoint coordinate. Each line will have 5 fields. The first two
   will represent first fragment, the following two will represent the
   second and the fifth field will correspond to number of contacts
   between these two fragments.  The fields can be separated by space
   or tab. Only the fragment pairs with non-zero contact counts are
   listed in this file.  One example file would look like below
   (excluding the header which is not a part of input):
 
| chr1 |  fragmentMid1 | chr2|fragmentMid2|contactCount|
|---|---|---|---|---|
| 1 | 15000 | 1|35000|23
|1|15000|1|55000|12
|...|...|...|...|...


##### -o, --outdir
A full path to an output directory of your choice. If it is not already created, it creates if for you.


##### -r, --resolution
Numerical value indicating resolution of fixed-size dataset being analyzed. If non-fixed size data being studied, set -r 0.

----------


#### Optional Arguments

##### -t, --biases
*Accepts* - a fullpath to a bias file generated by ICE or Knight-Ruiz normalization for Fit-Hi-C with the following format:

| chr | midpoint |bias|
|---|---|---|
| 1 | 20000 | 1.061
| ... |... |...

*Default* - None

*Description* - Bias files help Fit-Hi-C accurately generate statistical significance estimates. If you have it, use it! 

##### -p, --passes
*Accepts* - Number of spline passes. 

*Default* - 1

*Description* - Increasing it beyond 2 is unlikely to affect Fit-Hi-C's output significantly. If you don't understand what spline fit means then *you have not read the paper!*


##### -b, --noOfBins
*Accepts* - integer representing number of equal occupancy bins you would like Fit-Hi-C to bin your data with 

*Default* - 100

*Description* - used for spline fitting


##### -m, --mappabilityThres
*Accepts* - integer representing the minimum number of hits per locus that has to exist to call it mappable

*Default* - 1

*Description* - Increasing it leads to more stringent requirements for treating an interaction as reasonable. Decreasing it leads to less stringent requirements. If you have extremely high resolution data, it may help to bump this up.

##### -l, --lib
*Accepts* - String representing prepending information for output files

*Default* - fithic

*Description* - Name of the library that is to be analyzed.

##### -U, --upperbound
*Accepts* - Integer representing upper bound for the intrachromosomal interactions to be considered in base pairs.

*Default* - -1 (no limit)

*Description* - Highly recommended to bound the intrachromosomal interactions being considered. 

##### -L, --lowerbound
*Accepts* - Integer representing lower bound for the intrachromosomal interactions to be considered in base pairs.

*Default* - -1 (no limit)

*Description* - Highly recommended to bound the intrachromosomal interactions being considered. 

##### -v, --visual
*Accepts* - no argument

*Default* - None (no plots)

*Description* - Use if plots of spline fitting are desired. Unfortunately, different matplotlib versions are unstable in different ways. If you're getting an error, I suggest trying to run Fit-Hi-C without this option and see if that helps.

##### -x, --chromosome_region
*Accepts* - 'interOnly', 'intraOnly', 'All'

*Default* - intraOnly

*Description* - interOnly is used if you would only like to analyze interchromosomal interactions. intraOnly is used if youd would only like to analyze intrachromosomal interactions. All is used if you would like to analyze inter and intrachromosomal interactions. While you may now be thinking, "Why would I ever not choose 'All'? More analysis is better!" It is not this simple. Since you are adding significantly more interactions when you analyze interchromosomal and intrachromosomal interactions in tandem, qvalues will be depressed across the board. In addition, few to no datasets are at a high enough resolution to find significanct interchromosomal interactions. 


##### -bL, --biasLowerBound
*Accepts* - float value of lower bound for bias values

*Default* - 0.5

*Description* - bias values below this number will be discarded 

##### -bU, --biasUpperBound
*Accepts* - float value of upper bound for bias values

*Default* - 2

*Description* - bias values above this number will be discarded 

----------


#### Other Arguments
##### -V, --version
*Accepts* - No arguments

*Default* - None

*Description* - Prints version number. Check to make sure this is the latest version based on version.log file here

##### -h, --help
*Accepts* - No arguments

*Default* - None

*Description* - prints help message with all options


----------


### Output
Each step of Fit-Hi-C, the number of which is user-defined through the
-p flag, generates two output files. For step N and library name
prefix denoted by ${PREFIX} the two output files will have the
following names:

 1. ${PREFIX}.fithic_passN.txt 
 2.  ${PREFIX}.spline_passN.significances.txt.gz

The first file will report the results of equal occupancy binning in
five fields. An example of which is shown below:

| avgGenomicDist| contactProb | stdErr |numLocusPairs|CCtotal|
|---|---|---|---|---|
| 20077 | 2.38e-05 |2.11e-06|210|19574
| 20228  | 1.88e-05 |1.44e-06|268|19662
|...|...|...|..|...|

The second file will have the exact same lines as in the input file
that contains the list of mid-range contacts. This input file had 5
fields as described above. The output from each step will append two
more columns to these fields, namely p-value and q-value.

| chr1 |  fragmentMid1 | chr2|fragmentMid2|contactCount|p-value|q-value|
|---|---|---|---|---|---|---|
| 1 | 15000 | 1|35000|23|1.000000e+00   |1.000000e+00
|1|15000|1|55000|12|2.544592e-02|1.202603e-01
|...|...|...|...|...|..|..|

## Utilities

These utilities are provided as part of Fit-Hi-C (`/fithic/utils/`) to aid in certain common pre-processing/post-processing steps. They are as follows:

- HiCKRy.py (Pre-processing. Generates --bias calculation)
- HiCPro2FitHiC.py (Pre-processing. Generates --interactions, --fragments, and --bias inputs)
- createFitHiCFragments-fixedsize.py (Pre-processing. Generates --fragments input)
- createFitHiCFragments-nonfixedsize.sh (Pre-processing. Generates --fragments input)
- validPairs2FitHiC-fixedSize.sh (Pre-processing. Generates --interactions input)
- createFitHiCHTMLout.sh (Post-processing. Generates HTML page describing Fit-Hi-C run)

### HiCKRy
Regardless of the implementation, we strongly recommend the use of a normalization method in order to have meaningful results for further analysis. The only way for Fit-Hi-C to utilize data from Hi-C normalization is through the bias files. As long as the bias value are scaled to have an average of 1 and high values represent loci with higher overall raw counts, Fit-Hi-C will be able to use them in significance assignment. 

HiCKRy is an in-house version of Hi-C contact map normalization using the Knight-Ruiz algorithm for fast matrix balancing. It takes three arguments:
```
-i,--interactions       Path to the interactions file to generate bias values. Required.
-f, --fragments            Path to the interactions file to generate bias values. Required.
-o, --output            Full path to output the generated bias file to. Required.
-x, --percentOfSparseToRemove     Percent of sparse low contact count loci to remove. The default value is 0.05.
```
It then outputs a bias file in the format of Fit-Hi-C's `-t` input option.

### HiCPro2FitHiC
HiC-Pro is a common Hi-C mapping tool used to extract information from the raw reads after the Hi-C assay is run. The following script enables the generation of Fit-Hi-C input directly from HiC-Pro's output.

It takes the following arguments:
```
-i MATRIX, --matrix MATRIX     Input matrix file with raw contact frequencies. Required.
-b BED, --bed BED     BED file with bins coordinates. Required.
-s BIAS, --bias BIAS     The bias file provided after IC normalization.
-o OUTPUT, --output OUTPUT     Output path.
-r RESOLUTION, --resolution RESOLUTION     Resolution of the matrix.
```
The output is the contact maps and fragments file in the format of Fit-Hi-C.

### createFitHiCFragments-fixedsize
Generates the fragments file if using a fixed-size resolution with your Hi-C data.

The script takes the following arguments:
```
--chrLens         Path to a file describing chromosome lengths of the model organism. Required.
--resolution      Resolution of dataset being studied. Required.
--outFile         Full path to the output file desired.
```

Output is a fragments file in the format of Fit-Hi-C.

### createFitHiCFragments-nonfixedsize

A bash script to generate the fragments file if the Hi-C data is RE-digested. Note, order of arguments is critical.

```
bash createFitHiCFragments-nonfixedsize.sh [outputFile] [RE] [fastaReferenceGenome]
        
[outputFile]               A desired output file path. Required.
[RE]                       Either the name of the restriction enzyme used, or the cutting position using “^”. For example, A^AGCTT for HindIII. Required.
[fastaReferenceGenome]     A reference genome in fasta format. Required.
```

### validPairs2FitHiC-fixedSize

A bash script to generate the contact maps input for Fit-Hi-C from a valid pairs file. Note, order of arguments is critical.

```
bash validPairs2FitHiC-fixedSize.sh [resolution] [libraryName] [validPairsFile]        

[resolution]         The resolution of the dataset being studied. Required.
[libraryName]        The prefix of the file generated. Required.
[validPairsFile]     A textfile containing the validPairs. Required.
```

### createFitHiCHTMLout

A bash script to generate an HTML report of the Fit-Hi-C run. Note, works best if Fit-Hi-C was run with `--visual` option.

```
bash createFitHiCHTMLout.sh [Library Name] [No. of passes] [Fit-Hi-C output folder]
        
[Library Name]            The library name (-l option) used during Fit-Hi-C’s run
[No. of passes]            The number of spline passes conducted by the Fit-Hi-C run
[Fit-Hi-C output folder]    Path to the output folder for that Fit-Hi-C run (-o option)
```

### createFitHiCContacts-hic.sh

A bash script to create Fit-Hi-C contacts from .hic files.

```
bash createFitHiCContacts-hic.sh [Juicer's dump command] [chr1] [chr2] [Output file name]
        
[Juicer's dump command]    Full path to the output of Juicer's dump command
[chr1]                     Chromosome 1 of the argument used in Juicer's dump command
[chr2]                     Chromosome 2 of the argument used in Juicer's dump command
[Output file name]          Name of output file
```

## Citing Fit-Hi-C
If Fit-Hi-C was used in your analysis, please issue the following citation:

doi: 10.1101/gr.160374.113 

or

Ferhat Ay, Timothy L. Bailey, William S. Noble. 2014. "Statistical confidence estimation for Hi-C data reveals regulatory chromatin contacts." Genome Research. 24(6):999-1011, 2014. doi: 10.1101/gr.160374.113.

## License
Copyright (c), 2012, University of Washington

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

