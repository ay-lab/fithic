# fithic

Please use the Google Group for discussions/bug reports/analysis questions:
https://groups.google.com/forum/#!forum/fithic
___________________________________________________________
Developed by Ferhat Ay, Timothy Bailey and William Noble
January 19th, 2014
___________________________________________________________________

Fit-Hi-C is a tool for assigning statistical confidence estimates
to chromosomal contact maps produced by genome architecture
assays.

___________________________________________________________________

HOW TO INSTALL DEPENDENCIES 

In order to run the fit-hi-c software you need the following to be 
present on a Unix machine (tested on LinuxMint Maya and RedHat 5).

1- Python 2.7 or higher with the libraries below installed:
  a. Scipy 
  b. Numpy
  c. rpy2 (http://rpy.sourceforge.net/rpy2_download.html)

  Once these libraries are installed you can test them by typing
  "python" in a terminal window, followed by the import statements:

	import numpy as np
	from scipy import *
	import rpy2.robjects as ro

  The output should look like this:

	Python 2.7.3 (default, Apr 20 2012, 22:39:59) 
	[GCC 4.6.3] on linux2
	Type "help", "copyright", "credits" or "license" for more information.
	>>> import numpy as np
	>>> from scipy import *
	>>> import rpy2.robjects as ro
	>>> 

___________________________________________________________________

HOW TO EXTRACT THE SOFTWARE AND SAMPLE DATA SETS

1- Extract the fit-hi-c.tgz file in a directory named fit-hi-c

	tar xzvf fit-hi-c.tgz

Once the above dependencies are installed and the extraction 
is completed, then the software is ready for use.

___________________________________________________________________


HOW TO RUN THE SOFTWARE

The tar ball contains a runall script that will run fit-hi-c on 
four sample data sets with default parameters. 

To simply run this script type

	./runall

You can also run the python script fit-hi-c.py to read usage information.

	python fithic/fithic.py -h

----	
Usage: fithic.py [options]

Options:
  -h, --help            show this help message and exit
  -f FRAGSFILE, --fragments=FRAGSFILE
                        midpoints (or start indices) of the fragments are read
                        from FRAGSFILE
  -i INTERSFILE, --interactions=INTERSFILE
                        interactions between fragment pairs are read from
                        INTERSFILE
  -o OUTDIR, --outdir=OUTDIR
                        where the output files will be written
  -t BIASFILE, --biases=BIASFILE
                        OPTIONAL: biases calculated by ICE for each locus are
                        read from BIASFILE
  -p NOOFPASSES, --passes=NOOFPASSES
                        OPTIONAL: number of passes after the initial (before)
                        fit. DEFAULT is 1 (after)
  -b NOOFBINS, --noOfBins=NOOFBINS
                        OPTIONAL: number of equal-occupancy (count) bins.
                        Default is 100
  -m MAPPABILITYTHRESHOLD, --mappabilityThres=MAPPABILITYTHRESHOLD
                        OPTIONAL: minimum number of hits per locus that has to
                        exist to call it mappable. DEFAULT is 1.
  -l LIBNAME, --lib=LIBNAME
                        OPTIONAL: Name of the library that is analyzed to be
                        used for plots.
  -U DISTUPTHRES, --upperbound=DISTUPTHRES
                        OPTIONAL: upper bound on the intra-chromosomal
                        distance range (unit: base pairs). DEFAULT no limit.
  -L DISTLOWTHRES, --lowerbound=DISTLOWTHRES
                        OPTIONAL: lower bound on the intra-chromosomal
                        distance range (unit: base pairs). DEFAULT no limit.
  -v, --visual          OPTIONAL: use this flag for generating plots. DEFAULT
                        is False.
  -q, --quiet           OPTIONAL: use this flag for omitting plots. DEFAULT
                        behavior.
  -V, --version         fit-hi-c version 1.0.1.  A tool for assigning
                        statistical confidence estimates to intra-chromosomal
                        contact maps produced by genome architecture assays.
                        Released on January 19, 2014.  Method developed by
                        Ferhat Ay, Timothy Bailey and William Noble.
                        Implemented by Ferhat Ay (ferhatay@lji.org).
                        Copyright (c), 2012, University of Washington.  This
                        software is offered under an MIT license.  For
                        details: http://opensource.org/licenses/MIT

----	

In order to run fit-hi-c with different parameter settings or on different 
data sets follow the steps below when you are in the fit-hi-c directory. 
The below example will use HindIII library from Duan et al data.

1-  Locate your input files. There are two input files you need.

	- File with list of fragments. This will be passed with -f flag. 
		-f data/fragmentLists/Duan_yeast_HindIII	

	- File with list of contact counts. This will be passed with -i flag. 	

		-i data/contactCounts/Duan_yeast_HindIII

2- Choose a genomic distance range (mid-range) for confidence estimate
   assignments. The units are always in base pairs (bp) and NOT Kb.

	- Lower bound on mid-range distances. This will be passed with
	  -L flag.  The rule of thumb here is to avoid distances lower
	  than an average meta-fragment length. When 10 consecutive RE
	  fragments are used per meta-frament use at least 50000 bp. In
	  order to have no lower bound simply don't use this argument
	  or pass -1 with the appropriate flag.  In our example here a
	  resolution of 1 RE fragment is used to process
	  Duan_yeast_HindIII. For this example we use 20000 as lower
	  bound.  Default value is -1.

		-L 20000

	- Upper bound on mid-range distances. This will be passed with
	  -U flag.  This can be disable similar to lower bound by
	  passing -1.  For this example we use 200000 as upper
	  bound. Default value is -1.

		-U 200000

3- Choose the number of steps that will be used to assign confidence
   estimates.  For instance, if 3 is selected than initial spline fit
   (spline-1) plus 2 steps of refinement of the null model will be
   applied. Results from each step will be outputted in separate
   files.  This will be passed with -p flag.  Default value is 2.

		-p 3

4- Choose the number of equal-occupancy bins that will be used for
   spline fit.  This will be passed with -b flag. Default value is
   100.

		-b 200

5- Choose a prefix that will be added in front of each file name that
   is outputted by fit-hi-c. If you want the output files to be in a
   specific directory you can first give the directory then the
   prefix. This will be passed with -l flag.  For this example we want
   output files to be under outputs directory.  Default is "" meaning
   no prefix and write in current directory.

		-l "outputs/Duan_yeast_HindIII-customSettings"

6- Run the python script fit-hi-c by the parameters selected.

	python bin/fit-hi-c.py -f data/fragmentLists/Duan_yeast_HindIII 
		-i data/contactCounts/Duan_yeast_HindIII -L 20000 -U 200000
		-p 3 -b 200 -l "outputs/Duan_yeast_HindIII-customSettings"

___________________________________________________________________


HOW TO PREPARE INPUT FILES FOR fit-hi-c

fit-hi-c requires two main input files. As long as these files are
provided fit-hi-c can be applied on data sets from different
experiments (Hi-C, 5C, ChIA-PET) either using a resolution that is
dependent on restriction fragments or with fixed sized genomic
windows.

-- The first file contains a list of fragments/windows/meta-fragments.
   Each line will have 5 entries. The second and fifth fields can be
   any integer as they are not needed in most cases. The first field
   is the chromosome name or number, the third field is the coordinate
   of the midpoint of the fragment on that chromosome, the fourth
   field is the total number of observed mid-range reads (contact
   counts) that involve the specified fragment.  The fields can be
   separated by space or tab. All possible fragments need to be listed
   in this file.  One example file would look like below (excluding
   the header which is not a part of input):

"chr	extraField	fragmentMid	marginalizedContactCount	mappable? (0/1)"
1	0		15000		234				1
1	0		25000		0				0
...


-- The second file contains a list of mid-range contacts between the
   fragments/windows/meta-fragments listed in the first file
   above. Each fragment will be represented by its chromosome and
   midpoint coordinate. Each line will have 5 fields. The first two
   will represent first fragment, the following two will represent the
   second and the fifth field will correspond to number of contacts
   between these two fragments.  The fields can be separated by space
   or tab. Only the fragment pairs with non-zero contact counts are
   listed in this file.  One example file would look like below
   (excluding the header which is not a part of input):

"chr1	fragmentMid1	chr2	fragmentMid2	contactCount"
1	15000		1	35000				23
1	15000		1	55000				12
...

___________________________________________________________________

SAMPLE DATA SETS	

This tar ball includes four sample datasets. 

1- Duan_yeast_HindIII: Aggregate of all replicates generated using
   cross-linked DNA and HindIII digestion by Duan et al. This dataset
   is processed using the natural resolution of HindIII restriction
   fragments.

2- Duan_yeast_EcoRI: Aggregate of all replicates generated using
   cross-linked DNA and EcoRI digestion by Duan et al. This dataset is
   processed using the natural resolution of EcoRI restriction
   fragments.

3- Dixon_hESC_HindIII_hg18_combineFrags10_chr1: Aggregate of all
   replicates generated using human embryonic stem cell line and
   HindIII digestion by Dixon et al. This dataset is processed using
   meta-fragments that correspond to 10 consecutive HindIII
   restriction fragments. Data for only chromosome 1 is provided due
   to large file sizes for the whole-genome.

4- Dixon_mESC_HindIII_mm9_combineFrags10_chr1: Aggregate of all
   replicates generated using mouse embryonic stem cell line and
   HindIII digestion by Dixon et al. This dataset is processed using
   meta-fragments that correspond to 10 consecutive HindIII
   restriction fragments. Data for only chromosome 1 is provided due
   to large file sizes for the whole-genome.

For more data sets or processing your own data with fit-hi-c please
contact ferhatay@lji.org

___________________________________________________________________

OUTPUT FILES AND THEIR FORMAT

Each step of fit-hi-c, the number of which is user-defined through the
-p flag, generates two output files. For step N and library name
prefix denoted by ${PREFIX} the two output files will have the
following names:

1- ${PREFIX}.fithic_passN.txt 
2- ${PREFIX}.spline_passN.significances.txt.gz


The first file will report the results of equal occupancy binning in
five fields:

"avgGenomicDist	contactProbability	standardError	noOfLocusPairs	totalOfContactCounts"
20077	2.38e-05	2.11e-06	210	19574
20228	1.88e-05	1.44e-06	268	19662
...

The second file will have the exact same lines as in the input file
that contains the list of mid-range contacts. This input file had 5
fields as described above. The output from each step will append two
more columns to these fields, namely p-value and q-value.

"chr1	fragmentMid1	chr2	fragmentMid2	contactCount	p-value	q-value"
10	100695	10	127796	11	1.000000e+00	1.000000e+00
10	104051	10	229415	12	2.544592e-02	1.202603e-01
10	104051	10	231999	15	1.506105e-03	9.463644e-03
...


COPYRIGHT
___________________________________________________________________

This software is offered under an MIT license. For details:
http://opensource.org/licenses/MIT

Copyright (c), 2012, University of Washington

Permission is hereby granted, free of charge, to any person 
obtaining a copy of this software and associated documentation 
files (the "Software"), to deal in the Software without restriction, 
including without limitation the rights to use, copy, modify, merge, 
publish, distribute, sublicense, and/or sell copies of the Software, 
and to permit persons to whom the Software is furnished to do so, 
subject to the following conditions:

The above copyright notice and this permission notice shall be 
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES 
OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND 
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS 
BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN 
ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN 
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


CONTACT
___________________________________________________________________

For any problem or request about the software, please contact
Ferhat Ay <ferhatay@lji.org> or Arya Kaul <akaul@lji.org>

Additionally, please use the FitHiC google group as a means to
communicate to other users of the software: 
https://groups.google.com/forum/#!forum/fithic


