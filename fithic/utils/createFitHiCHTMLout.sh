#!/bin/bash
#set -o nounset
#set -o pipefail
#set -o errexit

#DESCRIPTION:
#bash createFitHiCHTMLout.sh [Library Name] [No. of passes] [Fit-Hi-C output folder]
#        
#        [Library Name]            The library name (-l option) used during Fit-Hi-C’s run
#        [No. of passes]            The number of spline passes conducted by the Fit-Hi-C run
#        [Fit-Hi-C output folder]    Path to the output folder for that Fit-Hi-C run (-o option)



libName=$1 # name of the library provided to fithic python script 
noOfPasses=$2 # assumes this is 1 unless otherwise specified
output=$3 #path to directory of output folder from Fit-Hi-C run
html=$output/$libName.results.html
echo \<html\> > $html
echo \<head\>\<title\> Fit-Hi-C Report for $libName \</title\>\</head\> >>$html
echo \<body\> >> $html
echo \<center\>\<h1\>Fit-Hi-C Report for $libName \</h1\>\</center\>  >>$html

echo \<p\> Displaying results of the Fit-Hi-C run for $libName.\</p\> >>$html


echo \<p\> For detailed log of the run please see this \<a href=$output/$libName.fithic.log\> log file\</a\> \</p\> >>$html

echo \<p\> Final significance values \(p-values, q-values\) after specified number of refinement steps can be downloaded from \<a href=$output/$libName.spline_pass$noOfPasses.significances.txt.gz\> this link\</a\> \</p\> >>$html

echo \</br\> >>$html
if [[ $noOfPasses -gt 1 ]]; then
	echo \<p\> \<h3\> Comparison of spline fits and number of significant contacts from the first and the last spline fit iteration. \</h3\>\</p\> >>$html

	echo \<table  border=\"0\"\> >> $html
	echo \<tbody\> >> $html
	echo \<tr\> >> $html
	echo \<th scope=\"col\"\> Spline fits \</th\> >> $html
	echo \<th scope=\"col\"\> Number of significant contacts \</th\> >> $html
	echo \</tr\> >> $html

	echo \<tr\> >> $html
	png=$output/$libName.spline_comparison.png
    img=$( base64 $png )
	echo "<td> <a href=\"$png\"><img src=\"data:image/png;base64,$img\" width=500> </a></td>" >> $html
	#echo \<td\> \<a href=\"$png\"\>\<img src=\"$png\" width=500\> \</a\>\</td\> >> $html
	png=$output/$libName.spline_FDR_comparison.png
    img=$( base64 $png )
	echo "<td> <a href=\"$png\"><img src=\"data:image/png;base64,$img\" width=500> </a></td>" >> $html
	echo \</tr\> >> $html
	echo \</tbody\>\</table\> >> $html
fi

echo \</br\> >>$html
echo \<p\> \<h3\> Relevant figures and text files listed separately for each pass of the spline fit. \</h3\>\</p\> >>$html
echo \<table  border=\"1\"\> >> $html
echo \<tbody\> >> $html
echo \<tr\> >> $html
echo \<th scope=\"col\"\> Spline Pass \</th\> >> $html
echo \<th scope=\"col\"\> Spline Fit \</th\> >> $html
echo \<th scope=\"col\"\> Number of significant contacts \</th\> >> $html
echo \</tr\> >> $html
for p in `seq 1 $noOfPasses`; do

	echo \<tr\> >> $html
	echo \<td\> Pass $p \</td\> >> $html
	png=$output/$libName.spline_pass$p.png
    img=$( base64 $output/$libName.spline_pass$p.png )
	echo "<td> <a href=\"$png\"><img src=\"data:image/png;base64,$img\" width=450> </a></td>" >> $html
	png=$output/$libName.spline_pass$p.qplot.png
    img=$( base64 $output/$libName.spline_pass$p.qplot.png )
	echo "<td> <a href="$png"><img src=\"data:image/png;base64,$img\" width=450> </a></td>" >> $html
	echo \</tr\> >> $html
	
	echo \<tr\> >> $html
	echo \<td\> Pass $p \</td\> >> $html
	txt=$output/$libName.fithic_pass$p.txt
	echo \<td\> \<a href=\"$txt\"\> Equal occupancy bin statistics \</a\>\</td\> >> $html
	png=$output/$libName.spline_pass$p.significances.txt.gz
	echo \<td\> \<a href=\"$png\"\> Significance estimates \(p-, q-values\) [gzipped] \</a\>\</td\> >> $html
	echo \</tr\> >> $html
	
	echo \<tr\> \</tr\> >> $html
	echo \<tr\> \</tr\> >> $html
done
 echo \</tbody\>\</table\> >> $html
 echo \</body\>\</html\> >> $html


