#!/bin/bash

# Copyright 2016 Mathieu Lupien

# This file is part of C3D.

# C3D is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# C3D is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with C3D.  If not, see <http://www.gnu.org/licenses/>.

# Cross Cell-type Correlation in DNaseI hypersensitivity (C3D)
# Written by: Tahmid Mehdi
# Princess Margaret Cancer Centre - University Health Network, December 18, 2016

# Creates a custom track by combining outputs of C3D for a single sample or multiple samples. 
# This track can be viewed on the UCSC genome browser or IGV.

# parameters
anchor="$1"
outDirectory="$2"
references="$3"
assembly="$4"

# an array of genes from the last column of the anchor file
genes=( $(cut -d$'\t' -f5 $anchor ) )
# if references is not passed, then there's only 1 sample
if [ "$references" = "" ]; then
	# iterate through genes & make custom track for each gene
  for i in "${!genes[@]}"; do
  	cat $outDirectory/${genes[$i]}.anchor.bedGraph $outDirectory/${genes[$i]}.bedGraph > $outDirectory/${genes[$i]}.tracks.txt 
	done
else # multiple samples
	# an array of bed files for each sample
  files=( $(cut -d ' ' -f1 $references ) )
	# an array of sample names
	sample=( $(cut -d ' ' -f2 $references ) )
	
	for i in "${!sample[@]}"; do
  	refFiles[$i]="$outDirectory/${sample[$i]}/ref.bed"
  done
	bedtools multiinter -i ${refFiles[*]} | awk -v numSamples="${#files[@]}" '$4 == numSamples' | awk '{print $1"\t"$2"\t"$3}' > $outDirectory/union.bed
	echo "track name=\"Sample Intersections\" description=\" \" visibility=1 color=0,0,0 db=$assembly" | cat - $outDirectory/union.bed > $outDirectory/temp.txt && mv $outDirectory/temp.txt $outDirectory/union.bed
	
	for i in "${!genes[@]}"; do
		# for each gene, write the browser header & anchor coordinates to a file
    cat $outDirectory/${sample[0]}/${genes[$i]}.anchor.bedGraph > $outDirectory/${genes[$i]}.1.txt
    for j in "${!sample[@]}"; do
      geneFiles[$j]="$outDirectory/${sample[$j]}/${genes[$i]}.bedGraph"
    done
		# concatenate the BEDGRAPHs for each sample & the union.bed to the anchor file
    cat $outDirectory/${genes[$i]}.1.txt ${geneFiles[*]} $outDirectory/union.bed > $outDirectory/${genes[$i]}.2.txt
		
		# just copy the concatenated BEDGRAPHs to a better named file
		cp $outDirectory/${genes[$i]}.2.txt $outDirectory/${genes[$i]}.tracks.txt
		rm $outDirectory/${genes[$i]}.2.txt $outDirectory/${genes[$i]}.1.txt
	done
fi
