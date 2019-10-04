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

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Cross Cell-type Correlation in DNaseI hypersensitivity (C3D)
# Written by: Tahmid Mehdi & Paul Guilhamon
# Princess Margaret Cancer Centre - University Health Network, December 18, 2016

# Cross Cell-type Correlation in DNaseI hypersensitivity (C3D)
# Takes correlations between open regions of chromatin based on DNaseI hypersensitivity signals
# Regions with high correlations are candidates for 3D interactions
# Performs association tests on each candidate & adjusts p-values
# Produces interaction landscapes and tracks in PDF format

# an array for parameters
typeset -A config
# default values
config=(
    	[reference]=""
    	[db]=""
    	[testAnchors]=""
	[outDirectory]=""
	[matrix]=""
    	[window]="500000"
    	[correlationThreshold]="0.7"
	[pValueThreshold]="0.05"
    	[qValueThreshold]="0.1"
	[correlationMethod]="pearson"
    	[figures]="n"
	[figureWidth]="500000"
	[zoom]="0"
	[colours]="#bdd7e7,#6baed6,#3182bd,#08519c"
	[tracks]="n"
	[sampleName]=""
	[assembly]="hg19"
)
# read parameters from config file
while read line
do
	if echo $line | grep -F = &>/dev/null
    	then
		eval expanded_line="$line"
        	varname=$(echo "$line" | cut -d '=' -f 1)
        	config[$varname]=$(echo $expanded_line | cut -d '=' -f 2-)
    	fi
	if echo $line | grep -F 'module load' &>/dev/null
    	then
    		eval $line
	fi
done < $1

# if one of these parameters is empty, they'll get their default values
if [ -z "${config[window]}" ]; then
        config[window]="500000"
fi
if [ -z "${config[correlationThreshold]}" ]; then
        config[correlationThreshold]="0.5"
fi
if [ -z "${config[pValueThreshold]}" ]; then
        config[pValueThreshold]="0.05"
fi
if [ -z "${config[qValueThreshold]}" ]; then
        config[qValueThreshold]="0.05"
fi
if [ -z "${config[correlationMethod]}" ]; then
        config[correlationMethod]="pearson"
fi
if [ -z "${config[figureWidth]}" ]; then
        config[figureWidth]=${config[window]}
fi
if [ -z "${config[zoom]}" ]; then
        config[zoom]="0"
fi
if [ -z "${config[colours]}" ]; then
        config[colours]="#bdd7e7,#6baed6,#3182bd,#08519c"
fi
if [ -z "${config[tracks]}" ]; then
        config[tracks]="n"
fi
if [ -z "${config[sampleName]}" ]; then
        config[sampleName]=""
fi
if [ -z "${config[assembly]}" ]; then
        config[assembly]="hg19"
fi
# the index of the sample
trackNumber=1
# the number of total samples
numSamples=1
# overwrite parameters for multiple samples
if [ "$2" = "-ref" ]; then
        config[reference]="$3"
        config[matrix]=""
fi
if [ "$2" = "-matrix" ]; then
        config[matrix]="$3"
fi
if [ "$4" = "-out" ]; then
        config[outDirectory]="$5"
fi
if [ "$6" = "-sample" ]; then
        config[sampleName]="$7"
fi
if [ "$8" = "-track" ]; then
        trackNumber="$9"
fi
if [ "${10}" = "-numSamples" ]; then
        numSamples="${11}"
fi

# make output directory
mkdir -p ${config[outDirectory]}

# timestamp function
timestamp() {
date +"%Y-%m-%d_%H-%M-%S"
}

# if the anchor file does not have 5 fields, create a new formatted one
anchorCols=$(awk '{print NF}' ${config[testAnchors]} | sort -nu | tail -n 1)
if [ "$anchorCols" -ne "5" ]; then
        awk '{print $1"\t"$2"\t"$3"\t.\t"$1"_"$2"-"$3}' ${config[testAnchors]} > ${config[outDirectory]}/anchors_temp.bed
else
	cat ${config[testAnchors]} > ${config[outDirectory]}/anchors_temp.bed
fi

# if matrix is missing, map bedgraphs to reference
if [ "${config[matrix]}" = "" ]; then
        # Map all background files to one reference sample/peak catalogue
        echo "$(timestamp): Mapping peak files"
        cut -f 1-3 ${config[reference]} > ${config[outDirectory]}/ref.bed

        counter=1
        cat ${config[db]} | while read i; do
                echo "Mapping ${i} to file $counter.map.bed"
                mapBed -a ${config[outDirectory]}/ref.bed -b ${i} -c 4 -o max -null 0 | awk 'BEGIN{ OFS="\t" }{ print $1, $2, $3, $4 }' > ${config[outDirectory]}/$counter.map.bed
                counter=$((counter + 1))
        done
else # if matrix is given
	# create a file for the reference sample/peak catalogue from the matrix
	tail -n +2 ${config[matrix]} | awk '{print $1}' | awk -F'[:-]' '{print $1"\t"$2"\t"$3}' | sort -k1,1 -k2,2n > ${config[outDirectory]}/ref.bed
fi
# filter the anchor file for anchors within the reference catalogue
intersectBed -wa -a ${config[outDirectory]}/anchors_temp.bed -b ${config[outDirectory]}/ref.bed | sort -u -k1,1 -k2,2n > ${config[outDirectory]}/anchors.bed
config[testAnchors]="${config[outDirectory]}/anchors.bed"
rm ${config[outDirectory]}/anchors_temp.bed 
# Run R script
Rscript $DIR/c3d.R \
--refmapdir "${config[outDirectory]:-NULL}" \
--outdir "${config[outDirectory]:-NULL}" \
--anchor "${config[testAnchors]:-NULL}" \
--bg "${config[db]:-NULL}" \
--window "${config[window]:-NULL}" \
--rcut "${config[correlationThreshold]:-NULL}" \
--pcut "${config[pValueThreshold]:-NULL}" \
--qcut "${config[qValueThreshold]:-NULL}" \
--cormethod "${config[correlationMethod]:-NULL}" \
--signalmat "${config[matrix]:-NULL}" \
--figures "${config[figures]:-NULL}" \
--figwidth "${config[figureWidth]:-NULL}" \
--zoom "${config[zoom]:-NULL}" \
--colour "${config[colours]:-NULL}" \
--tracks "${config[tracks]:-NULL}" \
--sample "${config[sampleName]:-NULL}" \
--tracknum "${trackNumber:-NULL}" \
--numsample "${numSamples:-NULL}" \
--assembly "${config[assembly]:-NULL}" \
--date "$(timestamp)" \
--wdir "${DIR:-NULL}"
