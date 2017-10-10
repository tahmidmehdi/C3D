# C3D

**C**ross **C**ell-type **C**orrelation in **D**NaseI hypersensitivity

Developed by [Tahmid Mehdi](https://github.com/tahmidmehdi)

## Instructions

C3D calculates correlations between open regions of chromatin based on DNase I hypersensitivity signals. Regions with high correlations are candidates for 3D interactions. It also performs association tests on each candidate and adjusts p-values.

At minimum, C3D requires:

1. BED file of accessible genomic regions
1. BED file of potential chromatin loop anchors (for example, promoters of genes you are interested in)
1. BEDGRAPH files of DNase I hypersensitivity signals for each biological sample you want to calculate correlations across
1. An output directory for results

The file paths and names of these files and directories must be provided in a configuration file. Many optional parameters can also be set in the file.

For each anchor, C3D will compute correlations between open regions, that overlap the anchor, and regions distal to it, based on DNase I hypersensitivity signals for the provided biological samples. There is an option to calculate correlations between each anchor region and all other open regions genome-wide.

Finally, each correlation is tested for significance and multiple test corrections are performed. p-values and q-values are reported.

## Usage

Usage: `sh c3d <CONFIG>`

### Configuration file parameters

| Parameter | Description |
|-----------|-------------|
| `testAnchors` | BED file of anchor regions. |
| `outDirectory` | Output directory for results. |
| `reference` | **Optional if `matrix` provided**. BED file of accessible genomic regions. |
| `db` | **Optional if `matrix` provided**. Text file containing a list of full paths and names of BEDGRAPH files for calculating correlations across.
| `matrix` | **Optional if `reference` and `db` provided**. Tab-delimited text file of signal data with regions as rows and samples as columns. It must have row names in `chr:start-end` format. If not provided, a matrix will be generated using `reference` and `db`. If you want to run C3D multiple times on the same `reference` and `db`, pass `reference` and `db` on the first run. For any subsequent runs, pass `matrix` with the `signalMatrix.txt` file created during the first run. |

### Configuration file options

| Option | Description |
|--------|-------------|
| `window` | The maximum number of bps from an anchor an open region can be to be considered distal. Correlations are only calculated between distal regions unless window is set to genome. If set to genome, it will compute correlations between each anchor region and every region in the reference file. Defaults to 500000. |
| `correlationThreshold` | Only interaction candidates with correlations above or equal to this will be shown in figures. Defaults to 0.7. |
| `pValueThreshold` | Only interaction candidates with p-values below or equal to this will be shown. Defaults to 0.05. |
| `qValueThreshold` | Only interaction candidates with q-values below or equal to this will be shown. Defaults to 0.1. |
| `correlationMethod` | Correlation coefficient (`pearson`, `spearman`, `kendall`). Defaults to `pearson`. |
| `figures` | `y` if an arc diagram should be generated for each anchor. Figures will be shown in the same order as anchors from the anchor file. Leave this parameter out if you do not want figures or have many anchors. |
| `figureWidth` | Number of flanking bps from anchor to display on figures. Use a comma-separated list to assign different widths to different figures. For example, if you have 3 figures and `figureWidth=500000,100000,400` in the configuration file, then the figures will show 500000, 100000, and 400 bps around the first, second, and third anchors respectively. Defaults to `window`. |
| `zoom` | Number of flanking bps from anchor to show for zoomed figures. You can also pass a comma-separated list of numbers if you want different zoom regions for each figure. For example, if you have 3 figures and you pass `10000,50000,75000`, then the zoom regions for figures 1,2, and 3 will be 10000,50000, and 75000 respectively. Use 0 if you do not wish to zoom. Defaults to no zoom. |
| `colours` | 4 colours for the loops in the arc diagrams. Loops are colour-coded based on q-values. For example, if blue,red,green,purple is passed then loops of interactions with q-values>0.05 will be blue, between 0.05 and 0.01 will be red, between 0.01 and 0.001 will be green and below 0.001 will be purple. Hexadecimal digits accepted. Defaults to shades of blue. |
| `tracks` | `y` if files with tracks should be generated. These files can be uploaded to the UCSC Genome Browser or the Integrative Genomics Viewer (IGV) to view interaction candidates. Leave this parameter out if you do not want figures or have many anchors. |
| `assembly` | Reference genome for tracks. Defaults to `hg19`. |

#### Recommended testAnchors file format

Do not include the header in the `testAnchors` file.

```shell
CHR   START     END       STRAND ID
chr12 103348963 103351963 +      ASCL1_ENST00000266744.3
chr7  5567229   5570229   -      ACTB_ENST00000464611.1
chr7  155604467 155607467 .      SHH_ENST00000297261.2
```

### Configuration file options for multiple samples

If you want to run C3D on multiple samples, do not include `reference` or `matrix`. Instead, include one of the following:

| Option | Description |
|--------|-------------|
| `references` | **Optional if `matrices` provided**. A file containing a list of reference BED files where each line is of the form `<reference.bed> <sample name>`. |
| `matrices` | **Optional if `references` and `db` provided**. A file containing a list of text files with matrices where each line is of the form `<signalMatrix.txt> <sample name>`. |

All other parameters are still available.

## Output

The results file is named `results_<timestamp>.txt` and the results are formatted as a table with the following columns:

| Column Name | Description |
|-------------|-------------|
| `COORD_1` | Genomic coordinates of open region in an anchor |
| `COORD_2` | Genomic coordinates of open region that is distal to an anchor |
| `R_<correlationMethod>` | Correlation score between `COORD_1` and `COORD_2` |
| `ID` | ID of anchor that COORD_1 overlaps |
| `p_Value` | P-value from association test between `COORD_1` and `COORD_2` |
| `q_Value` | Adjusted p-value after performing multiple test correction |

Arc diagrams can be viewed in `figures_<timestamp>.pdf` (if `figures=y` is in the configuration file).
Tracks can be found in files called `<anchor ID>.tracks.txt`, `<anchor ID>.anchor.bedGraph`, and `<anchor ID>.bedGraph` for each anchor (if `tracks=y` is in the configuration file).

## Installation

### Environment

* [Bash](http://www.gnu.org/software/bash/manual/html_node/Installing-Bash.html) (version >= 4.0)
* [R](http://cran.r-project.org) (version >= 3.3.1)
* [bedtools](http://bedtools.readthedocs.io/en/latest/content/installation.html) (version >= 2.19.0)

#### R packages

* GenomicRanges
* Sushi
* data.table
* preprocessCore
* dynamicTreeCut

To install them, start R and enter the following commands:

```R
source("https://bioconductor.org/biocLite.R")
biocLite(c("GenomicRanges", "Sushi", "preprocessCore"))
install.packages(c("data.table", "dynamicTreeCut"))
```

### Download C3D and example data

C3D and its documentation are available in this repo. You can download them via `git clone`.
Set the permissions for the 3 files that start with `c3d` to executable. i.e. `chmod +x c3d*`.

For this example, we will predict genomic interactions involving the promoters of the GATA2 gene in K562 (chronic myelogenous leukemia) cells by calculating correlations across DNase-seq signals for 79 ENCODE cell lines.

All of the data for this example is located in the `example` folder.

1. The BED file `test_anchors.bed` contains anchors for GATA2. This file was generated by

```shell
Rscript make_anchors_from_genes.R genes.txt 1000 0 test_anchors.bed
```

2. `signalMatrix.txt` is formatted as a matrix whose rows correspond to DNase I Hypersensitive Sites (DHSs) and columns correspond to cell lines. Each entry is the maximum DNase-seq signal over the given DHS in the given cell line. This result is generated by the remaining steps.
3. Download a BEDGRAPH file for each of the 79 cell lines: `sh download_79_ENCODE_cell_lines.sh`
4. Download 126 BEDGRAPH files from ENCODE and 53 BEDGRAPH files from Roadmap for the db (requires [UCSC command line bioinformatic utilities](https://github.com/ENCODE-DCC/kentUtils)):

```shell
sh download_126_encode_bedgraphs.sh
sh download_convert_53_roadmap_bigwigs.sh
```

5. Create a file that lists all of the full paths and file names of the BEDGRAPH files. Each file name should be on a separate line.
6. Download a BED file of DHSs in K562

```shell
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwDnase/wgEncodeUwDnaseK562PkRep1.narrowPeak.gz
gunzip wgEncodeUwDnaseK562PkRep1.narrowPeak.gz
```

7. Edit `config.txt` by removing the `matrix=` line, and adding

```shell
reference=/path/to/wgEncodeUwDnaseK562PkRep1.narrowPeak
db=/path/to/list_of_bedgraphs.txt
```

### Configure and run C3D

Edit `config.txt` by setting `testAnchors`, `outDirectory`, and `matrix` with appropriate paths and file names. Run C3D with via `sh c3d config.txt`

C3D will output the results to `results_<timestamp>.txt`. Arc diagrams for each anchor can be found in `figures_<timestamp>.pdf`. Tracks will be generated in files ending in `.tracks.txt` and `.bedGraph`.

Files ending with `.tracks.txt` can be uploaded to the [UCSC Genome Browser](https://genome.ucsc.edu/cgi-bin/hgGateway) to visualize the correlations. Alternatively, tracks ending in `.bedGraph` can be uploaded to IGV.

Instructions for downloading and using IGV can be found [here](https://software.broadinstitute.org/software/igv/download). All of these files will be located in the specified outDirectory.
