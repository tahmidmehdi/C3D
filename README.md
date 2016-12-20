# C3D
Cross Cell-type Correlation in DNaseI hypersensitivity

INSTRUCTIONS

C3D - Cross Cell-type Correlation in DNase I hypersensitivity

calculates correlations between open regions of chromatin based on DNase I hypersensitivity signals. Regions with high correlations are candidates for 3D interactions. It also performs association tests on each candidate and adjusts p-values. At minimum, C3D requires a BED file of accessible genomic regions, a BED file of potential chromatin loop anchors (for example, promoters of genes you are interested in), BEDGRAPH files of DNase I hypersensitivity signals for each biological sample you want to calculate correlations across, and an output directory for results. The file paths and names of these files and directories must be provided in a configuration file. In addition, many optional parameters can also be set in the file. For each anchor, C3D will compute correlations between open regions, that overlie the anchor, and regions distal to it, based on DNase I hypersensitivity signals for the provided biological samples. There is also an option to calculate correlations between each anchor region and all other open regions genome-wide. Then, each correlation is tested for significance and multiple test corrections are performed; p-values and q-values are reported.

RUNNING C3D

Usage: sh c3d < config file >

Parameters for config file

	reference (Optional if matrix provided): A BED file of accessible genomic regions.

	anchor (Required): A BED file of anchor regions.

		The recommended format of the anchor file is as follows (Do not include the header):

		CHR	START	END	STRAND	ID

		chr12   103348963       103351963       +       ASCL1_ENST00000266744.3
		chr7    5567229 5570229 -       ACTB_ENST00000464611.1
		chr7    155604467       155607467       .       SHH_ENST00000297261.2

	db (Optional if matrix provided): A TEXT file containing a list of full paths and names of BEDGRAPH files for calculating correlations across.

	outDirectory (Required): Output directory for results.

	matrix (Optional if reference and db provided): A tab-delimited TEXT file of signal data with regions as rows and samples as columns. It must have row names in chr:start-end format. If not provided, a matrix will be generated using reference and db. If you want to run C3D multiple times on the same reference and db, pass reference and db on the first run. For any subsequent runs, pass matrix with the signalMatrix.txt file created during the first run. 

	window (Optional): The maximum number of bps from an anchor an open region can be to be considered distal. Correlations are only calculated between distal regions unless window is set to genome. If set to genome, it will compute correlations between each anchor region and every region in the reference file. Defaults to 500000.

	correlationThreshold (Optional): Only interaction candidates with correlations above or equal to this will be shown in figures. Defaults to 0.7.

	pValueThreshold (Optional): Only interaction candidates with p-values below or equal to this will be shown. Defaults to 0.05.

	qValueThreshold (Optional): Only interaction candidates with q-values below or equal to this will be shown. Defaults to 0.1.

	correlationMethod (Optional): Correlation coefficient (pearson, spearman, kendall). Defaults to pearson.

	figures (Optional): y if an arc diagram should be generated for each anchor. Figures will be shown in the same order as anchors from the anchor file. Leave this parameter out if you do not want figures or have many anchors.

	figureWidth (Optional): Number of flanking bps from anchor to display on figures. Use a comma-separated list to assign different widths to different figures. For example, if you have 3 figures and 500000,100000,400 is passed then the figures will show 500000, 100000, and 400 bps around the first, second, and third anchors respectively. Defaults to window.

	zoom (Optional): Number of flanking bps from anchor to show for zoomed figures. You can also pass a comma-separated list of numbers if you want different zoom regions for each figure. For example, if you have 3 figures and you pass 10000,50000,75000 then the zoom regions for figures 1,2, and 3 will be 10000,50000, and 75000 respectively. Use 0 if you do not wish to zoom. Defaults to no zoom.

	colours (Optional): 4 colours for the loops in the arc diagrams. Loops are colour-coded based on q-values. For example, if blue,red,green,purple is passed then loops of interactions with q-values>0.05 will be blue, between 0.05 and 0.01 will be red, between 0.01 and 0.001 will be green and below 0.001 will be purple. Hexadecimal digits accepted. Defaults to shades of blue.

	tracks (Optional): y if files with tracks should be generated. These files can be uploaded to the UCSC Genome Browser or the Integrative Genomics Viewer (IGV) to view interaction candidates. Leave this parameter out if you do not want figures or have many anchors.

	assembly (Optional): Reference genome for tracks. Defaults to hg19.
Multi-sample Parameters

If you want to run C3D on multiple samples, do not include reference or matrix. Instead, include one of the following:

	references (Optional if matrices provided): A file containing a list of reference BED files where each line is of the form: 

		<reference.bed> <sample name>

	matrices (Optional if references and db provided): A file containing a list of TEXT files with matrices where each line is of the form: 

		<signalMatrix.txt> <sample name>

All other parameters are still available.

OUTPUT

The results file is named results_< timestamp >.txt and the results are formatted as a table.

	COORD_1	Genomic coordinates of open region in an anchor	
	COORD_2	Genomic coordinates of open region that is distal to an anchor
	R_<correlationMethod>	Correlation score between COORD_1 and COORD_2
	ID	ID of anchor that COORD_1 overlaps
	p_Value	P-value from association test between COORD_1 and COORD_2
	q_Value	Adjusted p-value after performing multiple test correction

Arc diagrams can be viewed in a file called figures_< timestamp >.pdf (if figures=y is in config file).
Tracks can be found in files called < anchor ID >.tracks.txt, < anchor ID >.anchor.bedGraph and < anchor ID >.bedGraph for each anchor (if tracks=y is in config file).

TUTORIAL (CASE EXAMPLES)

Setting up the software environment

	Bash (version >= 4.0), R (version >= 3.3.1) and bedtools (version >= 2.19.0) must be installed on the system to run C3D properly. Instructions for downloading Bash, R and bedtools can be found from http://www.gnu.org/software/bash/manual/html_node/Installing-Bash.html, http://cran.r-project.org and http://bedtools.readthedocs.io/en/latest/content/installation.html, respectively.

	The following R packages must be installed: GenomicRanges, Sushi, data.table, preprocessCore and dynamicTreeCut. To install them, start R and enter the following commands:
		source("https://bioconductor.org/biocLite.R")
		biocLite(c("GenomicRanges", "Sushi", "preprocessCore"))
		install.packages(c("data.table", "dynamicTreeCut"))
		
Download C3D and example data

	C3D and its documentation are available from https://github.com/mlupien/C3D.
	For this example, we will predict genomic interactions involving the promoters of abelson tyrosine-protein kinase 1 (ABL1) and the breakpoint cluster region protein (BCR) in K562 (chronic myelogenous leukemia) cells by calculating correlations across DNase-seq signals for 79 ENCODE cell lines. All of the data for this example is located in https://github.com/mlupien/C3D/tree/master/example.
	
    First, download the BED file called test_anchors.bed. This contains anchors for ABL1 and BCR.
    
        Note: this file was generated by entering the following command in the terminal:
        
            Rscript make_anchors_from_genes.R genes.txt 1000 0 test_anchors.bed
    
    Now, download signalMatrix.txt. This file is formatted as a matrix whose rows correspond to DNase I Hypersensitive Sites (DHSs) and columns correspond to cell lines. Each entry is the maximum DNase-seq signal over the given DHS in the given cell line.
    
        Note: signalMatrix.txt was generated by following these steps:
    
        1.  download a BEDGRAPH file for each of the 79 cell lines. This can be done by running download_79_ENCODE_cell_lines.sh from https://github.com/mlupien/C3D/tree/master/example.

            sh download_79_ENCODE_cell_lines.sh
	
	    This script requires wget. If it is not installed, instructions for downloading wget can be found from https://www.gnu.org/software/wget/ (http://gnuwin32.sourceforge.net/packages/wget.htm for windows).
        
            Aside: You can also download 126 BEDGRAPH files from ENCODE and 53 BEDGRAPH files from Roadmap for the db.
            
            sh download_126_encode_bedgraphs.sh
            sh download_convert_53_roadmap_bigwigs.sh
            
        One of these scripts requires gzip and UCSC command line bioinformatic utilities. If they are not installed, instructions for downloading gzip and UCSC command line bioinformatic utilities can be found from http://www.gzip.org/ (http://gnuwin32.sourceforge.net/packages/gzip.htm for windows) and https://github.com/ENCODE-DCC/kentUtils, respectively. 

	    2.  Create a file that lists all of the full paths and file names of the BEDGRAPH files. Each file name should be on a separate line. For example, the first three lines may be formatted as follows:
       
            /path/to/wgEncodeUwDnaseA549Aln_2Reps.norm5.rawsignal.bedgraph
            /path/to/wgEncodeUwDnaseAg04449Aln_2Reps.norm5.rawsignal.bedgraph
            /path/to/wgEncodeUwDnaseAg04450Aln_2Reps.norm5.rawsignal.bedgraph

        3.  download a BED file of DHSs in K562.
    
            wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwDnase/wgEncodeUwDnaseK562PkRep1.narrowPeak.gz
            gunzip wgEncodeUwDnaseK562PkRep1.narrowPeak.gz
        
        This command requires gzip.
        
        4.  Edit the config.txt located in https://github.com/mlupien/C3D/tree/master/example by replacing the line starting with 'matrix=' with:
            
            reference=</path/to/>wgEncodeUwDnaseK562PkRep1.narrowPeak
            db=</path/to/list_of_bedgraphs.txt>
            
        Set reference, db, anchor and outDirectory with appropriate paths and file names.
        
Configure and run C3D

    Download and edit config.txt by setting anchor, outDirectory and matrix with appropriate paths and file names. Run C3D with the following command:
    
        sh c3d config.txt

C3D will output the results to results_< timestamp >.txt. Arc diagrams for each anchor can be found in figures_< timestamp >.pdf. Tracks will be generated in files ending in .tracks.txt and .bedGraph. Files ending with .tracks.txt can be uploaded to the UCSC Genome Browser (https://genome.ucsc.edu/cgi-bin/hgGateway) to visualize the correlations. Alternatively, tracks ending in .bedGraph can be uploaded to IGV. Instructions for downloading and using IGV can be found at https://software.broadinstitute.org/software/igv/download. All of these files will be located in the specified outDirectory.
