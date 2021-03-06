# FRAGTE2
FRAGTE2: an enhanced algorithm to pre-select closely related genomes for bacterial species demarcation

# Please Cite
If you use FRAGTE in your publication, please cite: 

Zhou, Y., et al., A completeness-independent method for pre-selection of closely related genomes for species
 delineation in prokaryotes. BMC Genomics, 2020. 21(1): p. 183.
 
 If you use FRAGTE2 in your publication, please cite: 
 
Zeng, J., Wang, Y., Wu, Z. & Zhou, Y. FRAGTE2: An Enhanced Algorithm to Pre-Select Closely Related Genomes for Bacterial Species Demarcation. Frontiers in Microbiology 13, doi:10.3389/fmicb.2022.847439 (2022).

# Version
2.0

# Developer
Yizhuang Zhou (zhouyizhuang3@163.com)

# Affiliation
Guilin Medical University

# Support platform
FRAGTE2 was developed and maintained on Linux platform.

# Prerequisite
1. Perl5 with CPAN, FindBin, File::Basename, POSIX, Time::Local, Cwd, Getopt::Long and Statistics::R modules
(./autobuild_auxiliary will use CPAN to check and install the other 7 Perl modules)

# Installation
FRAGTE2 was developed in Perl language with embedded R programs. It just needs to install the obove 7 Perl modules. If these Perl modules have been installed, you directly use FRAGTE2. Otherwise, please type:
./autobuild_auxiliary

# Description
FRAGTE2 program is in bin directory; all tested data are in Data directory including assembly_accession, species_taxid and organism_name separatted by Tab; other scripts are in Scripts directory.
For FRAGTE1, please refer to https://github.com/Yizhuangzhou/FRAGTE. 

# Usage
1. Link genomes

perl Scripts/Linking_Genome.pl [GenomeInfo][outdir][output][len]

[GenomeInfo] has 7 fields separatted by Tab, including:
Assembly_accession,species_taxid,organism_name,strain_name,assembly_level,total genome size (in bp) and file for genome
  
[outdir] the directory to save files for linked genomes

[output] the output file containing information for each processed genome, including 7 fields separatted by Tab as follows:
Assembly_accession,species_taxid,organism_name,Average Size,assembly_level,total genome size (in bp) and file for genome

[len] the fragment length used for training for the Na??ve Bayesian Classifier, generally 10000
  
2. Run FRAGTE2

Examples including:

 perl bin/FRAGTE2.pl --Qfile  <Qfile> --Rfile <Rfile>
  
 perl bin/FRAGTE2.pl --Qfile  <Qfile> --Rfile <Rfile> --QPrefix <Prefix> --RPrefix <Prefix>
  
 perl bin/FRAGTE2.pl --Qfile  <Qfile> --Rfile <Rfile> --QPrefix <Prefix> --RPrefix <Prefix> --Outdir <outdir> 
 
 perl bin/FRAGTE2.pl --Afile  <Afile>
 
 perl bin/FRAGTE2.pl --Afile  <Afile> --APrefix <Prefix>
 
 perl bin/FRAGTE2.pl --Afile  <Afile> --APrefix <Prefix> --Outdir <outdir> 
## Note
(1)The input file for "--Qfile" or "--Rfile" has 7 fields separatted by Tab, including: 
  Assembly_accession,species_taxid,organism_name,Average Size,assembly_level,total genome size (in bp)
  and file for genome
  
(2) file for genome,the file with absolute path containing sequences in fasta format.
  
(3) "--Qfile" and "--Rfile" should be used together. Otherwise, "--Afile" is used instead which includes both 
  queries and references.
  
(4) the output file named "*_Pairs_byFRAGTE2.xls" for pairs sieved by FRAGTE2 includes 5 fields separatted by
  Tab as follows: Query ID, Reference ID, PCCD, GSCq and GSCr.
  
(5) the output file named "*_Pairs_byFRAGTE2.log" for pairs unsieve by FRAGTE2 also includes 5 fields separatted
  includes 5 fields separatted by Table as follows: Query ID, Reference ID, PCCD, GSCq and GSCr.
  
# Description for the direcory Data
  
Genomes_4Simulated.xls: genome accessions for 6230 complete genomes for generating simulated genomes
  
RealAssembly_Queries.xls: the 61,914 query accessions for real genomes
  
RealAssembly_References.xls: the 5680 reference accessions for real genomes
  
MAGs.xls: 3032 accessions for MAGs
  
IntraMAGs.xls: the intraspecific pairs for MAGs, which is used for assessing sieving performance on MAGs. 

# Support
If you need some other scripts and other materials or encounter bugs, problems, or have any suggestions, please contact me via zhouyizhuang3@163.com
