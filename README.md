# autoMLSA 

The goal of this project is to enable users to generate high quality and robust **M**ulti-**L**ocus **S**equence **A**lignment phylogenetic trees with minimal effort. The default settings for the software work in many situations, and the software as a whole works with minimal user intervention. 

An appropriate set of reference gene sequences (nucleotide or amino acid), the desired NCBI or local BLAST databases, an e-mail address, and a unique run ID are required as input to the software.

# Installation

No installation is required if dependencies are installed. See these optional steps to save some time running the program if you are going to be the sole person using the software.

For questions about your PATH, see [this link](https://superuser.com/questions/284342/what-are-path-and-other-environment-variables-and-how-can-i-set-or-use-them) for a decent basic explanation.

1. Edit line 42 to include your e-mail address. Alternatively, set your EMAIL environment variable to your e-mail address. Your e-mail address will only be used to query the NCBI databases using the E-utilites. This is to allow NCBI to contact you if there is a problem with your requests or if you are abusing the system.

   Environment variables:

   ```bash
   export EMAIL youremail@yourschool.edu
   ```

   or

   ```tcsh
   setenv EMAIL youremail@yourschool.edu
   ```

   Make the change:

    ```perl
    my $email = '';
    ```
    to
    ```perl
    my $email = 'youremail@yourschool.edu';
    ```

2. Edit line 56 to the path of your blast executables if they are not in your PATH.

3. If you wish to include autoMLSA.pl in your path, but don't want to move the ./scripts folder to the same location, put the full path to the scripts folder in line 60. For example:

    * Move autoMLSA.pl to a location in your PATH (e.g. /usr/bin or /usr/local/bin or ~/bin)
    * Add the full path to the unpacked scripts directory to line 60 (eg. ~/libs/autoMLSA/scripts)

4. Add the full path to the Gblocks and/or noisy to lines 58/59 if the programs are not in your PATH.

5. You can also change the optional defaults starting on line 74 if you are interested in changing them. For example:
Change your 'prog' from 'tblastn':

    ```perl
    my %defaults = (
                     'prog'         = 'tblastn',
                     'nr'           = 1,
                     'complete'     = 1,
                     'evalue'       = '1e-5',
                     'local_evalue' = '1e-5',
                     'target'       = 500,
                     'coverage'     = 50
                   );
    ```
    
    To 'blastn':
    
    ```perl
    my %defaults = (
                     'prog'         = 'blastn',
                     'nr'           = 1,
                     'complete'     = 1,
                     'evalue'       = '1e-5',
                     'local_evalue' = '1e-5',
                     'target'       = 500,
                     'coverage'     = 50
                   );
    ```
    
    If you are interested in making blastn the default rather than tblastn.

6. If you plan on using RAxML, add the path to your RAxML executables to line 40 of autoMLSA-raxml.pl if they are not in your PATH.

7. Add the path to RAxML in ./scripts/ProteinModelSelection.pl and, if you so choose, change the type of executable as well.

## Dependencies

Required 

* BLAST+ version 2.2.31+
* Perl version 5.10.1 or higher
* NCBI Entrez Direct scripts - Run edirect-dl.pl from scripts folder. [See here](https://www.ncbi.nlm.nih.gov/books/NBK179288/) for more information
* Perl Modules
  * BioPerl (Bio::Perl) v1.7.000 [(1.007000)](http://search.cpan.org/~cjfields/BioPerl-1.007000/)
  * ~~BioPerl Eutilities (Bio::DB::EUtilities)~~ **NO LONGER REQUIRED**

Multiple Sequence Alignment software
* [MAFFT v7](http://mafft.cbrc.jp/alignment/software/) **recommended**

Phylogenetic Tree Inference software
* [RAxML v8](http://sco.h-its.org/exelixis/web/software/raxml/index.html) **recommended**

Optional

* [Gblocks](http://molevol.cmima.csic.es/castresana/Gblocks.html) - No longer recommended. See [Tan et al. 2015](http://sysbio.oxfordjournals.org/content/64/5/778) as to why.
* [noisy](http://www.bioinf.uni-leipzig.de/Software/noisy/) **recommended** See [Tan et al. 2015](http://sysbio.oxfordjournals.org/content/64/5/778) for benchmarking of alignment filtering software.

# Pipeline Description

This pipeline uses a set of input marker sequences for a BLAST search, finds the corresponding gene or protein sequences in the target BLAST databases (nt and custom local databases), determines an appropriate organism name using information present in the NCBI taxonomy and nucleotide databases, determines which genomes have all genes of interest, aligns each gene individually, filters identical alignments with a log file to help determine which alignments are identical, determines the appropriate substitution model for each gene, and then concatenates the genes into an alignment while saving the partition data.

The output from the main script is a concatenated alignment file, the appropriate substitution model to provide to RAxML, and a partitioned dataset file that is formatted properly for input into RAxML. The RAxML wrapper provided along with the main script (autoMLSA-raxml.pl) will allow you to generate a phylogenetic tree based on the MLSA dataset. Alternatively, you can use these files to do your own custom analysis using RAxML or other phylogenetic tree generation software.

# Usage

This software is designed to be as straightforward to use as possible, with as minimal input as possible. Once the dependencies are installed, running autoMLSA.pl --help should provide a lengthy help message to remind you of all of the possible options with this software.

Minimally, two pieces of data are required to complete a full run:

1. A set of housekeeping or marker genes sufficient for MLSA phylogenetics. Ideally, a set of marker genes for your species of interest is already known. If not, the easiest option for determining a set of genes is to examine the literature of a related species, genus, or family, and use the same set of genes for your organism.

2. The NCBI taxonomy term for your clade of interest. Use the [Taxonomy](https://www.ncbi.nlm.nih.gov/guide/taxonomy/) database at NCBI to find the proper term if necessary.

You will need all of the gene sequences (nt or aa) from a reference organism, and the sequences can be contained in a single multi-FASTA file, or multiple single FASTA files. In general, best practice would include renaming the FASTA header to the gene of interest. For example:

`>gi|84028820|sp|P0AES6.2|GYRB_ECOLI RecName: Full=DNA gyrase subunit B; AltName: Full=Type IIA topoisomerase subunit GyrB`

Should be:

`>gyrB gi|84028820|sp|P0AES6.2|GYRB_ECOLI RecName: Full=DNA gyrase subunit B; AltName: Full=Type IIA topoisomerase subunit GyrB`

This is not completely necessary, but will help with naming files in a meaningful way throughout the process. The FASTA header ID (everything before the first space) is used to name the output files pertaining to each marker gene used.

If you are only interested in sequences found in NCBI (nt) then this will be all you need. You will run the program like this:

`autoMLSA.pl -email yourname@yourinstitution.edu -runid sample_run_id -threads 4 -entrez_query "YourTaxonomyTerm[ORGANISM]" -- markers.fasta`

Output will be sample\_run\_id.concat.derep and sample\_run\_id.partition, and the substitution model to supply to RAxML is printed to the screen and present in the log file at the end of the run.

By default, a protein query is expected, and a TBLASTN analysis is done. This helps eliminate any annotation discrepencies with genomes present in the NCBI nt database. Alignment software can more accurately align protein sequences and accurate alignments are preferred in most cases, with the cost being less phylogenetic resolution. Empirically, the default cutoffs for coverage and evalue are reasonable. Depending on the size of your clade of interest, you will have to increase the target cutoff to a larger value, especially for clades with several well represented species.

## Generating a phylogenetic tree using autoMLSA-raxml.pl

The files required for autoMLSA-raxml.pl will be printed at the end of the autoMLSA.pl run. You will use them like so:

`autoMLSA-raxml.pl -runid sample_run_id -partition sample_run_id/sample_run_id.partition -model PROTGAMMALG -alignment sample_run_id/sample_run_id.concat.derep -SSE3 -PTHREADS 8 -initial`

After that completes, you would run the mlsearch and bootstrap programs:

`autoMLSA-raxml.pl -runid sample_run_id -partition sample_run_id/sample_run_id.partition -model PROTGAMMALG -alignment sample_run_id/sample_run_id.concat.derep -SSE3 -PTHREADS 16 -mlsearch`

`autoMLSA-raxml.pl -runid sample_run_id -partition sample_run_id/sample_run_id.partition -model PROTGAMMALG -alignment sample_run_id/sample_run_id.concat.derep -SSE3 -PTHREADS 8 -boostrap`

By default, 100 ML searches are performed, and the bootstrap analysis completes based upon the autoMRE bootstopping criterion.

After both the ML search and bootstrap analysis complete, you can generate the final trees with this command:

`autoMLSA-raxml.pl -runid sample_run_id -partition sample_run_id/sample_run_id.partition -model PROTGAMMALG -alignment sample_run_id/sample_run_id.concat.derep -SSE3 -final`

The trees will take the form of:

`RAxML_bipartitions.sample_run_id.final.nwk` 
and 
`RAxML_bipartitionsBranchLabels.sample_run_id.final.nwk`

If you are interested in re-rooting the tree, you would provide the full strain name with: `-outgroup "Escherichia_coli_K12"` or a group of strains with `-outgroup "Escherichia_coli_K12,Escherichia_coli_DSM_30083"`

The .nwk files can be viewed with any phylogenetic tree software, and we recommend using [iTOL online](http://itol.embl.de/).

## Using a local BLAST database

NCBI has removed the ability to search the WGS database using a remote BLAST search from the command-line. Therefore, if you require a search of the WGS database, you will need to download the genomes from the WGS database and include then as a local BLAST database. Use the download\_genomes.pl script to do so. You will need to update the path (or include the programs in your PATH) to the edirect command-line utilties, which should be installed already prior to running the autoMLSA.pl script. To download the genomes, for example, the Pseudomonas genomes, do `autoMLSA/scripts/download_genomes.pl -type mlsa -download Pseudomonas`. The genomes will be donwloaded to a file called Pseudomonas.fasta. You can then combine this file with your own genome sequence files, or generate a BLAST DB using makeblastdb, to include these sequences in your local MLSA run. You will want to turn off searching of the nr/nt database if you use this approach, as you will be downloading the genomes present in nt as well. Do this by including -nonr as a flag for your autoMLSA run.

When using your own BLAST database, the FASTA files have to be formatted properly prior to BLAST DB generation, excluding any FASTA files that have been downloaded from NCBI (if the FASTA files have a recognizable accession number, they will be detected). FASTA files containing the genome sequence need to be formatted properly to be recognized and understood by autoMLSA. This process can be accomplished with the fasta\_format.pl script located in the autoMLSA/scripts folder. The formatting script expects each genome to be contained in an individual FASTA file, with the format of strainName.fasta. 

For example, if you were using the Escherichia coli strain K12, the genome file would be named K12.fasta. You would then provide the species name with the -genome flag:

`autoMLSA/scripts/fasta_format.pl -genome 'Escherichia coli' -overwrite K12.fasta`

If you had several genomes with the appropriate file names in a particular directory, you would do this:

`autoMLSA/scripts/fasta_format.pl -genome 'Escherichia coli' -overwrite *.fasta`

You must be careful with this rename script, as it will overwrite each file in turn, so preferably you would make a copy of each genome to convert to autoMLSA format.

The autoMLSA format is straightforward. For example:

`>gnl|K12|contig1 [Escherichia coli K12]`

The strain identifier is in the second field, and is used to link all of the sequences from this organism together. The full strain name is in square braces in the FASTA description field.

After all of the genome sequences are formatted properly using the fasta\_format.pl script (or manually), then you must concatenate all of the genome sequences into one file. `cat *.fasta * all.fas` works well on \*nix based systems and under cygwin. Then, using the makeblastdb tool provided with the NCBI BLAST suite, you would create the local bast DB as below:

`makeblastdb -dbtype prot -in all.fas`

You would provide the full path to this database like so:

`autoMLSA.pl -email yourname@yourinstitution.edu -runid sample_run_id -threads 4 -entrez_query "YourTaxonomyTerm[ORGANISM]" -local_db /path/to/all.fas -- markers.fasta`

If you have several local BLAST DBs, you must name them all uniquely, otherwise the results will be overwritten, and sequences from only one of the BLAST DBs will be contained in the final output.

## Understanding the intermediate files

There are several intermediate files produced by the script, some of which are useful to understand what type of processing is going on during the run. Each input sequence has several output files, depending on the number of BLAST databases used as targets. BLAST output to the nt database takes the form of FASTA\_ID.out and local databases have the form of FASTA\_ID\_vs\_localdbname.local.out. Each BLAST output also has a corresponding key (.key) file that provides the mapping of IDs to assemblies and genome names, as well as the FASTA formatted output (.fas). After the FASTA formatted sequences are extracted from the output, sequences from all target databases are combined into a single file (.all.fas) and sequences from genomes that are missing one or more of the gene sequences are removed from the analysis (.all.fas.sorted). 

Then, each sorted FASTA file is aligned using MAFFT L-INS-i by default (.all.aln), and trimmed with Gblocks if desired (all.aln-gb). The sorted/trimmed alignments are then concatenated (sample\_run\_id.concat) and dereplicated to remove identical concatenated alignments (.concat.derep, .concat.derep.log keeps a log of which sequences are dereplicated). If a certain sequence is dereplicated that you are interested in keeping, you can extract the line(s) into a text file to rereplicate the identical sequences. You give this information to the script with the flag -rerep /path/to/dereplicated/sequences.txt. The replicated concatenated alignment will have the name sample\_run\_id.concat.derep.rerep. The best substitution model is then determined for each partition with the intermediate RAxML files resulting in the model\_test folder.

## Checking on your progress

Detailed progress messages are printed both to STDERR and to a log file (present in the ./log folder). Run failures due to scripts or external programs failing are somtimes handled with an appropriate error message in the log file, but some errors are still unhandled and the error messages will only print to STDERR. Generally, unhandled errors include unexpected characters in the alignment (e.g. something other than ACGT in nucleotide alignments, X in protein alignments), problems with write access, unusual characters in ID names, and problems with RAxML runs. Problems with RAxML must be handled based on information in the RAxML\_info files found in the model\_test folder and the trees folders.

Two steps where manual restarts would commonly be required are when rereplication is desired, or when a proper genome name cannot be determined automatically from NCBI (files named .manual). Otherwise, the program will run to completion without any additional input after the command is entered.

# Detailed explanation of options

**-help|h**

Print a brief help message and exits.

**-man**

Print a verbose help message and exits.

**-quiet**

Turns off progress messages. Use noquiet to turn messages back on.

**-email REQUIRED**

Enter your email.  NCBI requires this information to continue. Can also set EMAIL environment variable.

**-runid REQUIRED**

Name of the run. Allows for continuation of runs and keeps a log of the progress.

**-config** (/path/to/config.txt)

Allows for input of a config file.  Command line parameters will supercede any parameters in the config file.  Useful for re-running an analysis with different parameters.  Will be generated prior to start of remote BLAST search.

**-log|nolog** [logging on]

Using -nolog turns off logging. Logfile is in the format of runid.log, where runid is provided at command line. Placed in the ./logs/ dir.

**-prog** [tblastn]

Sets blast program to run.  blastn and tblastn are valid options.

**-complete|nocomplete** [on]

Option searches for complete records only.  Adds "AND complete genome" to searches that have an entrez query (supplied by -entrez\_query 'text').

**-nr|nonr** [on]

Searches the nr/nt database.

**-local\_db** (/path/to/blastdb1 /path/to/blastdb2)

Allows a search of one or more local blastdbs in addition to the remote databases.

**-threads** [1]

Can run the local blast and mafft alignments in multithreaded manner, if desired. (Can specify threaded options for other alignment programs using -align\_params)

**-align\_prog** [mafft-linsi]

The linsi algorithm of mafft is listed by default.  Any alignment program that prints to STDOUT is supported natively.  See body of script to insert specific code for other algorithms.

**-align\_params**

Supply a string in quotes that will be added to the end of the alignment program for program specific commands. (e.g. '--legacygappenalty' for mafft)

**-trimmer** (gblocks OR noisy)

Run Gblocks or noisy to trim out gappy parts of the alignment.

**-trimmer_params**

Can flag for trimmer to run with additional parameters using -trimmer\_params='-b3=10 -b4=5 -b5=h' syntax at the commandline. noisy will automatically set seqtype parameter.

**-rereplicate** (filename)

Supply a filename with tab delimited names to rereplicate prior to model selection and tree generation.  (Set the flag with no filename on first run to stop the run at this step.  Then pull out the genomes you want to replicate from the dereplicate log.

**-cleanup**

Cleans up files and allows for addition of more genes to the analysis. This will delete previous results, so start a new analysis with a new runid if you are interested in saving your results.  Automatically cleans up old files when new genes are added to the analysis.  Previous alignments are not usable as some genomes may not contain all genes once more genes are added to the analysis.

**-relaxed**

Decreases match and mismatch scores, as well as gap costs to allow for more distant matches for nucleotide searches.

##PARAMETERS

Except for evalue, by default the local values are set to equal the 'nr/nt' parameters (e.g. -local\_target = -target).

**-evalue, -local\_evalue**

Sets the e-value cutoffs.

**-target, -local\_target**

Sets the limit on number of alignments returned. 

**-coverage and -local\_cov**

Sets the threshold for coverage to include.  By default, includes all hits, regardless of coverage.  Values between 25 and 100 are valid [ex: 50 for alignment of 50% of the query and subject sequences]. 

**-entrez\_query**

Sets an entrez query limiter to blast search.  Useful for searching only within a species, genus, etc.

Example: For looking only at Actinobacteria, use 'Actinobacteria[ORGANISM]'

##ADVANCED

**-debug\_cleanup**

Removes all files generated after BLAST search output.  Allows for restarting of jobs that may have failed and cannot continue due to unusual circumstances.

**-checkpoint** (fasta OR model)

Stops script prior to filtering for genes found in all genomes or prior to alignment.  Useful for acquiring gene sequences without requiring completion of the entire pipeline, or doing the BLAST searches using a slower machine, saving the alignments for a more powerful machine.

# Description of Included Files

| Filename | Description |
| -------- | ----------- |
| autoMLSA.pl | Main program |
| README.md | This README file |
| LICENSE | Description of GPL License |
| scripts/ProteinModelSelection.pl | Finds best protein model for each protein alignment |
| scripts/autoMLSA-concat.pl | Concatenates alignments; renames accessions to species names |
| scripts/autoMLSA-derep.pl | Removes identical aligned sequences |
| scripts/autoMLSA-fasta\_rename.pl | Use to rename non-concatenated sequences **NOT WORKING** |
| scripts/autoMLSA-filter.pl | Removes genomes that are missing one or more input genes/proteins |
| scripts/autoMLSA-raxml.pl | Use this to generate the phylogenetic tree after concatenation |
| scripts/autoMLSA-rerep.pl | Use this to add essential sequences back into the alignment so they appear in the final tree |
| scripts/autoMLSA-searchio.pl | This script is used to parse tab-delimited BLAST output for sequence and accession data |
| scripts/auto\_edirect.pl | This script runs the NCBI edirect command-line tools to find species, strains, and assembly information |
| scripts/edirect-dl.pl | Use this to download the edirect tools if you do not have them installed already |
| scripts/fasta\_format.pl | Use this to format user-generated FASTA data |
| scripts/download\_genomes.pl | Use this to download a genome dataset from NCBI to use in BLAST search. |

# Citations

Please cite:

```
Davis II EW, Weisberg AJ, Tabima JF, Grunwald NJ, Chang JH. (2016) Gall-ID: tools for genotyping gall-causing phytopathogenic bacteria. PeerJ 4:e2222. https://doi.org/10.7717/peerj.2222
```

if you use this software. Also, please cite these other publications as integral parts of autoMLSA:
```
Stamatakis A. (2014) RAxML version 8: a tool for phylogenetic analysis and post-analysis of large phylogenies. Bioinformatics 30(9):1312–1313. https://doi.org/10.1093/bioinformatics/btu033

Katoh K, Standley DM. (2013) MAFFT Multiple Sequence Alignment Software Version 7: Improvements in Performance and Usability. Mol. Biol. Evol. 30(4):772–780. https://doi.org/10.1093/molbev/mst010

Dress AWM, Flamm C, Fritzsch G, Grünewald S, Kruspe M, Prohaska SJ, & Stadler PF. (2008). Noisy: identification of problematic columns in multiple sequence alignments. Algorithms for Molecular Biology : AMB, 3, 7. http://doi.org/10.1186/1748-7188-3-7
```

# History

v2.1.0 - 2016-12-12 - Added -relaxed command to allow more dis-similar BLAST results when using nucleotide queries. Added support for failed elink attempts.

v2.0.0 - 2016-12-06 - New version dependent on NCBI edirect. No longer requires BioPerl EUtilities. WGS no longer supported. INCOMPATIBLE with previous versions.

v1.0.1 - 2016-10-10 - Dev branch started.

v1.0 - 2016-04-20 - First revision released to GitHub.

# Credits

* [Edward Davis](mailto:davised.dev@gmail.com)
* [Dr. Jeff Chang](mailto:changj@science.oregonstate.edu)

# License

Please see the LICENSE file for more information.
