#!/usr/bin/env perl
#######################################################################
#
# COPYRIGHT NOTICE
#
# autoMLSA - A program to automate Multi-Locus Sequence Analysis
# Phylogenetic Trees using NCBI databases and local datasets
#
# Copyright (C) 2015
# Edward Davis
# Jeff Chang
#
# This file is a part of autoMLSA.
#
# autoMLSA is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# autoMLSA is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more detail.
#
# You should have received a copy of the GNU General Public License
# along with autoMLSA.  If not, see <http://www.gnu.org/licenses/>.
#
#######################################################################

use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use Bio::SeqIO;

use File::Spec;
use Cwd 'abs_path';

my $version = '2.1.0';
my $date = 'December 12, 2016';
my $email = '';    #Can set to default
if (defined $ENV{'EMAIL'}) {
    $email ||= $ENV{'EMAIL'};
}

#Setting up paths properly

my $incommand = join( " ", $0, @ARGV );
my ( $svol, $sdir, $sfile ) = File::Spec->splitpath($0);
$sdir .= "scripts/";

#Global Path variables

#Change to /path/to/blastdir if blast executables are not in PATH
my $localblastdir = '';

my $gblockspath = 'Gblocks';    #Change to /path/to/Gblocks if not in PATH
my $noisypath = 'noisy';        #Change to /path/to/noisy if not in PATH
my $scriptdir   = ''
  ; #can change this to '/path/to/scriptdir/' if all scripts below are in the same directory outside of ./scripts
my $proteinmodelpath =
  ( $scriptdir ? $scriptdir : $sdir ) . 'ProteinModelSelection.pl';
my $searchiopath = ( $scriptdir ? $scriptdir : $sdir ) . 'autoMLSA-searchio.pl';
my $elinkpath    = ( $scriptdir ? $scriptdir : $sdir ) . 'auto_edirect.pl';
my $dereplicatepath = ( $scriptdir ? $scriptdir : $sdir ) . 'autoMLSA-derep.pl';
my $rereplicatepath = ( $scriptdir ? $scriptdir : $sdir ) . 'autoMLSA-rerep.pl';
my $mlsaconcatpath = ( $scriptdir ? $scriptdir : $sdir ) . 'autoMLSA-concat.pl';
my $filtergenomespath =
  ( $scriptdir ? $scriptdir : $sdir ) . 'autoMLSA-filter.pl';

#Set optional defaults

my %defaults = (
                 'prog'         => 'tblastn',
                 'remote'       => 'nt',
                 'complete'     => 1,
                 'evalue'       => '1e-5',
                 'local_evalue' => '1e-5',
                 'target'       => 500,
                 'coverage'     => 50
               );

#Other global variables

my $log = 1;
my $time;
my $runid;
my $config;
my $newconfig;
my %files;
my %inputs;
my @inputs;
my $blast_check;
my $runpath;
my $logpath;
my $logfile;

#Set some defaults

my %options = (
                'help'   => 0,
                'man'    => 0,
                'quiet'  => 0,
                'email'  => \$email,
                'runid'  => \$runid,
                'config' => \$config,
                'log'    => \$log,
                'version' => 0
              );


#Get command-line options

my $signal = GetOptions(
                         \%options,            'help|h',
                         'man',                'runid|r=s',
                         'email=s',            'config=s',
                         'quiet!',             'prog=s',
                         'evalue=f',           'target=i',
                         'coverage=f',         'local_evalue=f',
                         'local_target=i',     'local_cov=f',
                         'entrez_query=s',     'remote=s',
                         'complete!',
                         'local_db=s@',        'threads=i',
                         'align_prog=s',       'align_params=s',
                         'trimmer=s',            'rereplicate:s',
                         'cleanup',            'log!',
                         'clear_input',        'clear_dbs',
                         'debug_cleanup',      'concat',
                         'checkpoint=s',       'trimmer_params=s',
                         'version|v',          'relaxed'
                       );

#Print help statements if specified
pod2usage( -verbose => 1 ) if $options{help} == 1;
pod2usage( -verbose => 2 ) if $options{man} == 1;

#Print version statements if specified

if ( $options{version} ) {
    print STDOUT "$0 version $version.\nBuilt on $date.\n";
    exit(0);
}

#Die if GetOptions returns false [unknown option passed]
if ( !$signal ) {
    pod2usage(
        -verbose => 0,
        -message =>
          "Unknown option entered.  Run command with flag -h for available options.\n",
        -exitval => 3
    );
}

if ( !( defined $runid ) ) {
    pod2usage(
        -verbose => 0,
        -message =>
          "No runid given.  This is required to continue.  Please specify a runid using -runid and resubmit. Run the command with option -h for more information.\n",
        -exitval => 4
    );
} else {
    if ( -d "../$runid" ) {
        if ( -d $runid ) {
            $runpath = File::Spec->rel2abs($runid);
            $logpath = File::Spec->rel2abs("./log/");
        } else {
            $runpath = File::Spec->rel2abs('.');
            $logpath = File::Spec->rel2abs("../log/");
        }
    } else {
        $runpath = File::Spec->rel2abs($runid);
        $logpath = File::Spec->rel2abs("./log/");
    }
    $logfile = $logpath . "/$runid.log";
    logger("Command as submitted:\n");
    logger("$incommand\n\n");
    if ( !-d $runpath ) {
        logger("Generating folder $runid\n");
        system("mkdir $runid") == 0
          or die "Unable to make dir $runid. Check folder permissions.\n";
    } else {
        logger(
            "Dir found for run $runid.  Attempting to load previous results and continue.\n"
        );
        if ( -e "$runpath/$runid.config" ) {
            logger("Config file $runpath/$runid.config found.\n");
            $config = "$runpath/$runid.config";
        } else {
            logger("No config file found for $runid.\n");
        }
    }
    if ($config) {
        if ( -e "$config" ) {
            logger(
                "Loading values from $config. Command-line parameters will overwrite these.\n"
            );
            open CONFIG, "$config"
              or die
              "Unable to open config file.  Check permissions and runid and try again.\n";
            while (<CONFIG>) {
                my $line = $_;
                chomp($line);
                next if $line =~ /^#/;
                my ( $key, $value ) = split( '=', $line, 2 );
                my @junk;
                ( $value, @junk ) = split( '#', $value );
                if ( $key eq 'runid' ) {
                    next;
                } elsif ( $key eq 'inputs' ) {
                    if ( !$options{clear_input} ) {
                        my @inputs = split( ",", $value );
                        foreach my $input (@inputs) {
                            logger(
                                "Input file $input already found. Removing duplicate value.\n"
                              )
                              if ( defined $inputs{$input} );
                            $inputs{$input} = 1;
                        }
                    }
                } elsif ( $key eq 'email' ) {
                    $email = $value;
                } elsif ( $key eq 'local_db' ) {
                    if ( !$options{clear_dbs} ) {
                        my @local_dbs = split( ",", $value );
                        push( @{ $options{local_db} }, @local_dbs );
                    }
                } else {
                    $options{$key} //= $value;
                }
            }
            close CONFIG;
        }
    } else {
        system("touch $runpath/$runid.config") == 0
          or die "Unable to generate config file. Check permissions.\n";
    }
    $newconfig = "$runpath/$runid.config";
}

#Set some defaults
$options{complete} //= $defaults{complete};
#$options{remote}   //= $defaults{remote};

if ( scalar(@ARGV) == 0 && scalar( keys %inputs ) == 0 ) {
    pod2usage(
              -verbose => 1,
              -message =>
                "No inputs given.  Check options and inputs and try again.\n\n",
              -exitval => 5
             );
} else {
    if ( scalar(@ARGV) > 0 ) {
        $options{cleanup} = 1;
    }
    foreach my $input (@ARGV) {
        $input = abs_path($input);
        logger("Input file $input already found. Removing duplicate value.\n")
          if ( defined $inputs{$input} );
        $inputs{$input} = 1;
    }
}

if ( !$email ) {
    logger("Email is required. Please supply using -email flag.");
    die("\n");
}

@inputs = keys %inputs;

foreach my $infile (@inputs) {
    unless ( -e $infile ) {
        pod2usage(
            -verbose => 0,
            -msg =>
              "$infile does not exist. Check your parameters and the config file for paths and try again\n",
            -exitval => 6
        );
    }
}

$options{prog} //= $defaults{prog};

if ( defined $options{prog} ) {
    if ( $options{prog} =~ /^blastn$/ || $options{prog} =~ /^tblastn$/ ) {
        $options{dbfrom} = 'nucleotide';
    } else {
        pod2usage(
            -verbose => 0,
            -msg =>
              "Correct blast type is blastn or tblastn\nUse option -h for more information\n",
            -exitval => 7
        );
    }
}

if ( defined $options{threads} ) {
    my $check = `lscpu | grep -G '^CPU(s):' | grep -o -E '[0-9]+' `;
    if ( $? != 0 ) {
        logger(
            "Unable to check threads setting. Ensure you have enetered the proper number of threads or performance will be degraded (set to less than # of processors).\n"
        );
    } elsif ( $options{threads} > $check ) {
        logger(
            "Asking for more threads than cores are available.  Please check your parameters and try again."
        );
        die("\n");
    } else {
        logger(
            "Check to ensure you mean to use many threads ($options{threads}).\n"
          )
          if $options{threads} > 9;
    }
} else {
    $options{threads} = 1;
}

if ( !defined $options{align_prog} ) {
    logger("No alignment program chosen.  Setting to use mafft-linsi.\n");
    $options{align_prog} = 'mafft-linsi';
}

logger("Assuming blast executables in PATH\n") if $localblastdir eq '';
if ( $options{trimmer} ) {
    if ( $options{trimmer} =~ /[nN]oisy|[gG]blocks/ ) {
        if ( $options{trimmer} =~ /[gG]blocks/ ) {
            $options{trimmer} = 'Gblocks';
            if ( ! $options{trimmer_params} ) {
                $options{trimmer_params} = '-s=y -p=y';
            }
            logger("Assuming Gblocks is in PATH\n") if $gblockspath eq 'Gblocks';
            `which $gblockspath`;
            if ( $? != 0 ) {
                logger(
                    "There was a problem finding Gblocks.  Check your PATH to ensure Gblocks is present or provide the complete path to blast in the script."
                );
                die("\n");
            }
        } elsif ( $options{trimmer} =~ /[nN]oisy/ ) {
            $options{trimmer} = 'noisy';
            if ( $options{prog} =~ /^blastn$/ ) {
                $options{trimmer_params} = '--seqtype N';
            } elsif ( $options{prog} =~ /^tblastn$/ ) {
                $options{trimmer_params} = '--seqtype P';
            }
            $options{trimmer_params} .= ' -s';

        }
    } else {
        logger("Options for --trimmer are noisy or gblocks. Check your settings and try again.\n");
        die("\n");
    }
}

if ( !defined($options{remote}) && !defined( $options{local_db} ) )
{
    logger(
        "No BLAST searches scheduled to run. Please select at least one (-remote nt or supply a local BLAST db to -local_db.)"
    );
    die("\n");
}

#Set defaults
$options{evalue}       //= $defaults{evalue};
$options{local_evalue} //= $defaults{local_evalue};
$options{target}       //= $defaults{target};
$options{coverage}     //= $defaults{coverage};
$options{local_target} //= $options{target};
$options{local_cov}    //= $options{coverage};

#check for blast executable
$blast_check = (
                   $localblastdir ne ''
                 ? $localblastdir . '/' . $options{prog}
                 : $options{prog}
               );
my @blast_version = `$blast_check -version`;

if ( grep( /2.2.31|2.[1-9]?[3-9].[\d]+/, @blast_version ) ) {
    logger("Found BLAST version 2.2.31 or greater\n");
} else {
    logger(
        "BLAST version 2.2.31 or greater is REQUIRED!  Make sure the full path is provided or include it in your PATH!\n"
    );
    die("\n");
}

if ( $? != 0 ) {
    logger(
        "There was a problem running $blast_check.  Check your PATH to ensure the blast program is included, or provide the complete path to blast using /path/to/blastdir"
    );
    die("\n");
}

if ( $options{complete} == 1 ) {
    if ( defined $options{entrez_query} ) {
        if ( $options{entrez_query} !~ /complete genome/ ) {
            $options{entrez_query} .= ' AND complete genome';
        }
    } else {
        $options{entrez_query} = 'complete genome';
    }
}

#Parameter validation
if (%options) {
    my $evalue   = 0;
    my $target   = 0;
    my $coverage = 0;
    my $death    = 0;
    foreach my $value ( $options{evalue},
                        $options{local_evalue} )
    {
        if ( $value >= 10.1 ) {
            $evalue++;
            $death++;
        }
    }
    foreach my $value ( $options{target},
                        $options{local_target} )
    {
        if ( $value > 20000 || $value < 1 ) {
            $target++;
            $death++;
        }
    }
    foreach
      my $value ( $options{coverage}, $options{local_cov} )
    {
        if ( $value > 100 || $value < 25 ) {
            $coverage++;
            $death++;
        }
    }
    print STDERR
      "Valid evalues are numbers 10 or less. Check your parameters and try again.\n"
      if $evalue > 0;
    print STDERR
      "Valid alignment values are integers between 1 and 20000. Check your parameters and try again.\n"
      if $target > 0;
    print STDERR
      "Valid coverage values are numbers between 100 and 25. Check your parameters and try again.\n"
      if $coverage > 0;
    die("\n") if $death > 0;
}

my @db_names = ( ".nhr", ".nin", ".nsq" );

if ( $options{local_db} ) {
    foreach my $db ( @{ $options{local_db} } ) {
        my $dbpath = abs_path($db);
        foreach my $file (@db_names) {
            my $testpath = $dbpath.$file;
            if (! -e $testpath ) {
                print STDERR "Unable to find file $testpath. Entered parameter: $db.\nCheck to ensure you have generated a blast DB for this file.\n";
                die("\n");
            }
        }
        $db = $dbpath;
    }
    my %db_check = map { $_ => 1 } @{ $options{local_db} };
    @{ $options{local_db} } = sort keys %db_check;
}

if ( $options{rereplicate} ) {
    if ( -e "$options{rereplicate}" ) {
        $options{rereplicate} = abs_path( $options{rereplicate} );
    }
}

#export config file
open CONFIG, ">$newconfig" or die "Unable to open config file $newconfig\n";
$time = localtime();
print CONFIG "#Config file generated for run $runid at $time\n";
print CONFIG join( "=", 'runid', $runid ) . "\n";
print CONFIG join( "=", 'inputs', join( ",", @inputs ) ) . "\n";
if ( $options{local_db} ) {
    print CONFIG join( "=", 'local_db', join( ",", @{ $options{local_db} } ) )
      . "\n";
}
foreach my $key ( sort keys %options ) {
    my @skip =
      qw(config log runid help man db dbfrom cleanup local_db clear_dbs clear_input debug_cleanup checkpoint version);
    next if ( grep { $_ eq $key } @skip );
    if ( defined( $options{$key} ) ) {
        if ( $key =~ /email/ ) {
            print CONFIG join( "=", $key, $email ) . "\n";
        } else {
            print CONFIG join( "=", $key, $options{$key} ) . "\n";
        }
    }
}
close CONFIG;

logger("---------------------------------------------------\n");
logger("---------------------------------------------------\n\n");

$time = localtime();

logger("Blast searches started at $time\n\n");

foreach my $infile (@inputs) {
    for ( my $j = 0 ; $j < 2 ; $j++ ) {
        if ( !defined($options{remote}) ) {
            $j = 1;    #Skip remote search
        }
        if ( $j == 1 && !defined( $options{local_db} ) ) {
            $j = 2;
            next;      #skip local blast if no local_db found
        }

        #Set values for each blast type
        my ( $evalue, $target, $entrez_query, $threads );
        my ( @dbs, @outfiles );
        if ( $j == 0 ) {
            push( @dbs, $options{remote} );
            $evalue       = $options{evalue};
            $target       = $options{target};
            $entrez_query = $options{entrez_query};
        } elsif ( $j == 1 ) {
            #            $db = $options{local_db};
            @dbs     = @{ $options{local_db} };
            $evalue  = $options{local_evalue};
            $target  = $options{local_target};
            $threads = $options{threads};
        }
        #Old outfmt version
#        my $outfmt =
#          '\'7 qseqid sseqid sacc pident qlen length evalue qcovhsp staxids sscinames stitle sseq\'';
        my $outfmt =
          '\'7 qseqid sseqid saccver pident qlen length evalue qcovhsp stitle sseq\'';

        #Open file to get sequences for BLAST
        my $i = 0;
        my $str = Bio::SeqIO->new( -file   => "$infile",
                                   -format => 'fasta' );
        my @input_seqs;
        while ( my $input = $str->next_seq() ) {
            $input_seqs[$i] = $input->id;
            $input_seqs[$i] =~ s/[\\|:*?"<>]//g;

            if ( $options{prog} =~ /tblastn/ && $input->alphabet eq 'dna' ) {
                logger(
                    "DNA query when protein expected... translating sequence.\n"
                );
                $input = $input->translate;
            }
            if (    $options{prog} =~ /^blastn/
                 && $input->alphabet eq 'protein' )
            {
                logger("Protein query when DNA expected... exiting");
                die("\n");
            }
            foreach my $db (@dbs) {
                my $outfile = "$runpath/$input_seqs[$i]";
                if ( $j == 0 ) {
                    $outfile .= '.out';
                } elsif ( $j == 1 ) {
                    my ( $junk1, $junk2, $dbname ) = File::Spec->splitpath($db);
                    $outfile .= '_vs_' . $dbname . '.local.out';
                }

                push( @{ $files{ $input_seqs[$i] }{'out'} }, $outfile );
                if ( -s $outfile ) {
                    logger(
                        "Blast output for $input_seqs[$i] to $db already found. Skipping...\n"
                    );
                    next;
                }
                $options{cleanup} = 1;
                gene_cleanup($outfile);
                logger(   "Performing BLAST search for "
                        . $input->id . " to "
                        . $db
                        . " from file "
                        . $infile
                        . "..." );
                logger("\nUsing these parameters:\n");
                logger( "Blast program     => " . $options{prog} . "\n" );
                logger( "e-value cutoff    => " . $evalue . "\n" );
                logger( "# of alignments   => " . $target . "\n" );
                if ( $j == 0 ) {
                    logger("Entrez query      => $entrez_query\n");
                }
                if ( $j == 1 ) {
                    logger("# of threads      => $threads\n");
                }
                my $query = "$runpath/tmp";
                my $tmp_out =
                  Bio::SeqIO->new( -file   => ">$query",
                                   -format => 'fasta' );

                $tmp_out->write_seq($input);

                if ( $^O eq 'cygwin' ) {
                    $query = File::Spec->abs2rel($query);
                    $outfile = File::Spec->abs2rel($outfile);
                }

                my $command = join( " ",
                                    $blast_check, "-out",
                                    $outfile,     "-evalue",
                                    $evalue,      "-db",
                                    $db,          "-max_target_seqs",
                                    $target,      "-query",
                                    $query,       "-outfmt",
                                    $outfmt );
                if ( $blast_check !~ /tblastn/ ) {
                    $command .= " " . join( " ", "-task", "blastn" );
                    if ($options{'relaxed'}){
                        $command .= " " . join( " ", '-gapopen', '1',
                                                     '-gapextend', '1',
                                                     '-reward', '1',
                                                     '-penalty', '-2');
                    }
                }
                if ( $j == 0 ) {
                    $command .= " "
                      . join( " ",
                              "-entrez_query", "\'$entrez_query\'" );

                }
                if ( $j == 0 ) {
                    $command .= " -remote";
                }
                if ( $j == 1 ) {
                    $command .= " " . join( " ", "-num_threads", "$threads" );
                }

                `$command`;

                if ( $? != 0 ) {
                    `rm -f $runpath/tmp`;
                    `rm -f $outfile` if ( -e $outfile );
                    logger(
                        "BLAST command failed. Exiting script; resubmit and try again"
                    );
                    die("\n");
                } else {
                    `rm -f $query`;
                    logger("Unable to remove tmp file.\n") if $? != 0;
                    filecheck( "BLAST results", $outfile );
                }
            }
        }
    }
}

$time = localtime();
logger("\nFinished blast searches at $time.\n\n");
logger("---------------------------------------------------\n");
logger("---------------------------------------------------\n\n");

logger("Extracting sequences from BLAST searches\n");

if ( $options{debug_cleanup} ) {
    `rm -f $runpath/*.fas`;
    `rm -f $runpath/*.fas.accn*`;
    `rm -f $runpath/*key*`;
    `rm -f $runpath/*.aln*`;
    `rm -f $runpath/*concat*`;
    if ( -d "$runpath/model_test" ) {
        `rm -rf $runpath/model_test`;
    }
    `rm -f $runpath/*partition*`;
}

foreach my $sequence ( sort keys %files ) {
    foreach my $blastout ( @{ $files{$sequence}{'out'} } ) {
        my $fasfile = $blastout;
        $fasfile =~ s/.out/.fas/;
        if ( -s $fasfile ) {
            logger("FASTA file $fasfile already found. Skipping...\n");
            push( @{ $files{$sequence}{'fas'} }, $fasfile );
            next;
        } elsif ( -e $fasfile ) {
            `rm -f $fasfile`;
        }
        $options{cleanup} = 1;
        my $cov = 0;
        my $db  = '';
        if ( $blastout !~ /local/ ) {
            $cov = $options{coverage};
            $db  = 'nt';
        } elsif ( $blastout =~ /local/ ) {
            $cov = $options{local_cov};
            $db  = $blastout;
            $db  = $1 if $db =~ /_vs_([\S])+.local/;
        }

        my $command = join( " ",
                            $searchiopath, "-log", $logfile,
                            "-cov",        $cov,   "-accns",
                            $blastout,     ">",    $fasfile );
        my $exitcode = system($command);
        $exitcode = ( $exitcode >> 8 );

        if ( $exitcode == 255 ) {
            logger(
                "Search with query $sequence to database $db found no results. No sequences from $db will be included in the final analysis.\n"
            );
        } elsif ( $exitcode != 0 ) {
            die "Unable to generate fasta file from blast output!";
        } else {
            filecheck( "fasta file", $fasfile );
            push( @{ $files{$sequence}{'fas'} }, $fasfile );
        }
    }
}

$time = localtime();

logger("\nFinished extracting blast searches at $time.\n\n");
logger("---------------------------------------------------\n");
logger("---------------------------------------------------\n\n");

logger("Finding taxonomic and naming information.\n");
my %keys;
if ( -s "$runpath/all.keys" ) {
    open( my $keyfh, '<', "$runpath/all.keys" ) or die "Unable to open all.keys : $!\n";
    while(<$keyfh>){
        my $line = $_;
        my @data = split("\t",$line);
        $keys{$data[0]} = 1;
    }
}
my @searchfiles;
my @keyfiles;
foreach my $sequence ( sort keys %files ) {
    foreach my $fasout ( @{ $files{$sequence}{'fas'} } ) {
        my $keyfile = $fasout;
        $keyfile =~ s/\.fas$/.key/;
        if ( -s $keyfile ) {
            push(@keyfiles,$keyfile);
        }
        my $accnfile   = $fasout . ".accn.tmp";
        if ( -s $accnfile ) {
            push(@searchfiles,"$accnfile.dled");
            `mv $accnfile $accnfile.dled`;
        }
    }
}
chomp(@keyfiles);
my $kcompile = join(' ', @keyfiles);
if (@searchfiles) {
    my %searchids;
    my $files = join(" ",@searchfiles);
    my @tmpdata = `cat $files | sort | uniq`;
    foreach my $tmpid (@tmpdata) {
        if ( ! defined( $keys{$tmpid} ) ) {
            $searchids{$tmpid} = 1;
        }
    }
    if ( keys %searchids > 0 ) {
        open( my $searchfh, '>', "$runpath/all.accn" ) or die "Unable to open all.accn : $!\n";
        foreach my $accn (sort keys %searchids) {
            print $searchfh $accn;
        }
        close $searchfh;
        if ( ! -e "$runpath/keys.tmp" ) {
            system("touch $runpath/keys.tmp");
        }
        elink();
    }
} elsif ( -s "$runpath/all.accn" && ! -s "$runpath/keys.tmp" ) {
    #This occurs if the elink fails
    elink();
}
if ( -e "$runpath/all.keys" ) {
    `rm -f $runpath/all.keys`;
}
my $debug = 0;
if ($debug == 1) {
    print STDERR "Concatenating key files:\n";
    print STDERR "cat $kcompile $runpath/keys.tmp\n";
}

system("cat $kcompile $runpath/keys.tmp | sort | uniq > $runpath/all.keys") == 0 or die "Unable to compile key data : $!";

$time = localtime();

logger("\nFinished gathering taxonomic information at $time.\n\n");
logger("---------------------------------------------------\n");
logger("---------------------------------------------------\n\n");

my $file_debug = 1;
if ( $file_debug == 0 ) {
    foreach my $key ( keys %files ) {
        print STDERR "$key\n";
        foreach my $value ( @{ $files{$key}{'fas'} } ) {
            print "fas => $value\n";
        }
        foreach my $value ( @{ $files{$key}{'out'} } ) {
            print "out => $value\n";
        }
    }
}

$time = localtime();

logger(
    "Compiling and aligning fasta files started at $time; also trimming if desired\n"
);

if ( $options{cleanup} ) {
    &cleanup();
}

my $command = join( " ", $filtergenomespath, "-log", $logfile, "-key",
                    "$runpath/all.keys" );
my $skip_gfilter = 0;

foreach my $gene ( keys %files ) {
    if ( !-s "$runpath/$gene.all.fas.sorted" || !-s "$runpath/$gene.all.fas" ) {
        $skip_gfilter++;
    }
}


#concat & align files, one for each gene as input
if ( $skip_gfilter > 0 ) {
    &cleanup();
    foreach my $gene ( keys %files ) {
        my $fas_files = join( " ", @{ $files{$gene}{'fas'} } );
        if (    !-s "$runpath/$gene.all.fas.sorted"
             || !-s "$runpath/$gene.all.fas" )
        {
            `rm -f $runpath/$gene.all.fas` if -e "$runpath/$gene.all.fas";
            `rm -f $runpath/$gene.all.fas.sorted`
              if -e "$runpath/$gene.all.fas";
            system( "cat $fas_files" . ' > ' . "$runpath/$gene.all.fas" ) == 0
              or die
              "Unable to compile blast fasta files. Check permissions and try again.\n";
            $command .= " $runpath/$gene.all.fas";
        }
        push( @{ $files{$gene}{'cat'} }, "$runpath/$gene.all.fas" );
    }

    if ( defined $options{checkpoint} && $options{checkpoint} eq 'fasta'  ) {
        logger(
             "Halting before FASTA file compilation. Option -checkpoint fasta invoked.\n");
        exit();
    }

    logger("Finding genomes with copies of all genes\n");
    system("$command") == 0 or die "Unable to find genomes with all genes.\n";
}

foreach my $gene ( keys %files ) {
    filecheck( "genome filtered file", "$runpath/$gene.all.fas.sorted" );
    push( @{ $files{$gene}{'sorted'} }, "$runpath/$gene.all.fas.sorted" );
    $time = localtime();
    logger(
         "Beginning fasta alignment process for $gene.all.fas.sorted at $time\n"
    );
    if ( !( -s "$runpath/$gene.all.aln" ) ) {
        `rm -f $runpath/$gene.all.aln` if -e "$runpath/$gene.all.aln";
        die("Unable to remove old alignment file.  Check permissions and try again."
           )
          if $? != 0;
        my $command;
        if ( $options{align_prog} =~ /mafft|linsi|ginsi|einsi/ ) {
            $command = join(
                             " ",
                             $options{align_prog},
                             (
                                $options{threads}
                                ? "--thread $options{threads}"
                                : ''
                             ),
                             (
                                $options{align_params} ? $options{align_params}
                                : ''
                             ),
                             "--quiet",
                             "$runpath/$gene.all.fas.sorted > $runpath/$gene.all.aln"
                           );
        } elsif ( $options{align_prog} =~ /YOUR PROGRAM HERE/ ) {
            $command = '';    #PUT YOUR PROGRAM'S SPECIFIC COMMANDS HERE
        } else {
            $command = join( " ",
                            $options{align_prog}, $options{align_params},
                            "$runpath/$gene.all.fas > $runpath/$gene.all.aln" );
        }
        logger("Running command : $command\n");
        my $error = system("$command");
        if ( $error != 0 ) {
            `rm -f $runpath/$gene.all.aln`;
            die
              "Unable to align fasta files using $options{align_prog}.  Check to make sure this program is in your PATH or the full pathname is given.\n";
        }
    } else {
        logger("Alignment $gene.all.aln found. Skipping...\n");
        push( @{ $files{$gene}{'aln'} }, "$runpath/$gene.all.aln" );
        next;
    }
    push( @{ $files{$gene}{'aln'} }, "$runpath/$gene.all.aln" );
    $time = localtime();
    logger("Finished alignment for $gene.all.fas at $time\n");
}
$time = localtime();
logger("Finished all alignments using $options{align_prog} at $time\n\n");
logger("---------------------------------------------------\n");
logger("---------------------------------------------------\n\n");

#Perform Gblocks/noisy trim if desired

if ( $options{trimmer} ) {
    logger("Trimming alignments with $options{trimmer}...\n");
    foreach my $gene ( keys %files ) {
        foreach my $file ( @{ $files{$gene}{'aln'} } ) {
            if ( $options{trimmer} eq 'Gblocks' ) {
                if ( !-s "${file}-gb" ) {
                    `rm -f ${file}-gb` if -e "${file}-gb";
                    my $command = "$gblockspath $file $options{trimmer_params}";
                    logger("Running command : $command\n");
                    system("$command") == 256
                      or die "Unable to run Gblocks command properly : $command\n";
                } else {
                    logger("Gblocks file ${file}-gb already found. Skipping...\n");
                }

                filecheck( "Gblocks file", "${file}-gb" );
                push( @{ $files{$gene}{'trimmed'} }, "${file}-gb" );
            } elsif ( $options{trimmer} eq 'noisy' ) {
                my $trimmed = $file;
                $trimmed =~ s/\.aln/_out.fas/;
                if ( !-s "$trimmed" ) {
                    `rm -f $trimmed` if -e "$trimmed";
                    chdir("$runpath");
                    my $command = "$noisypath $options{trimmer_params} $file";
                    logger("Running command : $command\n");
                    system($command);
                    #my @output = `$command`; 
                    #logger(join("\n",@output));
                    if ( $? != 0 ) { 
                        logger("Unable to run noisy command properly : $command\n");
                        exit(-1);
                    }
                } else {
                    logger("noisy file $trimmed already found. Skipping...\n");
                }

                filecheck( "noisy file", $trimmed );
                push( @{ $files{$gene}{'trimmed'} }, $trimmed );
            }

        }
    }
    $time = localtime();
    logger("Finished trimming with $options{trimmer} at $time\n\n");
}
logger("---------------------------------------------------\n");
logger("---------------------------------------------------\n\n");

#Finding genes shared by all and concatenating

if ( -s "$runpath/$runid.concat" ) {
    logger(
        "Concatenated file $runpath/$runid.concat found.  Skipping concatenation.\n"
    );
} else {
    `rm -f $runpath/$runid.concat` if -e "$runpath/$runid.concat";
    logger(
        "Selecting isolates that contain all genes and concatenating sequences\n"
    );

    my @concat_files;
    foreach my $gene ( keys %files ) {
        if ( exists( $files{$gene}{'trimmed'} ) ) {
            foreach my $file ( @{ $files{$gene}{'trimmed'} } ) {
                push( @concat_files, $file );
            }
        } else {
            foreach my $file ( @{ $files{$gene}{'aln'} } ) {
                push( @concat_files, $file );
            }
        }
    }

    my $command = join( " ",
                        "cd $runpath;", $mlsaconcatpath,
                        "-k",           "$runpath/all.keys",
                        "-runid",       $runid,
                        "-log",         $logfile,
                        @concat_files,  ">",
                        "$runpath/$runid.concat" );

    system("$command") == 0 or die "Unable to concatenate alignment files\n";
}

${ $files{'output'} }{'concat'} = "$runpath/$runid.concat";
my $concat_out = "$runpath/$runid.concat";
filecheck( "concatenated alignments", $concat_out );

$time = localtime();
logger("Alignments concatenated. Finished at $time\n\n");
logger("---------------------------------------------------\n");
logger("---------------------------------------------------\n\n");

#logger("Concatenated file found at ${$files{'output'}}{'concat'}\n");

#Dereplicating sequences

my $noderep = 0;

if ( $noderep == 0 ) {
    logger("Removing identical sequences (dereplicating)...\n");

    if ( -e "$concat_out.derep" ) {
        logger("$concat_out.derep already found.  Skipping...\n");
    } else {
        my $command = join( " ", $dereplicatepath, $concat_out );
        system("$command") == 0 or die "Unable to dereplicate file. Exiting.\n";
    }
    filecheck( "dereplicated alignment", "$concat_out.derep" );
    ${ $files{'output'} }{'derep'} = "$concat_out.derep";
    $time = localtime();
    logger("Finished dereplicating files at $time.\n\n");
    logger("---------------------------------------------------\n");
    logger("---------------------------------------------------\n\n");
}

#If desired, rereplicate

if ( defined( $options{rereplicate} ) ) {
    if ( $options{rereplicate} eq '' ) {
        logger(
            "Stopping run so you can rereplicate the sequences you want.  Look in the $concat_out.derep.log file to see which sequences were removed.  Pull out those lines and supply to this script to rereplicate sequences.\n"
        );
        exit();
    } elsif ( -e "$options{rereplicate}" ) {
        if ( -e "$concat_out.derep.rerep" ) {
            logger("$concat_out.derep.rerep already found.  Skipping...\n");
        } else {
            logger("Rereplicating sequences...");
            my $command = join( " ",
                                $rereplicatepath, "-key", $options{rereplicate},
                                "$concat_out.derep" );
            system("$command") == 0
              or die "Unable to rereplicate file. Exiting.\n";
        }
    }
    filecheck( "rereplicated alignment", "$concat_out.derep.rerep" );
    ${ $files{'output'} }{'rerep'} = "$concat_out.derep.rerep";
    $time = localtime();
    logger("Finished rereplicating at $time\n\n");
    logger("---------------------------------------------------\n");
    logger("---------------------------------------------------\n\n");
}

#Begin with protein model selection.

if ( defined $options{checkpoint} && $options{checkpoint} eq 'model' ) {
    logger("Halting before model selection. Option --checkpoint model invoked.\n");
    exit();
}

my $proteinin;

if ( defined( ${ $files{'output'} }{'rerep'} ) ) {
    $proteinin = ${ $files{'output'} }{'rerep'};
} else {
    $proteinin = ${ $files{'output'} }{'derep'};
}

my ( $volume_junk, $path_junk, $proteinin_bare ) =
  File::Spec->splitpath($proteinin);

if ( $options{prog} eq 'tblastn' ) {

    logger(
        "Using concatenated and dereplicated/rereplicated file as input for protein model selection.\n"
    );

    if ( !-d "$runpath/model_test" ) {
        `mkdir $runpath/model_test`;
    }
    `cp $proteinin $runpath/model_test`
      unless ( -e "$runpath/model_test/$proteinin_bare" );
    `cp $runpath/mlsa.partition.tmp $runpath/model_test`
      unless ( -e "$runpath/model_test/mlsa.partition.tmp" );

    $command =
      join( " ", $proteinmodelpath, $proteinin_bare, $options{threads}, 'mlsa.partition.tmp' );

    if ( !-e "SPLIT_${proteinin_bare}_out" ) {
        system("cd $runpath/model_test;$command") == 0
          or die("Unable to split alignment files.");
    }
    my @splitaln = `cd $runpath/model_test;ls -1 ${runid}.concat*phy`;
    chomp(@splitaln);
    if ($? != 0) {
        logger("Unable to split alignment files.\n");
        logger("Check SPLIT_${proteinin_bare}_out in $runpath/model_test for more information!\n");
    }
    my %models;

    if ( -e "$runpath/model_test/model.log" ) {
        open MODELLOG, "$runpath/model_test/model.log"
          or die "Unable to open model log file : $!";
        while (<MODELLOG>) {
            my $line = $_;
            chomp($line);
            my ( $key, $value ) = split( "=", $line );
            $models{$key} = $value;
        }
        close MODELLOG;
    }

    foreach my $alignment (@splitaln) {
        my ( $vol, $dir, $file ) = File::Spec->splitpath($alignment);
        my $partition = $file;
        $partition =~ s/$proteinin_bare\.//;
        $partition =~ s/\.phy//;

        if ( defined( $models{$partition} ) ) {
            logger(
                 "Model for partition $partition already found. Skipping...\n");
            next;
        } else {
            my @partition_info = glob("$runpath/model_test/*$file*");
            if ( scalar(@partition_info) > 1 ) {
                `rm -f $runpath/model_test/*${file}_EVAL*`;
                `rm -f $runpath/model_test/ST_${file}_out`;
            }
        }
        open MODELLOG, ">>$runpath/model_test/model.log"
          or die "Unable to open model log file : $!";
        logger("Finding best model for partition $partition in file $file\n");
        my $command = join( " ", $proteinmodelpath, $alignment, $options{threads});
        my @proteinmodout = `cd $runpath/model_test;$command`;
        logger("Best Model for partition $partition : ");
        my $model = $1 if $proteinmodout[0] =~ /:(.*)$/;
        $model =~ s/ //;
        logger("$model\n");
        $models{"$partition"} = $model;
        print MODELLOG "$partition=$model\n";
        close MODELLOG;
    }

    open TMPMODEL, "$runpath/mlsa.partition.tmp"
      or die "Unable to open tmp model file : $!";
    open MODEL, ">$runpath/$runid.partition"
      or die "Unable to open partition file : $!";
    while (<TMPMODEL>) {
        my $line = $_;
        chomp($line);
        my ( $model_partition, $range ) = split( '=', $line );
        foreach my $key ( keys %models ) {
            if ( $model_partition =~ $key ) {
                print MODEL "$models{$key}, $key =$range\n";
                last;
            }
        }
    }
    close TMPMODEL;
    close MODEL;

    $time = localtime();
    filecheck( "protein model partition file", "$runpath/$runid.partition" );
    ${ $files{'output'} }{'partition'} = "$runpath/$runid.partition";
    logger("Protein model selection finished at $time\n\n");

} elsif ( $options{prog} = 'blastn' ) {
    logger("Generating partition file (with models) for nt analysis\n");

    open TMPMODEL, "$runpath/mlsa.partition.tmp"
      or die "Unable to open tmp model file : $!";
    open MODEL, ">$runpath/$runid.partition"
      or die "Unable to open partition file : $!";
    my $nt = 1;
    while (<TMPMODEL>) {
        my $line = $_;
        chomp($line);
        my ( $model_partition, $range ) = split( '=', $line );
        print MODEL "DNA, nt$nt =$range\n";
        $nt++;
    }
    close TMPMODEL;
    close MODEL;

    $time = localtime();
    filecheck( "nt model partition file", "$runpath/$runid.partition" );
    ${ $files{'output'} }{'partition'} = "$runpath/$runid.partition";
    logger("Partition file generation finished at $time\n");

}

logger("Use output files as input for autoMLSA-raxml.pl\n");
logger("Output files for run $runid:\n");
if ( -s $proteinin && -s "$runpath/$runid.partition" ) {
    logger(
            join( "\n",
                  "Alignment => " . $proteinin,
                  "Partition => $runpath/$runid.partition\n" )
          );
    logger("Model => PROTGAMMALG\n") if $options{prog} eq 'tblastn';
    logger("Model => GTRCAT\n")      if $options{prog} eq 'blastn';
} else {
    logger(
        "Delete $proteinin and $runpath/$runid.partition and try again.  Partition model and/or final alignment file are empty!\n"
    );
}

sub varcheck {

}

sub elink {
    my $retry = 1;
    while ($retry < 6) {
        my $command = join( " ", $elinkpath, "$runpath/all.accn", '-log', $logfile, '--email', $email, '>>', "$runpath/keys.tmp" );
        my $check = system($command);
        if ($check != 0) {
            if ($retry == 5) {
                logger("Unable to generate keyfile. Check your parameters and connection and try again.\n");
                exit(-1);
            }
            logger("Unable to generate keyfile on try $retry. Trying again.\n");
            $retry++;
        } else {
            last;
        }
    }
    `mv -f $runpath/all.accn $runpath/all.accn.dled`; 
}

sub logger {
    my $message = shift;
    print STDERR $message unless $options{quiet} == 1;
    if ( $log == 1 ) {
        unless ( -d $logpath ) {
            system("mkdir $logpath") == 0
              or die "Unable to make log dir. Check folder permissions.\n";
        }
        if ($runid) {
            open LOG, ">>$logfile" or die "$logfile is unavailable : $!";
            print LOG $message;
            close LOG;
        }
    }
}

sub filecheck {
    my $filetype = shift;
    my $filename = shift;
    if ( -s $filename ) {
        logger("Wrote $filetype to $filename\n");
    } else {
        die("Unable to generate $filetype! Check permissions and directory and try again.\n"
           );
    }
}

sub cleanup {
    `rm -f $runpath/*all.fas*`;
    `rm -f $runpath/*all.aln*`;
#    `rm -f $runpath/*all_*`;
    `rm -f $runpath/*concat*`;
    `rm -f $runpath/*partition*`;
    if ( -d "$runpath/model_test" ) {
        `rm -rf $runpath/model_test/`;
    }
    if ( -d "$runpath/trees" ) {
        `rm -rf $runpath/trees/`;
    }
}

sub gene_cleanup {
    my $gene = shift;
    $gene =~ s/\.out$/*/;
    if ( glob($gene) ) {
        `rm -f $gene`;
    }
}

#sub progress {
#    my $key = shift;
#    open PROGRESS, ">>$runid/$runid.progress" or die "$runid/$runid.progress is unavailable : $!";
#    print PROGRESS $key."\n";
#    close PROGRESS;
#}

__END__

=head1 NAME

autoMLSA.pl - Perform BLAST searches using provided input seqs, remove genomes without all genes, align and concatenate sequences

=head1 SYNOPSIS

autoMLSA.pl [options] -email (email@univ.edu) -runid (test_run) input[n].fasta input[n-1].fasta ... input[2].fasta input[1].fasta

=head1 OPTIONS

Defaults shown in square brackets.  Possible values shown in parentheses.

=over 8

=item B<-help|h>

Print a brief help message and exits.

=item B<-man>

Print a verbose help message and exits.

=item B<-quiet>

Turns off progress messages. Use noquiet to turn messages back on.

=item B<-email> (email@univ.edu) - B<REQUIRED>

Enter your email.  NCBI requires this information to continue. Can also set EMAIL environment variable instead.

=item B<-runid> (unique ID) - will generate directory for output - B<REQUIRED>

Name of the run. Allows for continuation of runs and keeps a log of the progress.

=item B<-config> (/path/to/config.txt)

Allows for input of a config file.  Command line parameters will supercede any parameters in the config file.  Useful for re-running an analysis with different parameters.  Will be generated prior to start of remote BLAST search.

=item B<-log|nolog> [logging on]

Using -nolog turns off logging. Logfile is in the format of runid.log, where runid is provided at command line. Placed in the ./logs/ dir.

=item B<-prog> [tblastn]

Sets blast program to run.  blastn and tblastn are valid options.

=item B<-complete|nocomplete> [on]

Option searches for complete records only.  Adds "AND complete genome" to searches that have an entrez query (supplied by -entrez_query 'text').

=item B<-remote> (nt) [off]

Searches the nt (or other valid e.g. refseq_genomic) remote database.

=item B<-local_db> (/path/to/blastdb1 /path/to/blastdb2)

Allows a search of one or more local blastdbs in addition to the remote databases.

=item B<-threads> [1]

Can run the local blast and mafft alignments in multithreaded manner, if desired. (Can specify threaded options for other alignment programs using -align_params)

=item B<-align_prog> [mafft-linsi]

The linsi algorithm of mafft is listed by default.  Any alignment program that prints to STDOUT is supported natively.  See body of script to insert specific code for other algorithms.

=item B<-align_params>

Supply a string in quotes that will be added to the end of the alignment program for program specific commands. (e.g. '--legacygappenalty' for mafft)

=item B<-trimmer> (gblocks OR noisy)

Run Gblocks or noisy to trim out gappy parts of the alignment.

=item B<-trimmer_params>

Can flag for trimmer to run with additional parameters using -trimmer_params='-b3=10 -b4=5 -b5=h' syntax at the commandline. noisy will automatically set seqtype parameter.

=item B<-rereplicate> (filename)

Supply a filename with tab delimited names to rereplicate prior to model selection and tree generation.  (Set the flag with no filename on first run to stop the run at this step.  Then pull out the genomes you want to replicate from the dereplicate log.

=item B<-cleanup>

Cleans up files and allows for addition of more genes to the analysis. This will delete previous results, so start a new analysis with a new runid if you are interested in saving your results.  Automatically cleans up old files when new genes are added to the analysis.  Previous alignments are not usable as some genomes may not contain all genes once more genes are added to the analysis.

=item B<-relaxed>

Decreases match and mismatch scores, as well as gap costs to allow for more distant matches for nucleotide searches.

B<PARAMETERS>

Except for evalue, by default the local values are set to equal the 'normal' parameters (e.g. -local_target = -target).

=over

=item B<-evalue> [1e-5] and B<-local_evalue> [1e-5]

Sets the e-value cutoffs.

=item B<-target> and B<-local_target> [500]

Sets the limit on number of alignments returned. 

=item B<-coverage> and B<-local_cov> [50]

Sets the threshold for coverage to include.  By default, includes all hits, regardless of coverage.  Values between 25 and 100 are valid [ex: 50 for alignment of 50% of the query and subject sequences]. 

=item B<-entrez_query>

Sets an entrez query limiter to blast search.  Useful for searching only within a species, genus, etc.

Example: For looking only at Actinobacteria, use 'Actinobacteria[ORGANISM]'

=back

B<ADVANCED>

=item B<-debug_cleanup>

Removes all files generated after BLAST search output.  Allows for restarting of jobs that may have failed and cannot continue due to unusual circumstances.

=item B<-checkpoint> (fasta OR model)

Stops script prior to filtering for genes found in all genomes or prior to alignment.  Useful for acquiring gene sequences without requiring completion of the entire pipeline, or doing the BLAST searches using a slower machine, saving the alignments for a more powerful machine. 

=back

=head1 DESCRIPTION

This program will read the given input file(s) and perform a blast search with 
the sequences.  The output includes alignment files for each input sequence
along with keyfiles needed to translate the headers to species names.

=cut
