#!/usr/bin/env perl
#######################################################################
#
# COPYRIGHT NOTICE
#
# autoMLSA-raxml.pl - A Perl wrapper for running RAxML according to the
# RAxML v 8 manual.  Use autoMLSA.pl to generate the input files and 
# other required settings prior to use.
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
use File::Spec;
use List::Util qw(max);
use Scalar::Util qw(looks_like_number);

my $incommand = join( " ", $0, @ARGV );
pod2usage( -verbose => 0 ) if !@ARGV;

my $raxmldir = '';
my ( $runid, $log, $align, $model, $boot, $ml, $part, $ft );
$log = 1;
my $destroy = 0
  ; #Change this value to 0 if you don't want to lose partial RAxML run information.  !!You will have to remove the files manually if RAxML fails for any reason (including using Ctrl+C to kill the process) for the script to continue properly!!  It can be useful to disable deletion if you know you will be killing the script but want to keep the partial results.

my @raxmlbins = (
                  'raxmlHPC',               'raxmlHPC-SSE3',
                  'raxmlHPC-AVX',           'raxmlHPC-PTHREADS',
                  'raxmlHPC-PTHREADS-SSE3', 'raxmlHPC-PTHREADS-AVX'
                );

my %options = (
                'help'      => 0,
                'man'       => 0,
                'quiet'     => 0,
                'runid'     => \$runid,
                'log'       => \$log,
                'alignment' => \$align,
                'model'     => \$model,
                'partition' => \$part,
                'boottype'  => 'x'
              );

my $signal = GetOptions(
                         \%options,     'help|h',
                         'man',         'quiet',
                         'runid=s',     'log!',
                         'alignment=s', 'model=s',
                         'partition=s', 'initial',
                         'bootstrap:s', 'mlsearch:i',
                         'AVX',         'PTHREADS=i',
                         'SSE3',        'numperexec=i',
                         'clean:s',     'rearrange=i',
                         'boottype=s',  'final',
                         'cluster',     'lewis',
                         'outgroup=s'
                       );

#Print help statements if specified
pod2usage( -verbose => 1 ) if $options{help} == 1;
pod2usage( -verbose => 2 ) if $options{man} == 1;

die("Unknown option passed. Check parameters and try again.\n") if !$signal;
my @essentials = ( $runid, $align, $model );
foreach my $item (@essentials) {
    die("Please check parameters and resubmit. Missing essential input (runid, alignment file, AND/OR model).\n"
       )
      if !$item;
}

my @progcheck = ( "initial", "bootstrap", "mlsearch", "final" );
my $progcount = 0;
foreach my $prog (@progcheck) {
    $progcount++ if defined( $options{$prog} );
}

die("Check parameters and resubmit.  A single program is required (initial, mlsearch, bootstrap, or final).\n"
   )
  if $progcount != 1;

if ( defined( $options{AVX} ) && defined( $options{SSE3} ) ) {
    die("AVX and SSE3 cannot be both given at once.  Choose one or the other as appropriate for your architecture.\n"
       );
}

die("Cannot find alignment file.  Check parameters and try again.\n")
  if ( !-e $align );
if ( defined($part) ) {
    die("Unable to find partition file. Check path and try again.\n")
      if ( !-e $part );
}
die("Valid bootstrap types (-boottype) are x [default] and b.  Check parameters and try again.\n"
   )
  if $options{boottype} !~ /^[xb]$/;

my ( $runpath, $logpath, $logfile );
if ( -d "../$runid" ) {
    $runpath = File::Spec->rel2abs('./trees');
    $logpath = "../log/";
} else {
    if ( !-d "$runid" ) {
        `mkdir "$runid"`;
    }
    $runpath = File::Spec->rel2abs("$runid/trees");
    $logpath = "./log/";
}
$logfile = $logpath . "$runid.log";

if ( !-d "$runpath" ) {
    `mkdir "$runpath"`;
    die("Unable to make directory $runpath. Check folder permissions and try again.\n"
       )
      if $? != 0;
}

if ( defined( $options{bootstrap} ) ) {
    $boot = 1;
    $options{numperexec} //= $options{bootstrap};
    if ( $options{bootstrap} eq '' ) {
        $options{bootstrap} = 'autoMRE';
        if ( $options{numperexec} ) {
            die("numperexec is not compatible with a bootstrap setting of autoMRE [default bootstrap setting]. The auto bootstrap settings must be run on one PC. Check settings and try again.\n"
               );
        } else {
            $options{numperexec} = $options{bootstrap};
        }
    } elsif (
             $options{bootstrap} =~ /^autoFC|^autoMR$|^autoMRE$|^autoMRE_IGN$/ )
    {
        #Nothing to do here, just checking for valid options.
    } elsif ( looks_like_number( $options{bootstrap} )
              && ( $options{bootstrap} % $options{numperexec} ) )
    {
        die("Numbers of bootstraps ($options{bootstrap}) is not divisible by numperexec ($options{numperexec}). Check your numbers and try again.\n"
           );
    } elsif ( looks_like_number( $options{bootstrap} )
              && ( $options{bootstrap} > 0 && $options{numperexec} > 0 ) )
    {
        #Another validity check
    } else {
        die("Bootstrap input of $options{bootstrap} is not valid.  Valid options are integers > 0 and autoFC, autoMR, autoMRE, autoMRE_IGN\n"
           );
    }
}
if ( defined( $options{mlsearch} ) ) {
    $ml = 1;
    if ( $options{mlsearch} == 0 ) {
        $options{mlsearch} = 100;
    }
    $options{numperexec} //= $options{mlsearch};
    die("Numbers of ML searches ($options{mlsearch}) is not divisible by numperexec ($options{numperexec}). Check your numbers and try again.\n"
       )
      if $options{mlsearch} % $options{numperexec};
}

if ( !$boot && !$ml && !$options{initial} && !$options{final} ) {
    die("No procedure to run provided! Use -bootstrap, -mlsearch, -initial, or -final to choose the program.  Use flag -h for more information\n"
       );
}

my $bin = 0;
$bin += 1 if defined( $options{SSE3} );
$bin += 2 if defined( $options{AVX} );
$bin += 3 if defined( $options{PTHREADS} );
my $raxmlbin = $raxmldir . $raxmlbins[$bin];

logger("Command as submitted :\n$incommand\n");

my $time = localtime();
open SEED, ">$runpath/seed.log" or die "$runpath/seed.log is unavailable : $!";
if ( $options{initial} ) {
    logger("Starting initial rerearrangement setting optimization at $time\n");
    for ( my $i = 0 ; $i < 5 ; $i++ ) {
        $ft = "ST";
        my $seed = int( rand(100000) );
        print SEED "$runid.$ft$i\t$seed\n";
        my $command = join( " ",
                            $raxmlbin, "-y", "-n", "$runid.$ft$i", "-p", $seed,
                            "-s", $align, "-m", $model, "-w", $runpath );
        $command .= join( " ", " -q", $part ) if defined($part);
        $command .= join( " ", " -T", $options{PTHREADS} ) if defined( $options{PTHREADS} );
        $command .= join( " ", " --asc-corr=lewis" ) if $options{lewis};
        if ( -s "$runpath/RAxML_info.$runid.$ft$i" ) {
            logger("Run found for $runid.$ft$i. Skipping...\n");
            next;
        }
        &submit( \$ft, \$i, \$command );
    }
    $time = localtime();
    logger("Initial Parsimony trees generated at $time\n");
    logger("Starting automatic rearrangement tests\n");
    for ( my $i = 0 ; $i < 5 ; $i++ ) {
        $ft = "AI";
        my $seed = int( rand(100000) );
        print SEED "$runid.$ft$i\t$seed\n";
        my $command = join( " ",
                            $raxmlbin, "-f d",
                            "-n",      "$runid.$ft$i",
                            "-t", "$runpath/RAxML_parsimonyTree.$runid.ST$i",
                            "-s", $align,
                            "-m", $model,
                            "-w", $runpath );
        $command .= join( " ", " -q", $part ) if defined($part);
        $command .= join( " ", " -T", $options{PTHREADS} )
          if defined( $options{PTHREADS} );
        $command .= join( " ", " --asc-corr=lewis" ) if $options{lewis};
        my $check = &filecheck( \$ft, \$i );
        next if ($check);
        &submit( \$ft, \$i, \$command );
    }
    $time = localtime();
    logger("Automatic rearrangment tests completed at $time\n");
    logger("Starting fixed rearrangement tests\n");
    for ( my $i = 0 ; $i < 5 ; $i++ ) {
        $ft = "FI";
        my $seed = int( rand(100000) );
        print SEED "$runid.$ft$i\t$seed\n";
        my $command = join( " ",
                            $raxmlbin, "-f d",
                            "-n",      "$runid.$ft$i",
                            "-t", "$runpath/RAxML_parsimonyTree.$runid.ST$i",
                            "-s", $align,
                            "-m", $model,
                            "-w", $runpath,
                            "-i", "10" );
        $command .= join( " ", " -q", $part ) if defined($part);
        $command .= join( " ", " -T", $options{PTHREADS} )
          if defined( $options{PTHREADS} );
        $command .= join( " ", " --asc-corr=lewis" ) if $options{lewis};
        my $check = &filecheck( \$ft, \$i );
        next if ($check);
        &submit( \$ft, \$i, \$command );
    }
    $time = localtime();
    logger("Fixed rearrangement tests completed at $time\n");

    #   $rearrange = &rearrangement();
} elsif ($ml) {
    $ft   = "ML";
    $time = localtime();
    logger("Started maximum likelihood searches at $time\n");
    my $rearrange = &likelihood("I");
    $rearrange = $options{rearrange} if defined( $options{rearrange} );
    if ( !$rearrange ) {
        die("Unable to determine appropriate rearrangement setting. Try rerunning with -initial to determine appropriate rearrangment setting, or supply a rearrangment setting with -rearrange.\n"
           );
    }
    &clean($ft);
    for ( my $i = 0 ; $i < $options{mlsearch} ; $i += $options{numperexec} ) {
        my $seed = int( rand(100000) );
        print SEED "$runid.$ft$i\t$seed\n";
        my $command = join( " ",
                            $raxmlbin, "-f d",
                            "-n",      "$runid.$ft$i",
                            "-p",      $seed,
                            "-s",      $align,
                            "-m",      $model,
                            "-w",      $runpath,
                            "-i",      $rearrange,
                            "-N",      $options{numperexec} );
        $command .= join( " ", " -q", $part ) if defined($part);
        $command .= join( " ", " -T", $options{PTHREADS} )
          if defined( $options{PTHREADS} );
        $command .= join( " ", " --asc-corr=lewis" ) if $options{lewis};
        my $check = &filecheck( \$ft, \$i );
        next if ($check);
        &submit( \$ft, \$i, \$command );
    }
    $time = localtime();
    logger("Finished maximum likelihood searches at $time\n");
} elsif ($boot) {
    $ft   = "boot";
    $time = localtime();
    logger("Started bootstrap searches at $time\n");
    my $rearrange = &likelihood("I");
    $rearrange = $options{rearrange} if defined( $options{rearrange} );
    if ( !$rearrange ) {
        die("Unable to determine appropriate rearrangement setting. Try rerunning with -initial to determine appropriate rearrangment setting, or supply a rearrangment setting with -rearrange.\n"
           );
    }
    &clean($ft);
    logger(
          "Performing "
            . (
              looks_like_number( $options{bootstrap} ) ? $options{bootstrap} : 1
            )
            . " bootstrap runs, with "
            . (
                looks_like_number( $options{numperexec} )
                ? $options{numperexec} . ' per execution'
                : "runs until $options{numperexec} criterion is satisfied"
              )
            . "\n"
          );
    for (
          my $i = 0 ;
          $i <
          ( looks_like_number( $options{bootstrap} ) ? $options{bootstrap} : 1 )
          ;
          (
            looks_like_number( $options{numperexec} ) ? $i +=
              $options{numperexec} : $i++
          )
        )
    {
        my $seed     = int( rand(100000) );
        my $bootseed = int( rand(100000) );
        print SEED "$runid.$ft$i\t$seed\tboot\t$bootseed\n";
        my $command = join( " ",
                            $raxmlbin,                "-f d",
                            "-n",                     "$runid.$ft$i",
                            "-p",                     $seed,
                            "-s",                     $align,
                            "-m",                     $model,
                            "-w",                     $runpath,
                            "-i",                     $rearrange,
                            "-N",                     $options{numperexec},
                            "-" . $options{boottype}, $bootseed );
        $command .= join( " ", " -q", $part ) if defined($part);
        $command .= join( " ", " -T", $options{PTHREADS} )
          if defined( $options{PTHREADS} );
        $command .= join( " ", " --asc-corr=lewis" ) if $options{lewis};
        my $check = &filecheck( \$ft, \$i );
        next if ($check);
        &submit( \$ft, \$i, \$command );
    }
    $time = localtime();
    logger("Finished bootstrap searches at $time\n");
} elsif ( $options{final} ) {
    $ft   = "ML";
    $time = localtime();
    logger(
        "Started final tree preparation (choosing best ML tree and applying boostrap support) at $time\n"
    );
    my $tree = &likelihood($ft);
    if ( !$tree ) {
        die("Unable to deteremine best ML tree.  Check to ensure you have run $0 -mlsearch and have results in the $runpath folder.\n"
           );
    }
    my $bootfile;
    my @bootstraps = `ls $runpath/RAxML_bootstrap.$runid.boot*`;
    chomp(@bootstraps);
    if ( scalar(@bootstraps) == 0 ) {
        die("Unable to find any bootstraps to apply to best ML tree. Check to ensure you have run $0 -bootstrap and have results in the $runpath folder.\n"
           );
    } elsif ( scalar(@bootstraps) > 1 ) {
        $bootfile = "$runpath/$runid.boot.all";
        `rm -f $bootfile` if ( -e $bootfile );
        logger("Detected multiple boostrap files.  Concatenating results.\n");
        my $command = "cat " . join( " ", @bootstraps ) . " > $bootfile";
        &submit( '', '', \$command );
        if ( -e $bootfile ) {
            logger("Concatenated bootstrap results found at $bootfile\n");
        } else {
            die("Unable to find concatenated bootstrap file. Check parameters and try again.\n"
               );
        }
    } else {
        $bootfile = $bootstraps[0];
        if ( -e $bootfile ) {
            logger("Bootstrap results found at $bootfile\n");
        } else {
            die("Unable to find concatenated bootstrap file. Check parameters and try again.\n"
               );
        }
    }
    $ft = 'final';
    &clean($ft);
    logger("Applying bootstraps to best tree.\n");
    my $i       = '.nwk';
    my $command = join( " ",
                        $raxmlbin,          "-f b", "-n",
                        "$runid.final.nwk", "-s",   $align,
                        "-m",               $model, "-t",
                        $tree,              "-z",   $bootfile,
                        "-w",               $runpath );
    $command .= join( " ", " -T", $options{PTHREADS} )
      if defined( $options{PTHREADS} );
    $command .= join( " ", " --asc-corr=lewis" ) if $options{lewis};
    $command .= join( " ", " -o",$options{outgroup} ) if $options{outgroup};
    my @outfiles = (
                     "$runpath/RAxML_bipartitions.$runid.final.nwk",
                     "$runpath/RAxML_bipartitionsBranchLabels.$runid.final.nwk",
                     "$runpath/RAxML_info.$runid.final.nwk"
                   );
    foreach my $outfile (@outfiles) {
        if ( -e $outfile ) {
            `rm -f $outfile`;
        }
    }
    &submit( \$ft, \$i, \$command );
    logger(   "Congratulations!!  Output files are found at these locations:\n"
            . join( "\n", @outfiles[0..1] )
            . "\n" );
    $time = localtime();
    logger("Run completed at $time\n");
}

close SEED;

sub likelihood {
    my $ft = shift;
    my $score;
    my @scores = `grep -a 'Final GAMMA' $runpath/*_info.$runid.*$ft*`;
    my %scores;
    if ( scalar(@scores) > 1 ) {
        foreach my $item (@scores) {
            my ( $file, $score ) = split( ":", $item );
            my @split = split( " ", $score );
            $scores{$file} = $split[6];
        }
        @scores = ( values %scores );
        my $max = max @scores;
        foreach my $key ( keys %scores ) {
            if ( $scores{$key} eq $max ) {
                if ( $ft eq "I" ) {
                    if ( $key =~ /FI/ ) {
                        logger(
                              "Best score found with rearrangement of -i=10\n");
                        $score = 10;
                        return $score;
                    } else {
                        my $score = `grep -a 'best rearrangement' $key`;
                        $score = $1 if $score =~ / ([\d]+)$/;
                        logger(
                            "Best score found using rearrangement -i=$score\n");
                        return $score;
                    }
                } elsif ( $ft eq "ML" ) {
                    $key =~ s/info/bestTree/;
                    logger("Best ML tree found : $key\n");
                    return '' if ( !-e $key );
                    return $key;
                } else {
                    die("Forgot to provide file type to likelihood subroutine. Fix your coding.\n"
                       );
                }
            }
        }
    } else {
        if ( $ft eq "ML" ) {
            my $file = "$runpath/RAxML_bestTree.$runid.ML0";
            logger("Best ML tree found : $file\n");
            return '' if ( !-e $file );
            return $file;
        }
    }

}

sub clean {
    my $ft = shift;
    if ( defined( $options{clean} ) ) {
        if ( $options{clean} =~ /force/ ) {
            logger("!!Deleting all $ft tree files.!!\n");
            `rm -f $runpath/*$runid.$ft*`;
        } elsif ( $options{clean} eq '' ) {
            print STDERR
              "Are you sure you want to delete ALL $ft tree files?  This cannot be undone [yes/no]\n";
            my $response = <STDIN>;
            if ( $response =~ /yes/i ) {
                logger("!!Deleting all $ft tree files.!!\n");
                `rm -f $runpath/*$runid.$ft*`;
            } else {
                logger(
                     "Not deleting tree files. Check parameters and resubmit.");
                die("\n");
            }
        }
    }
}

sub filecheck {
    my $ft    = shift;
    my $i     = shift;
    my @files = `ls $runpath/`;
    my $term;
    if ( $$ft eq 'ST' ) {
        $term = "info.$runid.$$ft$$i";
    } elsif ( $$ft eq 'ML' ) {
        $term = "parsimonyTree.$runid.$$ft$$i";
    } elsif ( $$ft ne 'boot' ) {
        $term = "log.$runid.$$ft$$i";
    } else {
        $term = "bootstrap.$runid.$$ft$$i";
    }
    if ( my @matched = grep $_ =~ /$term/, @files ) {
        foreach my $file (@matched) {
            chomp($file);
            if ( -s "$runpath/$file" ) {
                logger("Run found for $runid.$$ft$$i. Skipping...\n");
                return 1;
            } elsif ( -e "$runpath/RAxML_info.$runid.$$ft$$i" ) {
                `rm -f $runpath/*$runid.$$ft$$i*`;
            }
        }
    }
}

sub submit {
    my $ft      = shift;
    my $i       = shift;
    my $command = shift;
    logger("Submitting command:\n $$command\n");
    my $error = system("$$command");
    if ( $error != 0 ) {
        if ( !defined($ft) || !defined($i) ) {
            die("Unable to run command properly (Exited with code $error). Check parameters/files/permissions and try again.\n"
               );
        } else {
            `rm -f $runpath/*$runid.$$ft$$i*` if $destroy == 1;
            die("Unable to finish RAxML run ($runid.$$ft$$i) properly (Exited with code $error). Check parameters/files/permissions and try again.\n"
               );
        }
    }
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
            unless ( -e "$logfile" ) {
                system("touch $logfile") == 0
                  or die
                  "Unable to make log file $logfile Check folder permissions.\n";
            }
            open LOG, ">>$logfile" or die "$logfile is unavailable : $!";
            print LOG $message;
            close LOG;
        }
    }
}

__END__
=head1 NAME

raxml_wrapper.pl - Wrapper to run RAxML to generate maximum-likelihood based phylogenetic trees

=head1 SYNOPSIS

raxml_wrapper.pl -runid runid -alignment runid.aln -model PROTGAMMALG -partition runid.partition -SSE3 -PTHREADS threads -(program)

=head1 OPTIONS

Defaults shown in square brackets.  Possible values shown in parentheses.

=over 8

=item B<-help|h>

Print a brief help message and exits.

=item B<-man>

Print a verbose help message and exits.

=item B<-quiet>

Turns off progress messages.

=item B<-log|nolog> [logging on]

Using -nolog turns off logging. Logfile is in the format of runid.log, where runid is provided at command line. Placed in the ./logs/ dir.

=item B<-runid> (unique ID) - will generate directory for output - B<REQUIRED>

Name of the run. Allows for continuation of runs and keeps a log of the progress.

=item B<-alignment> (fasta or phylip alignment file) - B<REQUIRED>

Path to alignment file to use for input to RAxML.

=item B<-model> (protein or DNA substitution model) - B<REQUIRED>

Required, even if a partition file is provided.  The partition file will override this model, but it is required to tell RAxML if the input is protein (e.g. PROTGAMMALG) or DNA (e.g. GTRGAMMA).

=item B<-partition> (RAxML partition file)

Path to partition file to tell RAxML to perform a partitioned analysis.

B<PROGRAMS> - One (and only one) is required. Should be run in this order.

=over

=item B<-initial>

This program determines the best initial rearrangement setting. Helps speed up the mlsearch and bootstrap programs.

=item B<-mlsearch> [100] (integer {1 or greater})

This program performs a maximum likelihood search (-f d) with the number of specified starting trees.

=item B<-bootstrap> [autoMRE] (integer {1 or greater} or RAxML supported bootstopping criterion {autoMR, autoMRE, autoMRE_IGN, autoFC})

This program performs a bootstrapping analysis of the data, with the specified number of bootstraps or bootstopping criterion specified.

=item B<-final>

This program selects the best tree from the mlsearch and applies the bootstrapping support to it.

=back

B<PARAMETERS>

=over

=item B<-PTHREADS> [0] (integer > 1)

If the PTHREADS version is specified, the number of threads to use is also required.  Not useful for the -final program.

=item B<-SSE3> [off]

Computers made in the last 10 years should support this.  Use this if you are unsure if your computer supports AVX.

=item B<-AVX> [off]

Computers made in the last ~3-4 years should support this.  Only use this if you are sure your computer supports it.

=item B<-numperexec> [# of ML searches or bootstraps specified] (any factor of # of ML searches or bootstraps is supported)

Use this if you are interested in running multiple analyses on multiple computers.  You can submit the same command to multiple computers (through a system like SGE) and the script should handle it properly.  If you change the numperexec for different runs, the outcome is unexpected, so keep the numperexec the same for multiple runs.

=item B<-clean> (force)

If you are having issues with running your script and want to clear out old results, use this command.  Only recommended if you KNOW you want to delete previously generated data.  !!THIS CANNOT BE UNDONE SO DO NOT DO THIS LIGHTLY!!  Using -clean force will not prompt you in the script for deletion.

=item B<-rearrange> (integer)

If you know what rearrangement setting you want to use and have not done the initial program (or if you are generating a second tree from the same data and want to try again with different RAxML settings) you can supply it here.

=item B<-boottype> [x] (x or b)

Sets the bootstrapping type for RAxML.  By default, rapid bootstrapping is selected [x].

=item B<-lewis> [off]

Apply ascertainment bias correction with the standard correction by Paul Lewis.

=item B<-outgroup "OUTGROUP">

Supply an outgroup for rooting the final tree.

=back

=back

=head1 DESCRIPTION

This script is simply a wrapper for the most commomn RAxML commands to facilitate generation of high quality phylogenetic trees. Works well with output from remote_blast.pl.

=cut
