#!/usr/bin/env perl
#######################################################################
#
# COPYRIGHT NOTICE
#
# autoMLSA-searchio.pl - A script to parse tabular BLAST output 
# produced by autoMLSA.pl. Extracts the accession numbers and aligned
# nucleotide or amino acid sequence from the results.
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

use warnings;
use strict;

use Getopt::Long;
use Pod::Usage;
use Bio::SeqIO;
use Bio::SearchIO;

my $logging;
my $quiet = 0;
my $cov;
my %memory;
my $accns = 0;
my %local;
my %above;
my %below;
my %scinames;
my %check;
my $badbase =  qr/[JOUBZ]/;

my $signal = GetOptions ( 'log:s' => \$logging,
                          'quiet' => \$quiet,
                          'cov=f' => \$cov,
                          'accns' => \$accns
                          );

if ( !$signal ){
    pod2usage (
        -verbose => 1,
        -message => "Unknown option entered. Check parameters and try again.\n",
        -exitval => -1
    );
}

my $infile = shift;
my $tmpfile = $infile;
$tmpfile =~ s/.out/.fas.accn.tmp/;
my $keyfile = '';
if ($infile =~ /local/) {
    $keyfile = $infile;
    $keyfile =~ s/.out/.key/;
}

if (! $infile ) {
    logger("No infile specified. Check parameters and try again.\n");
    logger("Syntax is $0 blast.out\n");
    die("\n");
}

if (! -e $infile) {
    logger("Unable to find input file $infile. Check parameters and try again.\n");
    die("\n");
}

if (defined($logging)) {
    if ($logging eq '') { 
        $logging = 'searchIO.log';
    }
}

$cov //= 0;

open my $fh, "<", "$infile" or die "Unable to open input file : $!";
my $temp;
if ($accns) {
    open $temp, ">", "$tmpfile" or die "Unable to open temp file : $!"; 
}
my $keyfileout;
if ($keyfile) {
    open $keyfileout, ">", "$keyfile" or die "Unable to open key file : $!";
}

my $complete = 0;

if ($infile !~ /local|wgs/ ) {
    $complete = 1;
}

my @fields;

while ( <$fh> ) {
    my $line = $_;
    chomp($line);
    if ($line =~ /hits/) {
        my $hits = $1 if $line =~ /([0-9]+)/;
        if ($hits == 0) {
            logger("No hits found in $infile. Exiting...\n");
            exit(-1);
        }
    }
    if ($line =~ /# Fields/) {
        @fields = split(',',$line);
    }
    next if ($line =~ /^#/);
    if (! @fields ) {
        logger("Unable to determine the fields in the BLAST output. -f 7 with custom fields is expected and required.\n");
        logger("Check your BLAST output to ensure the proper format is found.\n");
        exit(-1);
    }
    my @data = split("\t",$line);
    
    my ($query,$subject,$sacc,$pident,$qlen,$length,$evalue,$qcov,$taxid,$sciname,$stitle,$sseq) = 'NULL';

    for (my $i = 0; $i < scalar(@fields); $i++) {
        if ($fields[$i] =~ /query id/) {
            $query = $data[$i];
        } elsif ($fields[$i] =~ /subject id/) {
            $subject = $data[$i];
        } elsif ($fields[$i] =~ /subject acc/) {
            $sacc = $data[$i];
        } elsif ($fields[$i] =~ /identity/) {
            $pident = $data[$i];
        } elsif ($fields[$i] =~ /query length/) {
            $qlen = $data[$i];
        } elsif ($fields[$i] =~ /alignment length/) {
            $length = $data[$i];
        } elsif ($fields[$i] =~ /evalue/) {
            $evalue = $data[$i];
        } elsif ($fields[$i] =~ /query coverage per hsp/) {
            $qcov = $data[$i];
        } elsif ($fields[$i] =~ /subject tax ids/) {
            $taxid = $data[$i];
        } elsif ($fields[$i] =~ /subject sci names/) {
            $sciname = $data[$i];
        } elsif ($fields[$i] =~ /subject title/) {
            $stitle = $data[$i];
        } elsif ($fields[$i] =~ /subject seq/) {
            $sseq = $data[$i];
        }
    }
    my ($accn,$matched) = get_accn($sacc);

#    if ($accn =~ /\.[0-9]+/) {
#        $accn =~ s/\.[0-9]+//;
#    }

    if ( $sseq =~ $badbase ) {
        logger("Bad residue (J, O, or U) found in $accn. Check BLAST output file ($infile) for more information.\n");
        next;
    }
    if ($matched == 0) {
        if ($stitle =~ /\[/ && $stitle =~ /\]/) {
            $stitle = $1 if $stitle =~ /\[([^\]]+)\]$/;
        }
        $local{$accn} = 1;
    }
    if ( $qcov > $cov ) {
        $above{$accn} = 1;
        if ( !exists( $memory{$accn} ) ) {
            $memory{$accn} = 1;
        } else {
            logger("Multiple hits found for $accn, skipping.\n");
            next;
        }
        print ">" . $accn . "\n" . $sseq . "\n";
        if (! $local{$accn}) {
            if ($accns) {
                print $temp $accn . "\n";
            }
        } else {
#Current format is Accession,AssemblyID,TaxID,SciName,GI,Master,GenBankName,Country,Source,Strain,CultureCollection,Year
#                  $match        0        1      2    3    4         5         6       7      8          9           10
            print $keyfileout join("\t",$accn,'NULL','NULL',$stitle, 'NULL',$accn,'NULL','NULL','NULL','NULL','NULL', 'NULL')."\n";
        }
    } else {
        $below{$accn} = 1;
    }
}
#print STDERR "These accessions were above the cutoff:\n".join("\n",keys %above)."\n";
#print STDERR "These accessions were below the cutoff:\n".join("\n",keys %below)."\n";
close $fh;
close $temp;
if ($keyfile) {
    close $keyfileout;
}

if ( ! %above ) {
    logger("No hits above coverage cutoff for $infile. No FASTA file produced.\n");
    exit(-1);
}

sub logger {
    my $message = shift;
    print STDERR $message if $quiet == 0;

    if ($logging) {
        open my $log, ">>", "$logging" or die "Unable to open log file : $!";
        print $log $message;
        close $log;
    }
}

sub get_accn {
    my $name = shift;
    my @split = split( /\|/, $name );
    my @accn_types = (
                qr/\A[A-Z][0-9]{5}(\.[0-9])?\Z/,
                qr/\A[A-Z]{2}_?[0-9]{6}(\.[0-9])?\Z/,
                qr/\A([A-Z]{2}_)?[A-Z]{4}[0-9]{8}([0-9]{1,2})?(\.[0-9])?\Z/,
                qr/\A[A-Z]{2}_[A-Z]{2}[0-9]{6}(\.[0-9])?\Z/
    );
    my $match;
    my $matched = 0;
    foreach my $data (@split) {
        foreach my $type (@accn_types) {
            if ( $data =~ $type ) {
                $match = $data;
                $matched = 1;
            }
        }
    }
    if (!$match) {
        if ($split[1]) {
            $match = $split[1];
        } else {
            $match = $name;
        }
    }
    return ($match,$matched);
}
