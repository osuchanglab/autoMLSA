#!/usr/bin/env perl
#######################################################################
#
# COPYRIGHT NOTICE
#
# autoMLSA-filter.pl - A script to automate filtering of a set of genes
# to include only genomes with complete sets of genes. Primarily uses
# the unique assembly IDs to identify genes in genomes with multiple
# replicons.
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
use Bio::SeqIO;
use Getopt::Long;

my $logging;
my $help;
my $quiet = 0;
my $keyfile;
my %headers;

#Get options
my $signal = GetOptions(
                         'log:s' => \$logging,
                         'key=s' => \$keyfile,
                         'help'  => \$help,
                         'quiet' => \$quiet
                       );

die("Unknown option entered.  Run command with option -h for available options.\n"
   )
  if !$signal;

if ( scalar(@ARGV) < 1 || $help ) {
    die( "Proper syntax is $0 alignment1.fas alingment2.fas\n",
         "More than two alignment files is required.\n" );
}
if ( !$keyfile ) {
    die("A keyfile must be supplied with -key! Check parameters and try again\n"
       );
} elsif ( !-e "$keyfile" ) {
    die("Unable to find keyfile $keyfile.\n Check parameters and try again\n");
}

if ( defined($logging) ) {
    if ( $logging eq '' ) {
        $logging = "filter.log";
    }
}

my @infiles = @ARGV;
my %accns;
my @hashArray;
my %wgs_seqs;
my %names;
my %all_g;
my $g = 0;
my $i = 0;

#Read keyfile

open KEYFILE, "$keyfile" or die "$keyfile unavailable : $!";

#Load hash with keyfile data
#Current format is Accession,AssemblyID,TaxID,SciName,GuessName,GI,Master
#Current format is Accession,AssemblyID,TaxID,SciName,GI,Master,GenBankName,Country,Source,Strain,CultureCollection,Year
#                  $match        0        1      2    3    4         5         6       7      8          9           10
while (<KEYFILE>) {
    my $line = $_;
    chomp($line);
    my ( $accn, @values ) = split( "\t", $line );
    $headers{$accn}{'assemid'}   = $values[0];
    $headers{$accn}{'taxid'}     = $values[1];
    $headers{$accn}{'sciname'}   = $values[2];
    $headers{$accn}{'gi'}        = $values[3];
    $headers{$accn}{'master'}    = $values[4];
    if ($values[0] eq 'NULL') {
        $names{$accn} = $values[2];
    } else {
        $names{$values[0]} = $values[2];
    }
}

close KEYFILE;

#Read input files
foreach my $input (@infiles) {
    my $in = Bio::SeqIO->new( -file   => "$input",
                              -format => 'fasta' );
    while ( my $seq = $in->next_seq() ) {

        #Set up individual hashes for each gene
        my $id      = $seq->id;
        
        if (! exists($headers{$id}) ) {
            logger("No keyfile information found for id $id, skipping\n");
            next;
        }

        my $assemid = $headers{$id}{'assemid'};
        
        if (!$assemid) {
            die("No assembly id/other id found for accession number $id\n");
        }

        if ( $assemid eq 'NULL' ) {
            $assemid = $headers{$id}{'master'};
        }
        my $header = $assemid;

        if ( defined( $hashArray[$g]{$header} ) ) {

            #	    print STDERR "$header already seen. Skipping...\n";
            next;
        }
        $accns{$header} = $id;
        $hashArray[$g]{$header} = $seq->seq;
        $i++;
    }

    #    print STDERR "$i genes counted for genome $g\n";
    $i = 0;
    $g++;
}
$g = 0;

#Debugging step
#foreach my $key (sort keys (%all_g)){
#    foreach my $value (@{$all_g{$key}}){
#	print STDERR $value."\t";
#    }
#    print "\n";
#}
my @save = keys( %{ $hashArray[0] } );
my %save;
my %remove;
foreach my $key (@save) {
    $save{$key} = 1;
}
my $count = keys %save;
logger("Starting search with $count keys (searching file $infiles[0])\n");
for ( $g = 1 ; $g < scalar(@hashArray) ; $g++ ) {
    foreach my $key ( keys(%save) ) {
        if ( !defined( $hashArray[$g]{$key} ) ) {
            $remove{$key} = 1;
            delete $save{$key};
        }
    }
    $count = keys %save;
    logger("$count keys remaining after searching gene $g (file $infiles[$g])\n");
}
$count = keys %save;

logger( "$count genomes have all " . scalar(@hashArray) . " genes\n" );
$count = keys %remove;
logger( "$count genomes did not have all " . scalar(@hashArray) . " genes\n" );
logger( "These genomes did not have all genes and were removed:\n" );
foreach my $id (sort keys %remove) {
    logger("$id => $names{$id}\n");
}
$g = 0;

while ( $g < scalar(@hashArray) ) {
    open OUTFILE, ">$infiles[$g].sorted";
    foreach my $gene_hash ( $hashArray[$g] ) {
        my $i = 0;
        foreach my $key ( sort sortHash keys(%save) ) {
            print OUTFILE ">$accns{$key}\n$gene_hash->{$key}\n";
            $i++;

        }

        #	print STDERR "$i hashes sorted\n";
    }
    close OUTFILE;
    $g++;
}

#Debugging step
#foreach my $array (@gene_hashes){
#    foreach my $gene (@{$array}){
#	print $gene."\t";
#    }
#    print "\n";
#}

sub sortHash {
    $a cmp $b;
}

sub sortHashwgs {
    my ( $junk1, $code1 ) = split( '\|', $a );
    my ( $junk2, $code2 ) = split( '\|', $b );
    $code1 cmp $code2;
}

sub logger {
    my $message = shift;
    print STDERR $message unless $quiet == 1;
    if ($logging) {
        open LOG, ">>$logging" or die "$logging is unavailable : $!";
        print LOG $message;
        close LOG;
    }
}

