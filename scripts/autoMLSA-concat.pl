#!/usr/bin/perl
#######################################################################
#
# COPYRIGHT NOTICE
#
# autoMLSA-concat.pl - A script to concatenate a series of genes used
# for a Multi-Locus Sequence Alignment.  Uses a keyfile to translate 
# the accession numbers to strain names and also generates a partition
# file to split the alignment apart later.
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
use Bio::SeqIO;
use Getopt::Long;
use Cwd 'abs_path';
use File::Spec;

my $keyfile;
my %headers;
my @genes;
my $i = 0;
my @lengths;
my $help;
my $runid;
my $logging;
my $quiet = 0;
my %idmatch;
my $strain = 1;

#Get options
my $signal = GetOptions(
                         'help|h'    => \$help,
                         'key|k=s'   => \$keyfile,
                         'runid|r=s' => \$runid,
                         'log:s'     => \$logging,
                         'quiet'     => \$quiet,
                         'strain!'   => \$strain
                       );

#Die if GetOptions returns false [unknown option passed]
die("Unknown option entered.  Run command with option -h for available options.\n"
   )
  if !$signal;

if ( scalar(@ARGV) < 2 || $help || !defined($keyfile) ) {
    die(
        "Proper syntax is $0 -k headers.keyfile alignment1.fas alingment2.fas\n > concat.fas\n",
        "More than two alignment files is acceptable.\n"
       );
}

my @infiles = @ARGV;

if ( defined($logging) ) {
    if ( $logging eq '' ) {
        $logging = "mlsa_concat.log";
    }
}

open KEYFILE, "$keyfile" or die "$keyfile unavailable : $!";

#Load hash with keyfile data
#Current format is Accession,AssemblyID,TaxID,SciName,GuessName,GI,Master,GenBankName,Host,Country,Source,Strain,CultureCollection
while (<KEYFILE>) {
    my $line = $_;
    chomp($line);
    my ( $accn, @values ) = split( "\t", $line );
    my $match = $values[0];
    if ( $match eq 'NULL' ) {
        $match = $values[5];
    }
    $idmatch{$accn} = $match;

    if ( !exists( $headers{$match} ) ) {
    
        #Choose appropriate header name
        my $sciname = $values[2];
        my $guess   = $values[3];
        my $gbname = '';
        my $strain = '';
        my $culture = '';
        if ($values[6]) {
            $gbname = $values[6];
        }
        if ($values[10]) {
            $strain = $values[10]; 
        }
        if ($values[11]) {
            $culture = $values[11];
        }
        my $header = $sciname;
        my @test = split( " ", $sciname );  #Test for quality of name from taxid
        my $test = @test;
        my @test2 = split( " ", $gbname );
        my $test2 = @test2;
        my $candidatus = 0;
        my $subsp = 0;
        if ($sciname =~ /candidatus/i ) {
            $candidatus = 1;
        }
        if ($sciname =~ /pv\.|bv\.|subsp\./i ) {
            $subsp = 2;
        }
        my $value = 2 + $candidatus + $subsp;
        if ( $test == $value ) {
            $header = $guess;
            if ( $gbname ) {
                if ($sciname ne $gbname) {
                    $header = $gbname;
                }
            }
        } 

        $headers{$match} = $header;
    }
}

close KEYFILE;

my @filenames;

foreach my $infile (@infiles) {

    my ($volume, $dir, $filename) = File::Spec->splitpath( $infile );

    $filename =~ s/\.[^\.]+//;

    $filenames[$i] = $filename;

    #Setup new SeqIO stream
    my $in = Bio::SeqIO->new( -file   => "$infile",
                              -format => 'fasta' );

    #Cycle through each sequence object in the SeqIO stream
    while ( my $seq = $in->next_seq() ) {
        my $sequence = $seq->seq;
        my $id       = $seq->id;
        my $match    = $idmatch{$id};
        $genes[$i]{$match} = $sequence;
        $lengths[$i]{ $seq->length } = 1;
    }
    $i++;
}

my %memory;

open( HEADERLOG, ">",
      ( defined($runid) ? "$runid.concat.log" : "mlsa_concat.log" ) )
  or die "Unable to open concatenation log file : $!";

foreach my $key ( sort keys %{ $genes[0] } ) {
    my $header = $headers{$key};
    $header =~ tr/ ()[]':/_{}{}__/;
    $header =~ s/,//;
    if ($strain) {
        if ($header =~ /substr\._/) {
            $header =~ s/substr\._//;
        }
        $header =~ s/strain_//;
        $header =~ s/str\._//;
    }
    if ( exists( $memory{$header} ) ) {
        $memory{$header}++;
        $header .= "_".$memory{$header};
    } else {
        $memory{$header} = 0;
    }
    print ">" . $header . "\n";
    print HEADERLOG join( "\t", $header, $key ) . "\n";
    for ( my $j = 0 ; $j < scalar(@infiles) ; $j++ ) {
        if ( exists( $genes[$j]{$key} ) ) {
            print $genes[$j]{$key};
        } else {
            &log("$key.\n");
            &log(   "Headers not equal in all "
                  . scalar(@infiles)
                  . " files... exiting" );
            die("\n");
        }
    }
    print "\n";
}

my $total_length = 0;

my $partfile = "mlsa.partition.tmp";
for ( my $k = 0 ; $k < scalar(@infiles) ; $k++ ) {
    if ( keys %{ $lengths[$k] } > 1 ) {
        die
          "Gene lengths not equal for gene $infiles[$k].  Re-align your sequences.\n";
    } else {
        if ( $k == 0 && -e $partfile ) {
            `rm -f $partfile`;
        }
        open PARTITION, ">>$partfile" or die "$partfile is unavailable : $!";
        my @length = keys %{ $lengths[$k] };
        &log("$infiles[$k] gene length is $length[0].\n");

    #	open GREP, "grep \'Best model according to BIC\' $infiles[$k].prottest |";
    #	my $model = <GREP>;
    #	chomp($model);
    #	close GREP;
    #	    &log $model."\n";
    #	$model = $1 if $model =~ /\: ([^\s]+$)/;
    #	    &log "\t".$model."\n";
        my $model = 'LGF'
          ;    #placeholder for now.  New system set up to get models from RAxML
        print PARTITION $model . ", " 
          . $filenames[$k]
          . " = "
          . ( $total_length + 1 ) . "-"
          . ( $total_length + $length[0] ) . "\n";
        close PARTITION;
        $total_length += $length[0];
    }
}

sub log {
    my $message = shift;
    print STDERR $message unless $quiet == 1;

    #    if ($logging) {
    #	open LOG, ">>$logging" or die "$logging is unavailable : $!";
    #	print LOG $message;
    #	close LOG;
    #   }
}
