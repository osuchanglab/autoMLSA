#!/usr/bin/env perl
#######################################################################
#
# COPYRIGHT NOTICE
#
# autoMLSA-fasta_rename.pl - Auxillary script to rename FASTA entries
# without concatenating them.  Useful if you use autoMLSA.pl to 
# download individual sequences without needing filtering and 
# concatenation.
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
use Bio::SeqIO;

my $infile;
my $keyfile;
my %headers;
my $strain = 0;
my $bv = 0;
my %seen;
my $log = "header_rename.log";
my $header = 0;

GetOptions( 'keyfile=s' => \$keyfile,
            'strain'    => \$strain ,
            'bv'        => \$bv,
            'header'    => \$header);
if ($header == 0) {
    $infile = shift;
    die "Cannot find infile: $infile\n"   if ( !-e $infile );
}
die "Cannot find keyfile: $keyfile\n" if ( !-e $keyfile );

open KEYFILE, "<$keyfile" or die "Unable to open keyfile : $!";
#Current format is Accession,AssemblyID,TaxID,SciName,GuessName,GI,Master
while (<KEYFILE>) {
    my $line = $_;
    chomp($line);
    my ( $accn, @values ) = split( "\t", $line );
    my $header = '';
    my $sciname = $values[2];
    my $guess   = $values[3];
    my @test    = split( " ", $sciname );   #Test for quality of name from taxid
    if ( scalar(@test) == 2 ) {
        $header = $guess;
    } elsif ( scalar(@test) == 3 && $test[0] =~ /candidatus/i ) {
        $header = $guess;
    } elsif ( $bv == 1 && $sciname =~ /bv\./ && scalar(@test) == 4 ) {
        $header = $guess;
    } else {
        $header = $sciname;
    }

    $header =~ tr/ ()[]':/_{}{}__/;
    $header =~ s/,//;
    if ($strain) {
        $header =~ s/strain_//;
        $header =~ s/str\._//;
    }

    if (! exists( $seen{$header} ) ) {
        $seen{$header} = 1;
    } else {
        $seen{$header}++;
        $header .= "_$seen{$header}";
    }

    $headers{$accn} = $header;
}
close KEYFILE;

if ($header > 0) {
    foreach my $accn (sort keys %headers) {
        print join("\t",$accn,$headers{$accn})."\n";
    }
    exit(0);
}

open my $logfh, ">", $log or die "Unable to open logfile $log : $!";

#Setup new SeqIO stream
my $in = Bio::SeqIO->new( -file   => "$infile",
                          -format => 'fasta' );
my $out = Bio::SeqIO->new( -fh     => \*STDOUT,
                           -format => 'fasta' );
#Cycle through each sequence object in the SeqIO stream
while ( my $seq = $in->next_seq() ) {
    my $id  = $seq->id;
    if ( exists( $headers{$id} ) ) {
        print $logfh join("\t",$id,$headers{$id})."\n";
        $seq->id("$headers{$id}");
        $out->write_seq($seq);
    } else {
        $out->write_seq($seq);
    }

    print "\n";
}
