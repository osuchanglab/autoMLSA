#!/usr/bin/perl
#######################################################################
#
# COPYRIGHT NOTICE
#
# autoMLSA-rerep.pl - A script to re-replicate essential taxa after
# de-replication of identical sequences in a concatenated alignment.
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

my $keyfile;

my $signal = GetOptions('key=s' => \$keyfile);

die "No infile specified.  Run with option -h for more information.\n" if scalar(@ARGV) == 0;

die("Unknown option given. Check parameters and try again.\n") if !$signal;

my $infile = shift;
my $outfile = "$infile.rerep";

my %replicates;

open KEYFILE, "$keyfile" or die "$keyfile unavailable : $!";
while (<KEYFILE>) {
    my $line = $_;
    chomp($line);
    my ($key, @others) = split("\t",$line);
    $replicates{$key} = \@others;
}
close KEYFILE;

my $in  = Bio::SeqIO->new(-file => "$infile",
			  -format => 'fasta');

open OUTFILE, ">$outfile" or die "$outfile unavailable : $!";
while ( my $seq = $in->next_seq() ) {
    if (defined($replicates{$seq->id})){
	print OUTFILE ">".$seq->id."\n".$seq->seq."\n";
	foreach my $id (@{$replicates{$seq->id}}){
	    print OUTFILE ">".$id."\n".$seq->seq."\n";
	}
    } else {
	print OUTFILE ">".$seq->id."\n".$seq->seq."\n";
    }
}
close OUTFILE;
