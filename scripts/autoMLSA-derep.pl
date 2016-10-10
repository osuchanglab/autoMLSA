#!/usr/bin/env perl
#######################################################################
#
# COPYRIGHT NOTICE
#
# autoMLSA-derep.pl - A script to automate dereplication of identical
# sequences after concatenation. Saves a log of dereplicated files that
# can be used to add in essential sequences (using autoMLSA-rerep.pl).
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
use Pod::Usage;

die "No infile specified.  Run with option -h for more information.\n" if scalar(@ARGV) == 0;

#initialize
my @infiles = @ARGV;
my $help = 0;
my $man = 0;

my $signal = GetOptions('help|h' => \$help,
			'man' => \$man);

pod2usage(-verbose => 1,
	  -output => *STDERR) if ($help == 1);
pod2usage(-verbose => 2,
	  -output => *STDERR) if $man == 1;

die("Unknown option given. Check parameters and try again.\n") if !$signal;

foreach my $infile (@infiles) {
    my $outfile = $infile;
    $outfile .= ".derep";
    my $logfile = $outfile.".log";
    my $in  = Bio::SeqIO->new(-file => "$infile",
			      -format => 'fasta');
    my %sequences;
    my %duplicates;
    while ( my $seq = $in->next_seq() ) {
	if (exists($sequences{$seq->seq})){
	    my $dupid = $sequences{$seq->seq};
	    push(@{$duplicates{$dupid}},$seq->id);
	} else {
	    $sequences{$seq->seq} = $seq->id;
	    $duplicates{$seq->id} = [];
	}
    }
    open OUTFILE, ">$outfile" or die "$outfile unavailable : $!";
    foreach my $seq (sort keys %sequences) {
	print OUTFILE ">".$sequences{$seq}."\n".$seq."\n";
    }
    close OUTFILE;
    open LOG, ">$logfile" or die "$logfile unavailable : $!";
    foreach my $dupid (sort keys %duplicates) {
	print LOG join("\t",$dupid,@{$duplicates{$dupid}})."\n";
    }
    close LOG;
}
