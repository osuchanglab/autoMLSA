#!/usr/bin/perl
#######################################################################
#
# COPYRIGHT NOTICE
#
# autoMLSA-manual.pl - A script to fix genome names that were unable to
# be gathered properly using the NCBI E-utilities.
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

my $logging;
my $help;
my $quiet = 0;
my @infiles;

my $signal = GetOptions( 'log:s' => \$logging,
                         'help'  => \$help, );

@infiles = @ARGV;

foreach my $manual (@infiles) {
    if ( !-e $manual ) {
        logger(
            "Unable to find input file $manual. Check settings and try again.");
        die;
    }
    open my $manin, "<", $manual
      or die
      "Unable to open manual file $manual. Check settings and try again.";
    while (<$manin>) {
          my $line = $_;
          chomp($line);
          my @data    = split( "\t", $line );
          my $file    = $data[0];
          my $acc     = $data[1];
          my $orig    = $data[2];
          my $replace = $data[3];
          my $infile  = "$file.old";
          `mv $file $infile`;
          open my $fh, "<", $infile
            or die "Unable to open input file $infile : $!";
          open my $out, ">", $file
            or die " Unable to open output file for writing $file : $!";

          while (<$fh>) {
              my $line = $_;
              if ( $line =~ /$acc/ ) {
                  $line =~ s/\t${orig}\t/\t$replace\t/g;
                  logger("Replaced $orig with $replace in file $file\n");
              }
              print $out $line;
          }
#          `rm -f $infile`;
    }
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
