#!/usr/bin/env perl

#######################################################################
#
# autoMLSA - A program to automate Multi-Locus Sequence Analysis
# Phylogenetic Trees using NCBI databases and local datasets
#
# This file is released with autoMLSA, but not as a part of autoMLSA.
# This file is derived from ProteinModelSelection.pl downloaded from
# http://sco.h-its.org/exelixis/resource/download/software/
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

#print $ARGV[0]." ".$ARGV[1]." ".$#ARGV."\n";

#adapt as required, modify to read ./raxmHPC if raxml is not in your path

$raxmlpath = '';    #Change to absolute path if raxml not found in PATH
@raxmlbins = ( "raxmlHPC", "raxmlHPC-SSE3", "raxmlHPC-PTHREADS",
               "raxmlHPC-PTHREADS-SSE3" );

#Change as desired, see raxmlbins above
$raxmlExecutable = $raxmlpath . $raxmlbins[3];

$alignmentName = shift;
$thread        = shift;
$partition     = shift;

if ( !$thread || $thread == 1 && $raxmlExecutable =~ /PTHREADS/ ) {
    print STDERR
      "You must supply a number of threads to use with this program when using PTHREADS version.\n";
    print STDERR "Check settings and try again.\n";
    exit(-1);
}

$raxmlExecutable .= " -T $thread" if $raxmlExecutable =~ /PTHREADS/;

$UNLIKELY = -1.0E300;

sub getLH {
    my $fileID = $_[0];
    open( CPF, $fileID );
    my @lines = <CPF>;
    close(CPF);
    my $numIT  = @lines;
    my $lastLH = pop(@lines);
    my $k      = index( $lastLH, '-' );
    my $LH     = substr( $lastLH, $k );
    return $LH;
}

sub getTIME {
    my $fileID = $_[0];
    open( CPF, $fileID );
    my @lines = <CPF>;
    close(CPF);
    my $numIT  = @lines;
    my $lastLH = pop(@lines);
    my $k      = index( $lastLH, '-' );
    my $TIME   = substr( $lastLH, 0, $k - 1 );
    return $TIME;
}

@AA_Models = (
               "DAYHOFF",  "DCMUT",  "JTT",       "MTREV",
               "WAG",      "RTREV",  "CPREV",     "VT",
               "BLOSUM62", "MTMAM",  "LG",        "MTART",
               "MTZOA",    "PMB",    "HIVB",      "HIVW",
               "JTTDCMUT", "FLU",    "DAYHOFFF",  "DCMUTF",
               "JTTF",     "MTREVF", "WAGF",      "RTREVF",
               "CPREVF",   "VTF",    "BLOSUM62F", "MTMAMF",
               "LGF",      "MTARTF", "MTZOAF",    "PMBF",
               "HIVBF",    "HIVWF",  "JTTDCMUTF", "FLUF",
               "LG4X"
             );

if ($partition) {
    print "Splitting up multi-gene alignment\n";
    $cmd =
        $raxmlExecutable
      . " -f s -m PROTCATJTT -p 12345 -s "
      . $alignmentName . " -q "
      . $partition
      . " -n SPLIT_"
      . $alignmentName
      . " \> SPLIT_"
      . $alignmentName . "_out";
    system($cmd);
    $count = 0;
    while ( open( CPF, $alignmentName . ".GENE." . $count ) ) {
        close CPF;
        print "PARTITION: " . $count . "\n";

        #print "perl ProteinModelSelection.pl ".$alignmentName.".GENE.".$count;
        system(   "perl ProteinModelSelection.pl "
                . $alignmentName
                . ".GENE."
                . $count );
        $count = $count + 1;
    }
} else {

 #print "Determining AA model data\n";
 #print "Computing randomized stepwise addition starting tree number :".$i."\n";
    $cmd =
        $raxmlExecutable
      . " -y -p 12345 -m PROTCATJTT -s "
      . $alignmentName
      . " -n ST_"
      . $alignmentName
      . " \> ST_"
      . $alignmentName . "_out";
    system($cmd);

    $numberOfModels = @AA_Models;

    for ( $i = 0 ; $i < $numberOfModels ; $i++ ) {
        $aa = "PROTGAMMA" . $AA_Models[$i];
        $cmd =
            $raxmlExecutable
          . " -f e -m "
          . $aa . " -s "
          . $alignmentName
          . " -t RAxML_parsimonyTree.ST_"
          . $alignmentName . " -n "
          . $AA_Models[$i] . "_"
          . $alignmentName
          . "_EVAL \> "
          . $AA_Models[$i] . "_"
          . $alignmentName
          . "_EVAL.out\n";

        #print($cmd);
        system($cmd);
    }

    for ( $i = 0 ; $i < $numberOfModels ; $i++ ) {
        $logFileName =
          "RAxML_log." . $AA_Models[$i] . "_" . $alignmentName . "_EVAL";

        #print $logFileName."\n";
        $lh[$i] = getLH($logFileName);
    }

    $bestLH = $UNLIKELY;
    $bestI  = -1;

    for ( $i = 0 ; $i < $numberOfModels ; $i++ ) {

        #print "Model: ".$AA_Models[$i]." LH: ". $lh[$i]."\n";
        if ( $lh[$i] > $bestLH ) {
            $bestLH = $lh[$i];
            $bestI  = $i;
        }
    }

    print "Best Model : " . $AA_Models[$bestI] . "\n\n";
    $bestModel = $AA_Models[$bestI];
}

