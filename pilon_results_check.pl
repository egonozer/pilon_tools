#!/usr/bin/perl

use strict;
use warnings;
$|++;

my $usage = "
pilon_results_check.pl

Looks for any leftover indels after pilon correction of an assembly. Could
be used beforehand, too.

Required:
  -b    alignment.bam file. Alignment file must be sorted and indexed.
  -f    reference fasta file.

";

use Getopt::Std;
our ($opt_b, $opt_f);
getopts('b:f:');

die $usage unless ($opt_b and $opt_f);
my $samfile = $opt_b;
my $reffile = $opt_f;

 my $min_h_leng = 3;
 my $min_frac = 20;

## Read in the reference sequence
my %refseq;
my ($id, $seq);
open (my $in, "<$reffile") or die "ERROR: Can't open $reffile: $!\n";
while (my $line = <$in>){
    chomp $line;
    if ($line =~ m/^>/){
        if ($id){
            $refseq{$id} = $seq;
        }
        $seq = "";
        $id = substr($line, 1);
        $id =~ s/\s.*$//;
        next;
    }
    $line =~ s/\s//g;
    $seq .= $line;
}
close ($in);
if ($id){
    $refseq{$id} = $seq;
}

my $sam_in;
if ($samfile =~ m/\.bam$/){
    open ($sam_in, "samtools view $samfile | ");
} else {
    die "ERROR: File must be a sorted and indexed bam file"
}

my %positions;
while (my $line = <$sam_in>){
    my @tmp = split("\t", $line);
    my ($rid, $flag, $ref, $start, $cigar) = @tmp[0,1,2,3,5];
    next if $cigar eq "*"; #unmapped reads not counted
    #my @cigars = split(/(?<=\D)/, $cigar);
    if ($cigar =~ m/[DI]/){ #read contains deletions or insertions
        my @cigars = split(/(?<=\D)/, $cigar);
        my $pos = $start;
        my $aleng = 0;
        foreach my $part (@cigars){
            if ($part =~ m/(\d+)([MID])/){
                my ($num, $op) = ($1, $2);
                #unless ($op eq "D"){
                #    $aleng += $num;
                #}
                if ($op =~ m/[ID]/){
                    my $opos = $pos;
                    $opos-- if $op eq "I";
                    push @{$positions{$ref}{$opos}}, ([$op, $rid, $flag, $num]);
                }
                next if $op eq "I"; #insertions relative to the reference don't change the position at all
                $pos += $num;
                #$pos = ($pos + $num) - 1;
            }
        }
    }
}
close ($sam_in);

my @dcorr;
my @icorr;
print "contig\tposition\tmut_leng\thp_base\thp_leng\ttotal_cov\tnum_del\tnum_ins\tdel+ins\tpct_del+ins\n";
foreach my $ref (sort keys %positions){
    foreach my $pos (sort {$a <=> $b} keys %{$positions{$ref}}){
        my ($dcount, $icount) = (0) x 2;
        my $mut_sum = 0;
        foreach (@{$positions{$ref}{$pos}}){
            my ($op, $rid, $flag, $num) = @{$_};
            $dcount++ if $op eq "D";
            $icount++ if $op eq "I";
            $mut_sum += $num;
            #print STDERR "$pos, $op, $rid, $flag\n";
        }
        my $sum = $dcount+$icount;
        next if $sum < 5;

        my $avg_mut = sprintf("%.0f", $mut_sum/$sum);

        ## check if this is a homopolymer position
        my $testpos = $pos;
        $testpos = $pos - 1 if $dcount > $icount;
        my $testseq = substr($refseq{$ref}, $testpos, 30);
        my $firstbase = substr($testseq, 0, 1);
        $testseq =~ m/^($firstbase+)/;
        my $hpleng = length($1);
        next unless ($hpleng >= $min_h_leng);
	
        my $tmpref = $ref;
        $tmpref =~ s/\|/\\|/g;
        $tmpref =~ s/\s.*//g;
        chomp(my $dline = `~/Applications/samtools-1.9/samtools depth -aa -r $tmpref:$pos-$pos $samfile`);
        (my $rdep) = $dline =~ m/(\d+)\D*$/;
        my $tcov = $dcount + $rdep;
        my $frac = sprintf("%.1f", 100 * ($sum / $tcov));

        next if $frac < $min_frac;

        print "$ref\t$pos\t$avg_mut\t$firstbase\t$hpleng\t$tcov\t$dcount\t$icount\t$sum\t$frac";
        if ($dcount > $icount){
            my $outpos = $pos;
            if ($avg_mut > 1){
                $outpos .= "-" . ($pos + ($avg_mut - 1));
            }
            print "\tD$outpos";
        } else {
            my $ostring = join("", ("$firstbase") x $avg_mut);
            print "\tI$pos$ostring";
        }
        print "\n";
    }
}
