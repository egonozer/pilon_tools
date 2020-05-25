#!/usr/bin/perl

use strict;
use warnings;

$|++;

my $usage = "
fasta_indel_inserter.pl

Will take a fasta sequence and add or delete bases according to the list given.

Required:
  -f    sequence file to be modified
  -c    list of changes, one per line

Changes to be made should take the following formats:

First column should be the name of the contig to change, second column the
change to be made. Separate columns with a tab. One change per line.

All positions are 1-based (i.e. the first base of a sequence is position 1)

Deletions should be denoted by a D and the position(s) to remove.
Examples:
\"contig_1  D15\": remove the 15th base on contig_1
\"contig_3  D200-203\": remove bases 200, 201, 202, and 203 on contig_3

Insertions should be denoted by an I, the position of the base before the
insertion, and the bases to insert.
Examples:
\"contig_1  I44A\": insert the base A after the base at position 44 on contig_1
\"contig_4  I19GGG\": insert the bases GGG after the base at position 19 on contig _4

Single base mutations should be denoted by an M, the position of the base, and
the replacement base.
Exampe:
\"contig_2  M87A\": replace the base at position 87 with an A


";

use Getopt::Std;
our ($opt_f, $opt_c);
getopts('f:c:');

die $usage unless ($opt_f and $opt_c);
my $file = $opt_f;
my $cfile = $opt_c;

open (my $cin, "<$cfile") or die "ERROR: Can't open $cfile: $!\n";
my %hash;
my $incount = 0;
while (my $line = <$cin>){
    chomp $line;
    my ($ctg, $mut) = split("\t", $line);
    if ($mut =~ m/^([IDM])(\d+-*\d*)(\D*)/){
        my ($type, $pstring) = ($1, $2);
        my $repl = $3 if $3;
        if ($type eq "D"){
            my ($start, $stop) = split("-", $pstring);
            $stop = $start unless $stop;
            for my $i ($start .. $stop){
                push @{$hash{$ctg}{$i}}, "#";
                $incount++;
            }
        } elsif ($type eq "I") {
            unless ($repl){
                print STDERR "WARNING:\"$line\" is missing a replacement sequence. Skipping\n";
                next;
            }
            push @{$hash{$ctg}{$pstring}}, "$repl";
            $incount++;
        } else {
            unless ($repl){
                print STDERR "WARNING:\"$line\" is missing a replacement base. Skipping\n";
                next;
            }
            push @{$hash{$ctg}{$pstring}}, "&$repl";
            $incount++;
        }
    }
}

open (my $in, "<$file") or die "ERROR: Can't open $file: $!\n";
my ($id, $seq);
my $outcount = 0;
while (my $line = <$in>){
    chomp $line;
    $line =~ s/\r//g;
    if ($line =~ m/^>/){
        if ($id){
            print ">$id\n";
            if ($hash{$id}){
                my @array = split(//, $seq);
                my $count = 0;
                foreach my $base (@array){
                    $count++;
                    if ($hash{$id}{$count}){
                        my @changes = sort @{$hash{$id}{$count}};
                        foreach my $change (@changes){
                            $outcount++;
                            if ($change eq "#"){
                                $base = "";
                            } else {
                                if ($change =~ m/^&/){
                                    $base = substr($change, 1, 1);
                                } else {
                                    $base .= $change;
                                }
                            }
                        }
                    }
                    print "$base";
                }
                print "\n";
            } else {
                print "$seq\n";
            }
        }
        $id = substr($line, 1);
        $seq = "";
        next;
    }
    $line =~ s/\s//g;
    $seq .= $line;
}
close ($in);
if ($id){
    print ">$id\n";
    if ($hash{$id}){
        my @array = split(//, $seq);
        my $count = 0;
        foreach my $base (@array){
            $count++;
            if ($hash{$id}{$count}){
                my @changes = sort @{$hash{$id}{$count}};
                foreach my $change (@changes){
                    $outcount++;
                    if ($change eq "#"){
                        $base = "";
                    } else {
                        if ($change =~ m/^&/){
                            $base = substr($change, 1, 1);
                        } else {
                            $base .= $change;
                        }
                    }
                }
            }
            print "$base";
        }
        print "\n";
    } else {
        print "$seq\n";
    }
}

print STDERR "changes to be made: $incount\nchanges made: $outcount\n";
