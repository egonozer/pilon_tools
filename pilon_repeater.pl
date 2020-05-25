#!/usr/bin/perl

use strict;
use warnings;

my $bwa = "/usr/bin/bwa";
my $samtools = "/usr/bin/samtools";
my $pilon = "/home/hauserlab/Applications/pilon/pilon-1.23.jar";

#bwa index pilon3.fasta ; samtools faidx pilon3.fasta; bwa mem -t 40 pilon3.fasta ~/lacie/reads/C_auris_reads/CA01_1.fastq.gz ~/lacie/reads/C_auris_reads/CA01_2.fastq.gz | samtools sort > aln4.bam ; samtools index aln4.bam ; java -Xmx16G -jar ~/Applications/pilon/pilon-1.23.jar --genome pilon3.fasta --frags aln4.bam --output pilon4 --fix all --mindepth 0.5 --changes --verbose --threads 40

my $usage = "
pilon_repeater.pl

Performs multiple rounds of pilon on an assembly until the number of changes
introduced by pilon falls below a given threshold. Will also stop if the
number of changes performed is static.

Required:
  -a    assembly file, in fasta format
  -1    path to forward Illumina read file, in fastq format (can be gzipped)
  -2    path to reverse Illumina read file, in fastq format (can be gzipped)
  
Options:
  -o    output file prefix
        (default: 'pilon')
  -m    maximum number of changes performed by pilon before stopping
        (default: 0)
  -d    minimum depth
        (default: 0.1)
  -t    threads
        (default: 40)
        
";

use Getopt::Std;
our ($opt_a, $opt_1, $opt_2, $opt_o, $opt_m, $opt_d, $opt_t);
getopts('a:1:2:o:m:');

die $usage unless ($opt_a and $opt_1 and $opt_2);
my ($afile, $rfile1, $rfile2) = ($opt_a, $opt_1, $opt_2);
my $pref    = $opt_o ? $opt_o : "pilon";
my $maxch   = $opt_m ? defined $opt_m : 0;
my $mindep  = $opt_d ? $opt_d : 0.1;
my $threads = $opt_t ? $opt_t : 40;

my $iteration = 0;
my $done = 0;
my $genome = $afile;
my $last_numchanges = 0;
while (!$done) {
    if (-e "aln$iteration.bam"){
        unlink("aln$iteration.bam");
        unlink("aln$iteration.bam.bai");
    }
    $iteration++;
    my $status = system("bwa index $genome >/dev/null");
    die "ERROR: bwa exited with status $status\n" if $status;
    $status = system("samtools faidx $genome");
    die "ERROR: samtools faidx exited with status $status\n" if $status;
    print STDERR "Running bwa...\n";
    $status = system("bwa mem -t $threads $genome $rfile1 $rfile2 2>/dev/null | samtools sort > aln$iteration.bam 2>/dev/null");
    die "ERROR: bwa exited with status $status\n" if $status;
    $status = system("samtools index aln$iteration.bam");
    die "ERROR: samtools index exited with status $status\n" if $status;
    open (my $lines, "java -Xmx16G -jar ~/Applications/pilon/pilon-1.23.jar --genome $genome --frags aln$iteration.bam --output $pref$iteration --fix all --mindepth $mindep --changes --verbose --threads $threads 2>&1 | ");
    my $symb = "|";
    while(<$lines>){
        $symb = spinner($symb);
        print STDERR "\rRunning pilon $symb";
    }
    print STDERR "\n";
    if (-e "$pref$iteration.changes"){
        unlink("$genome.amb");
        unlink("$genome.ann");
        unlink("$genome.bwt");
        unlink("$genome.fai");
        unlink("$genome.pac");
        unlink("$genome.sa");
        chomp(my $numchanges = `wc -l $pref$iteration.changes`);
        $numchanges =~ s/^(\d+).*/$1/;
        print STDERR "Iteration $iteration: $numchanges changes\n";
        $done = 1 if $numchanges <= $maxch;
        if ($iteration > 1){
            $done = 1 if $numchanges == $last_numchanges;
        }
        $last_numchanges = $numchanges;
        $genome = "$pref$iteration.fasta";
    } else {
        die "ERROR: Could not find file '$pref$iteration.changes'\n";
    }
}
print STDERR "Done\n";


sub spinner {
    my $symb = shift;
    if ($symb eq "|"){
        $symb = "\\";
    } elsif ($symb eq "\\"){
        $symb = "-";
    } elsif ($symb eq "-"){
        $symb = "/";
    } elsif ($symb eq "/") {
        $symb = "|";
    }
    return($symb);
}