# Pilon Tools  
## pilon_repeater.pl
Automatically performs multiple successive rounds of Pilon correction on an assembly file until no further changes are produced.  
**NOTE 1:** You must change the paths to bwa, samtools, and pilon in the first few lines of code to match the locations in your system.  
**NOTE 2:** The script is hardcoded to request 16 Gb of RAM for pilon. If you want to change this, change the `-Xmx16G` call on line 66 
to the amount of RAM you'd like to request.  

## pilon\_results_check.pl  
Looks for any leftover homopolymer errors after pilon correction of an assembly.   
**NOTE:** Must change paths to samtools on lines 55 and 121 to the path to samtools on your system. Minimum samtools version is v1.9.

## fasta\_indel_inserter.pl  
Uses the output of pilon\_results_check.pl to correct residual homopolymer errors in the corrected assembly.  
