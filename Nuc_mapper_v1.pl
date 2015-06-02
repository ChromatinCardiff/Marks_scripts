#!/usr/bin/perl
# Written: Mark Robinson, April 2015
# Last updated:
# USAGE: Nuc_mapper_vn.pl
#
# Maps reads from SAM to genome, producing an .sgr histogram of dyad positions and also a histogram
# of read lengths and alignment numbers.
#		
# Input: SAM file of aligned reads from MNase-seq 
#
# Output: Nucleosome map style .sgr histogram of dyads and histogram of read lengths (.txt)
#
# NB: can also be used to output SAMparser style files (e.g. for use with Nick's KDE smoothing script)
# 
################################ Load Modules ################################

use strict;
use warnings;
use Math::Round qw(:all);
use Data::Dumper;
use POSIX;
use Getopt::Long;
use Cwd qw();

################################ Settings ################################

# default option settings
my $bin = 10; # binning to use for sgr file
my $min_size = 120; # minimum read length to map
my $max_size = 180; # maximum read length to map
my $SAMparser = ''; # if SAMparser format files required set as 1, else set as ''
my $PE = ''; # if data is paired end sequenced - set as 1, else set as ''
my $outdir = Cwd::cwd();

# getting user specified options
GetOptions ("bin=i" => \$bin,
            "min=i" => \$min_size,
            "max=i"   => \$max_size,
            "SAMparser"  => \$SAMparser,
            "PE" => \$PE,
            "out=s" => \$outdir)
or die("Error in command line arguments\n");

# getting command line file names
my $SAMfile = shift(@ARGV) or die "Inusufficient arguments supplied - usage: $0 -options SAMfile.sam\n";

################################ Main Program ################################

print "\nstart:\t", `date`."\n";
print "Mapping MNase protected fragments between $min_size and $max_size basepairs to the genome\n";
print "Note: SAMparser output option turned on\n" if $SAMparser;
print "Running in paired end mode\n" if $PE;
print "Binning: $bin bp\n";

if ( $SAMfile !~ /^\.+/ && $SAMfile =~ /\/([A-Za-z\d]{1,10}_[A-Za-z\d]{1,10}_[A-Za-z\d]{1,10}).*\.sam/) {
    
    my $condition = $1;
    
    print "Found and processing $SAMfile \n";
    
    open(my $sam_in, '<', "$SAMfile") || die "Unable to open $SAMfile: $!";
    
    my (%samhash,$read_count,$mapped_reads,$unaligned_count,$size_read_count,%size_hash,%chrn_lengths,%SAMparser_hash);
    
    LINE: while(<$sam_in>) {
        
        chomp;

        # Initialise chromosome sizes from header
        if ($_ =~ /^@/) {
            my @header = split('\t');
            if ($header[1] =~ /^SN:(chr[\dFMRB]+)/) {
                my $chrn = $1;
                if ($header[2] =~ /^LN:(\d+)/) {
                    my $length = ceil($1/$bin);
                    $chrn_lengths{$chrn} = $length;
                }
            }
            next LINE;
        }
        
        # count all reads and split line
        $read_count++;
        my @line = split('\t');
        
        # Check for chromosome alignment and skip if unaligned
        my $chrn;
        if ($line[2] =~ /chr[\dMRFBIVX]+/) { # only includes 
            $chrn = $line[2];
        } elsif ($line[2] =~ /(NC_\d+)/) { # handle yeast refseq references in stead of chromosome names using a sub
            my $refseq = $1;
            $chrn = &refseqname($refseq);
        } elsif ($line[2] =~ /\*/) {
            $unaligned_count++;
            next LINE;
        }
        
        if (!defined $chrn) {
            print "Couldn't find chromosome name for read at line: $.\n";
            next LINE;
        }
        
        # Determine fragment size bin according to user specified bins
        my $read = $line[10];
        my $frag_size = nearest($bin,length($read));
        
        # change this to distance between paired reads if paired 
        $frag_size = nearest($bin,$line[8]) if $PE; 
        
        # increment count of size and mapped reads
        $size_hash{$frag_size}++;
        $mapped_reads++;
        
        # Check fragment length is within window
        if ($frag_size >= $min_size && $frag_size <= $max_size) {
            
            # Determine dyad position
            my $dyad_pos = nearest($bin,($line[3] + ($frag_size * 0.5)));
            
            # Increment the count for caluclated fragment size bin at dyad position
            $samhash{$chrn}{$dyad_pos}++;
            
            # Increment count of reads within window
            $size_read_count++; 
            
            # if SAMparser output wanted, print out
            if ($SAMparser) {
                
                push(@{$SAMparser_hash{$chrn}{"Start"}},$line[3]);
                push(@{$SAMparser_hash{$chrn}{"Frag"}},$frag_size);
                push(@{$SAMparser_hash{$chrn}{"Pos"}},$dyad_pos);
            }
        }
    }
    
    print "\nTotal reads input: $read_count\n";
    print "Reads mapping to the genome: $mapped_reads\n";
    print "Unaligned reads: $unaligned_count\n";
    print "Reads within specified size range: $size_read_count\n\n";			
    
    # open output files
    my $outfile = "$condition\_read_distribution_$bin.txt";
    open (my $sizeout, '>', "$outdir/$outfile") || die "Unable to open output file: $outfile\n";
    
    $outfile = "$condition\_$min_size\-$max_size"."bp_histogram_$bin.sgr";
    open (my $sgrout, '>', "$outdir/$outfile") || die "Unable to open output file: $outfile\n";
    
    # 3 bin smoothing average print sgr
    for my $chrn (sort keys %chrn_lengths) {
        for (0..$chrn_lengths{$chrn}) {
            # smooth sgr sum with 3 bin average and print
            my $pos = $_*$bin;
            my $sum_val = nearest(0.0001,( (($samhash{$chrn}{$pos-$bin}||0)+($samhash{$chrn}{$pos}||0)+($samhash{$chrn}{$pos+$bin}||0))/3 ));
            print ($sgrout "$chrn\t$pos\t$sum_val\n");
        }
    }
    
    # print size distribution histogram
    print ($sizeout "Reads mapped from $SAMfile\n");
    print ($sizeout "Total reads: $read_count\nReads aligned:$mapped_reads\n\n");
    print ($sizeout "Fragment size\tRaw reads\tFreq\n");
    for my $size (sort {$a <=> $b} keys %size_hash) {
        my $frag_count = $size_hash{$size};
        my $freq = nearest(0.00001,$frag_count/$mapped_reads);
        print ($sizeout "$size\t$frag_count\t$freq\n");
    }
    
    # print SAMparser files
    if ($SAMparser) {
        for my $chrn (sort keys %SAMparser_hash) {
            
            my $filename = "$chrn\_1_$min_size\-$max_size"."bp_$condition\_SAMparser.txt";
            open (my $sam_out, '>', "$outdir/$filename") || die "Unable to open $filename: $!\n";
            
            for (0..$#{$SAMparser_hash{$chrn}}) {
                
                my ($start,$frag,$pos) = ($SAMparser_hash{$chrn}{"Start"},$SAMparser_hash{$chrn}{"Frag"},$SAMparser_hash{$chrn}{"Pos"});
                print ($sam_out "$chrn\t$start\t$frag\t$pos\n");
            
            }
        }
    }

}
print "End: ".`date`."\n";