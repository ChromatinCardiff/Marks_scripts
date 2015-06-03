#!/usr/bin/perl
# Written: Mark Robinson, April 2015
# Last updated:
# USAGE: perl --options Nuc_dist_vn.pl SGR_FILE PEAKS_FILE
#
# Script to measure how tightly nucleosomes are positioned accross the population.
#
# Takes cumulative sum of dyads in a given window around peak and calculates gradient of 
# central region as a measure of "tightness". Allows a minimum distance between peaks to be set
# if you want to exclude overlapping peaks from analysis.
#
# Also will eliminate any peaks with a total read number below a specified limit (to deal with
# false peaks from KDE smoothing).
#
# Input: a dyad histogram .sgr (un-smoothed) and a peaks file (.sgr) for each condition.
#
# Output: adds gradient and overlap status (if wanted) as extra columns to peaks file.
#
################################ Load Modules ################################

use strict;
use warnings;
use Math::Round qw(:all);
use Data::Dumper;
use Cwd qw();
use Getopt::Long;
use File::Basename;

################################ Settings ################################

# default settings
my $outdir = Cwd::cwd(); # current directory
my $bin = 10; # binning of sgr and peaks files
my $peak_window = 5; # bins either side of peak to sum dyads
my $grad_window = 2; # which bins to calculate gradient over
my $exclusion_window = 60; # minimum distance (bp) between peaks (set as '' if don't want to exclude overlapping)
my $min_reads = 10; # minimum number of reads within window to count as peak

# getting user specified options
GetOptions ("bin=i" => \$bin,
            "pwind=i" => \$peak_window,
            "gwind=i" => \$grad_window,
            "ewind=i" => \$exclusion_window,
            "min=i" => \$min_reads,
            "out=s" => \$outdir)
or die("Error in command line arguments\n");

# get files from command line
my $sgr_file = shift(@ARGV) or die "Inusufficient arguments supplied - usage: $0 -options SGR_FILE PEAKS_FILE\n";
my $peaks_file = shift(@ARGV) or die "Inusufficient arguments supplied - usage: $0 -options SGR_FILE PEAKS_FILE\n";

################################ Main Program ################################

print "Start: ".`date`."\n";
print "Settings:\n";
print "Bins either side of peaks for summing dyads: $peak_window\n";
print "Window for calculating gradient from cumulative sum: $grad_window\n";

## Read all sgr file values into hash

my (%sgr_hash);

my @suffixes = (".sgr",".txt");
my ($filename, $dir, $suffix) = fileparse($sgr_file,@suffixes);
die "Error: $filename doesn't appear to be an .sgr or .txt file!\n" unless ($suffix);

if ($filename =~ /^([A-Za-z\d]{1,10}_[A-Za-z\d]{1,10}_[A-Za-z\d]{1,10}_[A-Za-z\d\-]{1,10})/) {
    
    my $condition = $1;
    
    open (my $sgr_in, '<', "$sgr_file") || die "Unable to open $sgr_file: $!\n";

    print "Found and processing: $filename \n";
    
    while (<$sgr_in>) {
        chomp;
        my @line = split('\t');
        my $chrn = $line[0];
        my $pos = $line[1];
        my $val = $line[2];
        $sgr_hash{$condition}{$chrn}{$pos} = $val;
    }
}

my ($peaks_filename, $peaks_dir, $peaks_suffix) = fileparse($peaks_file,".sgr");
die "Error: $peaks_filename doesn't appear to be an .sgr file!\n" unless ($peaks_suffix);

if ($peaks_filename =~ /^([A-Za-z\d]{1,10}_[A-Za-z\d]{1,10}_[A-Za-z\d]{1,10}_[A-Za-z\d\-]{1,10}).*(t\d+)/) {
    
    my $condition = $1;
    my $peak_thresh = $2;
    
    open (my $peaks_in, '<', "$peaks_file") || die "Unable to open $peaks_file: $!\n";
    
    # open output files
    my $out_filename = "$condition\_$peak_thresh\_dist.txt";
    open (my $peaks_out, '>', "$outdir/$out_filename") || die "Unable to create file $out_filename: $!\n";
    $out_filename = "$condition\_$peak_thresh\_e$exclusion_window\_overlapping.sgr";
    open (my $elim_out, '>', "$outdir/$out_filename") || die "Unable to create file $out_filename: $!\n";
    
    print "Found and processing: $peaks_filename \n";
    
    # If exclusion window defined loop through once to get hash of overlapping peaks
    my ($prev_chrn,$prev_pos,%excl_hash);
    if ($exclusion_window) {
        while (<$peaks_in>) {
            
            chomp;
            my @line = split("\t");
            my $chrn = $line[0];
            my $pos = $line[1];
            
            if ($prev_chrn && $prev_pos) {
                unless ($chrn ne $prev_chrn) {
                    # if distance less than limit add both peaks to exclusion hash
                    if ($pos-$prev_pos < $exclusion_window) {
                        $excl_hash{$condition}{$chrn}{$pos}++;
                        $excl_hash{$condition}{$prev_chrn}{$prev_pos}++;
                        print ($elim_out "$chrn\t$pos\t$exclusion_window\n");
                    }
                }
            } 
            
            $prev_chrn = $chrn;
            $prev_pos = $pos;
        }
        seek($peaks_in,0,0);
    }
    
    my ($line_count,$peaks_count,$inc_count,$skipped_count,$excl_count,$ends_count) = (0,0,0,0,0,0);
    
    PEAK: while (<$peaks_in>) {
        
        $line_count++;
                
        chomp;
        my @line = split('\t');
        my $chrn = $line[0];
        my $pos = $line[1];
        my $first_pos = $pos-($bin*$peak_window);
        my $last_pos = $pos+($bin*$peak_window);
        
        # check for header and either add on a gradient column or write new one
        if ($line_count == 1) {
            if ($line[0] eq "ChrN") {
                print ($peaks_out join("\t",@line)."\tGrad\n");
                next PEAK;
            } else {
                print ($peaks_out "ChrN\tPos\tHeight\tGrad\n");
            }
        }
        $peaks_count++;
        
        # if peak in exclusion hash, skip
        if (exists $excl_hash{$condition}{$chrn}{$pos}) {
            $excl_count++;
            next PEAK;
        }
        
        # cumulatively sum values accross bin window
        my ($sum,@sum_array);
        if ( defined $sgr_hash{$condition}{$chrn}{$first_pos} && defined $sgr_hash{$condition}{$chrn}{$last_pos} ) {
            for ( my $i = $first_pos, my $j = 0; $i <= $last_pos; $i += $bin ,$j++ ) {
                $sum += $sgr_hash{$condition}{$chrn}{$i};
                $sum_array[$j] = $sum;
            }
        } else {
            $ends_count++;
            next PEAK;
        }
        
        # skip peaks with low read numbers
        if ($sum < $min_reads) {
            $skipped_count++;
            next PEAK;
        } else {
            $inc_count++; 
        }
        
        # convert to frequency
        my @cfd = map { $_ / $sum } @sum_array;
        
        # calc gradient over 5 bins around peak
        my $grad = nearest(0.00001,($cfd[$peak_window+$grad_window]-$cfd[$peak_window-$grad_window])/(($grad_window*2)+1)); # $peak_window index will equal dyad position here
        
        # print peak back out with gradient added
        print ($peaks_out join("\t",@line)."\t$grad\n");
        
    }
    
    print "$peaks_count peaks processed for $peaks_filename\n";
    print "$inc_count peaks output\n";
    print "$skipped_count peaks skipped due to low read depth (<$min_reads)\n";
    print "$ends_count peaks skipped due to proximity to chromosome ends\n";
    
    if ($exclusion_window) {
        print "$excl_count overlapping peaks excluded (<$exclusion_window bp apart)\n";
    }
}

print "End: ".`date`."\n";