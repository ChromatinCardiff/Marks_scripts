#!/usr/bin/perl
# Written: Nick Kent,  Aug 2010
# Last updated: Mark Robinson, May 2015
# USAGE:- perl Nuc_caller_v1.pl

# Basically just Nick's dicty_peakmarker.pl, main alterations are outputting a 3 bin sum 
# for peak height value and preventing dual peak calling within 1 bin of each other.
#
# This script takes an .sgr file as an input, and calls peak centre/summit 
# bins above a single, but scalable, threshold.
# It then lists these bin positions with a y-axis value proportional to the scaled 
# read#frequency.
#
# The scaling value can be chosen to reflect differences in read depth between two
# experiments - either based on total depth or a SiteWriter-derived local depth.

################################ Load Modules ################################

use strict;
use warnings;
use Math::Round;
use Cwd qw();
use Getopt::Long;
use File::Basename;
use List::MoreUtils qw(any);

################################ Settings ################################

# Default settings
my $outdir = Cwd::cwd(); # current directory
my $thresh = 10; # Peak calling threshold
my $max_thresh = 200; # Maximum threshold over which peaks are ignored
my $scale_factor = 1.0; # If want to scale for differences in read depth
my $dual_bins = 2; # minimum number of bins between peaks

# getting user specified options
GetOptions ("thresh=i" => \$thresh,
            "max=i" => \$max_thresh,
            "scale=f" => \$scale_factor,
            "dual=i" => \$dual_bins,
            "out=s" => \$outdir)
or die("Error in command line arguments\n");

# get files from command line
my $file = shift(@ARGV) or die "Inusufficient arguments supplied - usage: $0 -options FILENAME\n";

################################ Main Program ################################

# print some useful shizz
print "Start: ".`date`."\n";
print "Settings:\n";
print "Scale factor set to: $scale_factor \n";
print "Maximum threshold set to: $max_thresh\n";
print "Minimum threshold set to: $thresh \n";

# check file name
my ($filename, $dir, $suffix) = fileparse($file,".sgr");

die "Error: $filename doesn't appear to be an .sgr file!\n" unless ($suffix);

# define outfile name from infile name
my $outfile = "$filename\_peak_t$thresh\.sgr";

# print out some useful info
print "Processing .sgr file: $filename\n";

open ( my $in, '<', "$file" ) || die "Unable to open $file: $!";

# define three new arrays to store required values from infile
my (%bins,%freq);

# loop through sgr and store positions and values for each chromosome
while(<$in>){
    
    chomp;
    my @line = split('\t');
    my $chrn = $line[0];
    push(@{$bins{$chrn}},$line[1]);
    push(@{$freq{$chrn}},$line[2]);
    
}

# open output file
open ( my $out, '>', "$outdir/$outfile" ) || die "Unable to open $outfile: $!";

## Peak Calling ##

my ($peak_count);

for my $chrn (sort keys %freq) {
    
    my ($count,@bin_tracker);
    
    LABEL: for (0..$#{$freq{$chrn}}) {
        
        $count = $_;
        
        # skip any positions within 2 bins of chromosome ends or over max threshold
        if (!defined $freq{$chrn}[$count-2] || !defined $freq{$chrn}[$count+2] || $freq{$chrn}[$count] >= $max_thresh) {
            push(@bin_tracker,"null");
            next LABEL;
        }
        
        # check if bin value greater than surrounding bins and above threshold
        if (($freq{$chrn}[$count]>$freq{$chrn}[$count-2]) &&
            ($freq{$chrn}[$count]>$freq{$chrn}[$count-1]) && 
            ($freq{$chrn}[$count]>$freq{$chrn}[$count+1]) &&
            ($freq{$chrn}[$count]>$freq{$chrn}[$count+2]) &&
            ($freq{$chrn}[$count]*$scale_factor>=$thresh) &&
            ($freq{$chrn}[$count]*$scale_factor<=$max_thresh)){
            
            # check neighbouring bins haven't been called
            if ( $count == 0 ) {
                
                push(@bin_tracker,"null");
                next LABEL;
            }
            
            if ($dual_bins >= 1) {
                if ( any { $bin_tracker[$count-$_] ne "null" } 1..$dual_bins ) {
 
                    push(@bin_tracker,"null");
                    next LABEL;
                }
            }
            
            # print to peaks file
            push(@bin_tracker,$bins{$chrn}[$count]);
            print($out "$chrn\t$bins{$chrn}[$count]\t".round(($freq{$chrn}[$count-1]+$freq{$chrn}[$count]+$freq{$chrn}[$count+1])*$scale_factor)."\n");
            $peak_count++;
            
        } else {
            push(@bin_tracker,"null");
        }
    }
}
print "Peaks called for $filename = $peak_count \n";

print "End: ".`date`."\n";