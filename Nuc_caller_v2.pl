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
my $grad_thresh; # minimum gradient below which peaks are treated as disorganised regions
my $window = 5; # bins either side of peak to include for gradient calc
my $help; # whether to print usage info

# getting user specified options
GetOptions ("thresh=i" => \$thresh,
            "max=i" => \$max_thresh,
            "scale=f" => \$scale_factor,
            "dual=i" => \$dual_bins,
            "help" => \$help,
            "window=i" => \$window,
            "grad=f" => \$grad_thresh,
            "out=s" => \$outdir)
or die("Error in command line arguments\n");

my $help_string = "\nUsage: $0 -options path/file.sgr\nOptional settings:\n";
$help_string.="--thresh|-t = minimum dyad frequency for calling peaks (default = 10)\n";
$help_string.="--max|-m = maximum dyad frequency of peaks (default = 200)\n";
$help_string.="--scale|-s = scaling factor (default = 1.0)\n";
$help_string.="--dual|-d = minimum number of bins between called peaks (default = 2)\n";
$help_string.="--out|-o = output directory (default = current working directory)\n";
$help_string.="--grad|-g = minimum gradient threshold below which region is treated as a disorganised region not a peak (default = off) \n";
$help_string.="--window|-w = bins over which to sum dyads for gradient calculation (default = 5) \n";
$help_string.="--help|-h = print help info (but then I guess you already knew that!)\n\n";

# print help info in flagged
if ($help) {
    print $help_string;
    exit;
}

# get files from command line
my $file = shift(@ARGV) or die "Inusufficient arguments supplied\n$help_string";

################################ Main Program ################################

# print some useful shizz
print "Start: ".`date`."\n";
print "Settings:\n";
print "Scale factor set to: $scale_factor \n";
print "Maximum threshold set to: $max_thresh\n";
print "Minimum threshold set to: $thresh \n";
print "Minimum peak spacing set to: $dual_bins\n";
if ($grad_thresh) {
    print "Gradient measuring activated - gradient threshold set to $grad_thresh\n";
}

# check file name
my ($filename, $dir, $suffix) = fileparse($file,".sgr");

die "Error: $filename doesn't appear to be an .sgr file!\n" unless ($suffix);

# define outfile name from infile name
my $outfile = "$filename\_peak_t$thresh\.sgr";
$outfile = "$filename\_peak_g$grad_thresh\_t$thresh\.sgr" if ($grad_thresh);

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
    
    my ($prev_peak,$count,@bin_tracker) = (0);
    
    LABEL: for (0..$#{$freq{$chrn}}) {
        
        $count = $_;
        
        # skip any positions within 2 bins of chromosome ends or over max threshold
        if (!defined $freq{$chrn}[$count-2] || !defined $freq{$chrn}[$count+2] || $freq{$chrn}[$count] >= $max_thresh) {
            push(@bin_tracker,0);
            next LABEL;
        }
        
        # check if bin value greater than surrounding bins and above threshold
        if (($freq{$chrn}[$count]>$freq{$chrn}[$count-2]) &&
            ($freq{$chrn}[$count]>=$freq{$chrn}[$count-1]) &&
            ($freq{$chrn}[$count]>=$freq{$chrn}[$count+1]) &&
            ($freq{$chrn}[$count]>$freq{$chrn}[$count+2]) &&
            ($freq{$chrn}[$count]*$scale_factor>=$thresh) &&
            ($freq{$chrn}[$count]*$scale_factor<=$max_thresh) &&
            (!defined $bin_tracker[$count])) {
            
            ## Dealing with overlapping peaks
                
            # check neighbouring bins upstream haven't already been called
            if ( $count == 0 ) {
                
                push(@bin_tracker,0);
                next LABEL;
            }
            
            if ($dual_bins >= 1) {
                
                if ( any { $bin_tracker[$count-$_] > 0 } 1..$dual_bins ) {
                    
                    push(@bin_tracker,0);
                    next LABEL;
                }
            
            }
                
            ## Dealing with disorganised regions
            
            if ($grad_thresh) {
                
                # sum dyads to get cumulative distribution
                my ($index,$sum,@grad_array) = (0,0);
                    
                for ($count-$window..$count+$window) {
                    
                    $sum += $freq{$chrn}[$_];
                    $grad_array[$index] = $sum;
                    $index++;
                
                }
                
                # normalise cumulative array
                my @cfd = map {$_ / $sum } @grad_array;
                    
                # measure gradient of right side of distribution
                my $grad = nearest( 0.00001, ($cfd[$window+2]-$cfd[$window-2]) / 5);
                    
                # if gradient below limit and less than 3x threshold, treat as disorganized region
                if ( $grad < $grad_thresh && $freq{$chrn}[$count] <= $thresh*3 ) {
                    
                    # back track to find upstream edge of region
                    my $value = $freq{$chrn}[$count];
                    my $start_pos = $count;
                    
                    until ( $value < ($thresh/2) || $start_pos <= $prev_peak+4 ) {
                        
                        $start_pos--;
                        $value = $freq{$chrn}[$start_pos];
                    }
                    
                    # find where coverage drops below the threshold again
                    $value = $freq{$chrn}[$count];
                    my $end_pos = $count;
                    
                    until ( $value < ($thresh/2) ) {
                        
                        $end_pos++;
                        $value = $freq{$chrn}[$end_pos];
                    }
                    
                    # find centre of region
                    my $peak_pos = $start_pos + (($end_pos-$start_pos)/2);
                    
                    # add appropriate values to bin tracker
                    push (@bin_tracker,0) for ($count..$end_pos);
                    $bin_tracker[$peak_pos] = $freq{$chrn}[$peak_pos];
                    
                    # print peak
                    #print ($out "$chrn\t$bins{$chrn}[$peak_pos]\t".$freq{$chrn}[$peak_pos]*$scale_factor."\n");
                    print ($out "$chrn\t$bins{$chrn}[$peak_pos]\t2\n");
                    $peak_count++;
                    
                    # skip to next position
                    $prev_peak = $peak_pos;
                    next LABEL;
                    
                }
            }
                
            
            # print to peaks file
            push(@bin_tracker,$freq{$chrn}[$count]);
            #print($out "$chrn\t$bins{$chrn}[$count]\t".$freq{$chrn}[$count]*$scale_factor."\n");
            print($out "$chrn\t$bins{$chrn}[$count]\t1\n");
            $peak_count++;
            $prev_peak = $count;
            
        } elsif (!defined $bin_tracker[$count]) {
            
            push(@bin_tracker,0);
            
        }
    }
}
print "Peaks called for $filename = $peak_count \n";

print "End: ".`date`."\n";