#!/usr/bin/perl
# Written: Mark Robinson, April 2015
# Last updated:
# USAGE: Nuc_changes_vn.pl
#
# Takes a peaks file with added nucleosome parameter data (height/dist/occ/size)
# for 2 or more conditions and measures change in each, outputting lists of altered positions.
#
# Input: Reference peaks file(s) and test peaks file(s) with matching parameter measures.
#
# Output: One .sgr style file of change in each input parameter (bp change for pos and size,
#		  log2 fold change for gradient and height) per test file/ref file combo
#
#		  Also outputs an info file for each combo with general info on parameter changes.
#
################################ Load Modules ################################

use strict;
use warnings;
use Math::Round qw(:all);
use Data::Dumper;
use Cwd qw();
use Getopt::Long;
use List::MoreUtils qw(any);

################################ Settings ################################

# Default settings
my $bin = 10; # binning of sgr and peaks files (bp)
my $common_window = 40; # threshold below which peaks are considered commonly positioned (bp)
my $pos_thresh = 10; # threshold over which peaks are considered re-positioned (bp)
my $size_thresh = 20; # threshold over which peaks are considered to have altered size (bp)
my $grad_thresh = 1.1; # threshold over which peaks considered re-distributed (fold change)
my $height_thresh = 1.3; # threshold over which peaks considered to have altered occupancy (fold change)
my $outdir = Cwd::cwd();
my $help;

# getting user specified options
GetOptions ("bin=i" => \$bin,
            "common=i" => \$common_window,
            "pos=i"   => \$pos_thresh,
            "size=i"  => \$size_thresh,
            "grad=f" => \$grad_thresh,
            "height=f" => \$height_thresh,
            "help" => \$help,
            "out=s" => \$outdir)
or die("Error in command line arguments\n");

my $help_string = "\nUsage: $0 -options ref_file test_file\nOptional settings:\n";
$help_string .= "--bin|-b = binning (default = 10)\n";
$help_string .= "--common|-c = bp window either side of peak to look for common peaks between conditions (default = 40)\n";
$help_string .= "--pos|-p = bp position change threshold for calling peaks with altered positioning (default = 10)\n";
$help_string .= "--size|-s = bp size change threshold for calling peaks with altered fragment size (default = 20)\n";
$help_string .= "--grad|-g = fold gradient change threshold for calling peaks with altered distribution (default = 1.1)\n";
$help_string .= "--height|-h = fold height change threshold for calling peaks with altered height (default = 1.3)\n";
$help_string .= "--out|-o = output directory (default = current working directory)\n";

if ($help) {
    print $help_string;
    exit;
}

# get files from command line
my $ref_file = shift(@ARGV) or die "Inusufficient arguments supplied $help_string";
my $test_file = shift(@ARGV) or die "Inusufficient arguments supplied $help_string";

################################ Main Program ################################

print "Start: ".`date`."\n";
print "Settings:\n";
print "Binning: $bin\n";
print "Window for common peaks: $common_window\n";
print "Threshold for position change: $pos_thresh\n";
print "Threshold for size change: $size_thresh\n";
print "Threshold for gradient change: $grad_thresh\n";
print "Threshold for height change: $height_thresh\n\n";

# convert fc thresholds to log2 fc
my $log_grad_thresh = log($grad_thresh)/log(2);
my $log_height_thresh = log($height_thresh)/log(2); 

## Read all reference peak paramaters into HoA

my (%ref_hash);

if ($ref_file !~ /^\.+/ && $ref_file =~ /\/([A-Za-z\d]{1,10}_[A-Za-z\d]{1,10}_[A-Za-z\d]{1,10})_.*\.txt/) {
		
    my ($ref_cond,@header,$ref_count) = ($1);
    
    open (my $ref_in, '<', "$ref_file") || die "Unable to open $ref_file: $!\n";

    print "Found and processing: $ref_file \n";
    
    while (<$ref_in>) {
                    
        chomp;
        my @line = split('\t');
        
        # store header
        if ($line[0] eq "ChrN") {
            for (0..$#line) {
                push(@header,$line[$_]);
            }
        } else {
            # add each peaks parameters to HoA (take keys from header)
            my $chrn = $line[0];
            my $pos = $line[1];
            
            for (2..$#line) {
                $ref_hash{$chrn}{$pos}{$header[$_]} = $line[$_];
            }
            
            $ref_count++;
        }
    }
    
    ## Loop through all test files and compare against ref


    if ($test_file !~ /^\.+/ && $test_file =~ /\/([A-Za-z\d]{1,10}_[A-Za-z\d]{1,10}_[A-Za-z\d]{1,10})_.*\.txt/) {
                            
        my ($cond,@params,$peaks_count,$common_count,%changes) = ($1);
        
        open (my $test_in, '<', "$test_file") || die "Unable to open $test_file: $!\n";
        
        print "Found and processing: $test_file \n";
        
        # open output files
        my $outfile = "$cond\_vs_$ref_cond\_pos_changes.sgr";
        open (my $pos_out, '>', "$outdir/$outfile") || die "Unable to open $outfile: $!\n";
        $outfile = "$cond\_vs_$ref_cond\_height_changes.sgr";
        open (my $height_out, '>', "$outdir/$outfile") || die "Unable to open $outfile: $!\n";
        $outfile = "$cond\_vs_$ref_cond\_unique_changes.sgr";
        open (my $uniq_out, '>', "$outdir/$outfile") || die "Unable to open $outfile: $!\n";
        
        $outfile = "$ref_cond\_vs_$cond\_pos_changes.sgr";
        open (my $ref_pos_out, '>', "$outdir/$outfile") || die "Unable to open $outfile: $!\n";
        $outfile = "$ref_cond\_vs_$cond\_height_changes.sgr";
        open (my $ref_height_out, '>', "$outdir/$outfile") || die "Unable to open $outfile: $!\n";
        
        $outfile = "$cond\_vs_$ref_cond\_info.txt";
        open (my $info_out, '>', "$outdir/$outfile") || die "Unable to open $outfile: $!\n";
        
        # if gradient parameter present in reference set, open files for gradient and occupancy
        my ($grad_out,$ref_grad_out,$occ_out,$ref_occ_out);
        if ( any { lc($_) eq 'grad' } @header  ) {
            
            $outfile = "$cond\_vs_$ref_cond\_grad_changes.sgr";
            open ($grad_out, '>', "$outdir/$outfile") || die "Unable to open $outfile: $!\n";
            $outfile = "$cond\_vs_$ref_cond\_occ_changes.sgr";
            open ($occ_out, '>', "$outdir/$outfile") || die "Unable to open $outfile: $!\n";
            
            $outfile = "$ref_cond\_vs_$cond\_grad_changes.sgr";
            open ($ref_grad_out, '>', "$outdir/$outfile") || die "Unable to open $outfile: $!\n";
            $outfile = "$ref_cond\_vs_$cond\_occ_changes.sgr";
            open ($ref_occ_out, '>', "$outdir/$outfile") || die "Unable to open $outfile: $!\n";
        }
        
        # if size paramaters present in reference set, open files for size change
        my ($size_out,$ref_size_out);
        if ( any { lc($_) eq 'size' } @header ) {
            
            $outfile = "$cond\_vs_$ref_cond\_size_changes.sgr";
            open ( $size_out, '>', "$outdir/$outfile") || die "Unable to open $outfile: $!\n";
            
            $outfile = "$ref_cond\_vs_$cond\_size_changes.sgr";
            open ( $ref_size_out, '>', "$outdir/$outfile") || die "Unable to open $outfile: $!\n";
        }
        
        PEAK: while (<$test_in>) {
            
            chomp;
            my @line = split('\t');
            
            # store order of parameters from header
            if ($line[0] eq "ChrN") {
                
                for (0..$#line) {
                    push(@params,$line[$_]);
                }
            
            
            # for each peak, check parameter, compare against ref & output
            } else {
                
                $peaks_count++;
                
                my $chrn = $line[0];
                my $pos = $line[1];
                
                ## Find common peaks ##
                
                # check if peak(s) within specified window for common peaks
                my (@pos_matches);
                
                for (my $ref_pos = $pos-$common_window; $ref_pos <= $pos+$common_window; $ref_pos += $bin) {
                    
                    if (exists $ref_hash{$chrn}{$ref_pos}) {
                        
                        push(@pos_matches,$ref_pos);
                    }
                }
                
                # if no matches found, print to unique peaks file and skip
                if (! @pos_matches) {
                    
                    print (	$uniq_out "$chrn\t$pos\t1\n");
                    
                    next PEAK;							
                }
                
                $common_count++;
                
                # find ref peak closest to current peak
                my ($diff,$ref_pos) = ($common_window);
                for (0..$#pos_matches) {
                    
                    my $temp_diff = $pos-$pos_matches[$_];
                    if ( abs($temp_diff) <= $diff ) {
                        
                        $diff = $temp_diff;
                        $ref_pos = $pos_matches[$_];
                    }
                }
                
                # mark reference peak as mapped (so can print unique peaks in WT later)
                $ref_hash{$chrn}{$ref_pos}{"mapped"}++;
                
                ## Positioning ##
                
                # track all position changes
                push(@{$changes{"Pos"}{"Values"}},$diff);
                $changes{"Pos"}{"Total"} += $diff;
                
                # print posioning diff if over threshold
                if ( abs($diff) > $pos_thresh ) {
                    
                    print ($pos_out "$chrn\t$pos\t$diff\n");
                    print ($ref_pos_out "$chrn\t$pos\t-$diff\n");
                }
                
                ## Height ##
                
                # calc log2 fold change in height
                my $height_fc = $line[2] / $ref_hash{$chrn}{$ref_pos}{"Height"};
                $height_fc = log($height_fc)/log(2);
                
                # store all changes to work out stats later 
                push(@{$changes{"Height"}{"Values"}},$height_fc);
                $changes{"Height"}{"Total"} += $height_fc;
                
                # if over threshold, print
                if ( abs($height_fc) > $log_height_thresh ) {
                    
                    print ($height_out "$chrn\t$pos\t$height_fc\n");
                    print ($ref_height_out "$chrn\t$pos\t-$height_fc\n");
                }
                
                # loop through other parameters 
                for (3..$#params) {
                    
                    ## Gradient ## 
                    if ( $params[$_] =~ /grad/i ) {
                        
                         # calc fold change
                         my $grad_fc = $line[$_] / $ref_hash{$chrn}{$ref_pos}{"Grad"};
                         $grad_fc = log($grad_fc)/log(2);
                         
                         # store all changes to work out stats later 
                         push(@{$changes{"Grad"}{"Values"}},$grad_fc);
                         $changes{"Grad"}{"Total"} += $grad_fc;
                         
                         # if over threshold, print out
                         if ( abs($grad_fc) > $log_grad_thresh ) {
                            
                            print ($grad_out "$chrn\t$pos\t$grad_fc\n");
                            print ($ref_grad_out "$chrn\t$pos\t-$grad_fc\n");
                         
                         } elsif ( abs($height_fc) > $log_height_thresh ) { ## v2 edit:  if no "significant" gradient change, but is height change, print to occ sgr
                    
                            print ($occ_out "$chrn\t$pos\t$height_fc\n");
                            print ($ref_occ_out "$chrn\t$pos\t-$height_fc\n");
                         }
                    }
                    
                    ## Size ##
                    if ( $params[$_] =~ /size/i ) {
                        
                        
                        # calc absolute change
                        my $size_change = $line[$_] - $ref_hash{$chrn}{$ref_pos}{"Size"};
                        
                        # store all changes to work out stats later
                        push(@{$changes{"Size"}{"Values"}},$size_change);
                        $changes{"Size"}{"Total"} += $size_change;
                        
                        # if over threshold, print out
                        if ( abs($size_change) > $size_thresh ) {
                            
                            print ($size_out "$chrn\t$pos\t$size_change\n");
                            print ($ref_size_out "$chrn\t$pos\t$size_change\n");
                            
                        }
                    }
                    
                }
            }
        }
            
        # print some info to file
        print ($info_out "Test condition: $cond\n");
        print ($info_out "Reference condition: $ref_cond\n");
        print ($info_out "$ref_count peaks found in $ref_file\n");
        print ($info_out "$peaks_count peaks found in $test_file\n");
        print ($info_out "$common_count peaks in common between conditions\n");
        print ($info_out "\nSettings:\nBinning: $bin\nWindow for common peaks: $common_window\n");
        print ($info_out "Threshold for position change: $pos_thresh\nThreshold for size change: $size_thresh\n");
        print ($info_out "Threshold for gradient change: $grad_thresh\nThreshold for height change: $height_thresh\n\n");
        print ($info_out "\nParameter\tAverage change\tStandard Deviation\n");
        
        
        ## Calculate SD for each paramter
        
        for my $param (keys %changes) {
            
            my $average = nearest(0.000001,$changes{$param}{"Total"}/$common_count); # average
            my $sos;
            $sos += ($average-$_)**2 foreach (@{$changes{$param}{"Values"}}); # sum of squares
            my $sd = nearest(0.000001,($sos/($common_count-1))** 0.5); # standard deviation
            
            # print info
            print ($info_out "$param\t$average\t$sd\n");
        }
        
        # print unique peaks from WT
        $outfile = "$ref_cond\_vs_$cond\_unique_changes.sgr";
        open (my $ref_uniq_out, '>', "$outdir/$outfile") || die "Unable to open $outfile: $!\n";
        
        for my $chrn (sort keys %ref_hash) {
            
            for my $pos (sort {$a <=> $b} keys %{$ref_hash{$chrn}}) {
                
                unless (exists $ref_hash{$chrn}{$pos}{"mapped"}) {
                    
                    print ($ref_uniq_out "$chrn\t$pos\t1\n");
                }
            }
        }
    }
}
print "End: ".`date`."\n";