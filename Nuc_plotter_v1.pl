#!/usr/bin/perl
# Written: Mark Robinson, April 2015
# Last updated: 
# USAGE: Nuc_plotter_vn.pl -options SGR_FILE PEAKS_FILE COORD_FILE ELIM_FILES(optional, comma seperated)
#
# Script to plot average dyad distributions for defined nucleosome categories for a given set of peak positions.
#
# Can also (optionally) eliminate sets of nucleosomes from plotting to observe influence. 
#
# Required Input: peak locations (.sgr), normalised dyad histograms (.sgr), and a gene coordinates file. 
#
# Optional Input: peak files (.sgr) of peaks to eliminate from analysis.  
#
# Output: a) files of the average normalised distribution of dyads around each specified nucleosome position (and standard error intervals)
#		  b) simplified version of above with just genic and intergenic averages
#		  c) if elimination files provided, will plot distribution of these changes accross nuc positions
#		  d) info file
#		  e) .sgr of peak categories determined to check in IGB
#
# NB: Numbering of output nucleosome positions will be set limits plus 3 averages in linear order.
#	  e.g. if all nucleosome limits set to 1, nucleosome order will be:
# 				0 = Upstream average
#				1 = -1 nucleosome
#				2 = +1 nucleosome
#				3 = Mid gene average
#				4 = TN
#				5 = TN+1
#				6 = Downstream average 
################################ Load Modules ################################

use strict;
use warnings;
use Math::Round qw(:all);
use Data::Dumper;
use List::Util qw(sum);
use List::MoreUtils 'pairwise';
use POSIX qw(strftime);
use Cwd qw();
use Getopt::Long;
use File::Basename;
use List::MoreUtils qw(any);

################################ Settings ################################

# Default settings
my $outdir = Cwd::cwd(); # current directory
my $bin = 10; # binning of sgr/peaks files
my $nuc_window = 50; # bp either side of peak to plot
my $overlap; # if want to remove overlapping peaks set min distance between peaks
my $extend = 1500; # bp up/down stream of gene to include in each gene (won't go into coding regions)
my $up_limit = 1; # max number of nucleosomes upstream from TSS to plot
my $start_limit = 3; # max number of nucleosomes downstream from TSS to plot
my $end_limit = 2; # max number of nucleosomes upstream from TTS to plot
my $down_limit = 1; # max number of nucleosomes downstream from TTS to plot
my $CI = 95; # % for confidence intervals (choose from 90, 95 or 99)
my $elim_list;

# getting user specified options
GetOptions ("bin=i" => \$bin,
            "window=i" => \$nuc_window,
            "overlap=i" => \$overlap,
            "add=i" => \$extend,
            "up=i" => \$up_limit,
            "start=i" => \$start_limit,
            "end=i" => \$end_limit,
            "down=i" => \$down_limit,
            "ci=i" => \$CI,
            "elim=s" => \$elim_list,
            "out=s" => \$outdir)
or die("Error in command line arguments\n");

# get files from command line
my $sgr_file = shift(@ARGV) or die "Inusufficient arguments supplied - usage: $0 -options SGR_FILE PEAKS_FILE COORD_FILE ELIM_FILES(optional)\n";
my $peaks_file = shift(@ARGV) or die "Inusufficient arguments supplied - usage: $0 -options SGR_FILE PEAKS_FILE COORD_FILE ELIM_FILES(optional)\n";
my $coord_file = shift(@ARGV) or die "Inusufficient arguments supplied - usage: $0 -options SGR_FILE PEAKS_FILE COORD_FILE ELIM_FILES(optional)\n";

my @elim_files = split(",",$elim_list);

################################ Main Program ################################

# Save settings as a string to print to info files later
my $settings = "Settings:\nBinning: $bin\nPlotting window: $nuc_window bp either side of peak\nExtension either side of gene: $extend bp\n";
$settings.="Number of peaks upstream plotted: $up_limit\nNumber of peaks at the 5' end plotted: $start_limit\nNumber of peaks at the 3' end plotted: $end_limit\n";
$settings.="Number of peaks downstream plotted: $down_limit\nConfidence intervals: $CI\%\n";
$settings.="Minimum peak spacing: $overlap bp\n" if ($overlap);

my $date = strftime "%d:%m:%Y", localtime;

print "\nStart:\t",`date`."\n";

################################ Gene Coordinates ################################

# check valid
my $crit_list = "909599";
die "$CI invalid arguement for CI, choose from 90,95 or 99\n" unless ($crit_list =~ $CI);

## convert CI to critical value (from z-table)
my $crit_val = 1.65 if ($CI == 90);
$crit_val = 1.96 if ($CI == 95);
$crit_val = 2.58 if ($CI == 99);

## Open coordinates file and add coordinates info to a hash and pos reference list
my (%coord_hash,%gene_ref);
    
open (my $coord_in, '<', "$coord_file") || die "Unable to open $coord_file: $!\n";

my ($coord_filename) = fileparse($coord_file);

print "Coordinates file used: $coord_filename \n";
    
while (<$coord_in>) {
    
    chomp;
    my @line = split('\t');
    my ($chrn,$start,$end,$id,$strand) = @line;
    
    # add gene info to coord hash
    $coord_hash{$id} = [$chrn,$start,$end,$strand];
    
    # add gene id to each bin covered by gene as ref
    ($start, $end) = (nearest($bin,$start),nearest($bin,$end));
    push(@{$gene_ref{$chrn}{$_}},$id) for (map {$_*$bin} $start/$bin..$end/$bin);
}


################################ Dyads ################################

my ($sgr_filename, $sgr_dir, $sgr_suffix) = fileparse($sgr_file,".sgr");
die "Error: $sgr_filename doesn't appear to be an .sgr file\n" unless ($sgr_suffix);

if ( $sgr_filename =~ /^(([A-Za-z\d]{1,10}_[A-Za-z\d]{1,10}_[A-Za-z\d]{1,10})_[A-Za-z\d-]{1,10})_.*/ ) {
    
    my $condition = $1;
    my $small_cond = $2;
    
    # open info file and print some info
    open (my $info_out, '>', "$outdir/Run_info_$condition\_$date.txt") || die "Unable to open file: $!\n";
    print ($info_out $settings);
    print ($info_out "\nInput dyad histogram file: $sgr_filename\n");
    
    my (%sgr_hash);
    
    open (my $sgr_in, '<', "$sgr_file") || die "Unable to open $sgr_file: $!\n";

    print "Found and processing: $sgr_filename$sgr_suffix\n";
    
    while (<$sgr_in>) {
        
        chomp;
        my @line = split('\t');
        my ($chrn,$pos,$freq) = @line;
        
        # add dyad freq to hash
        $sgr_hash{$chrn}{$pos} = $freq;
        
    }
    
    ################################ Elim ################################
    
    # find any elimination conditions matching current and store
    my ($elim_count,%elim_hash,@elim_cat,$elim_string) = (0);
    my @suffixes = (".txt",".sgr");
    
    for my $elim_file (@elim_files) {
        
        my ($elim_filename, $elim_dir, $elim_suffix) = fileparse($elim_file,@suffixes);
        
        die "Error: $elim_filename should be a .txt or .sgr file\n" unless ($elim_suffix);
        
        if ($elim_filename =~ /^$small_cond.*_([a-zA-Z]{1,10})_changes/) {
            
            my $change = $1;
            
            open (my $elim_in, '<', "$elim_file") || die "Unable to open $elim_file: $!\n";
            
            while (<$elim_in>) {
                
                chomp;
                my @line = split('\t');
                my ($chrn,$pos) = @line;
                $elim_hash{$chrn}{$pos}{$change}++; 
            }
            
            push(@elim_cat,$change);
            
            $elim_string.="\t$elim_filename\n";
            
            $elim_count++;
        }
    }
    
    print "$elim_count peak sets found for elimination from $condition\n";
    print ($info_out "\n$elim_count peak sets eliminated: ".join("\t",@elim_cat)."\n");
    print ($info_out "Elimination files used:\n$elim_string\n") if ($elim_string);
    
    ################################ Peaks ################################
    
    # will create HoA with following structure:
    #		%peaks_hash{$chrn}{$pos}[0] = [dyad values]
    #								[1] = [elimination categories]
    #								[2] = [nuc position categories]
    #								[3] = [termini relative orientation]
    #								[4] = [gene ids]
    #								[5] = [genic or intergenic]
    #								[6] = [orientation]
    
    # Note: Peaks are categorised twice: globally as either genic or intergenic and then in a gene-dependent manner as specific positions.
    #		This means intergenic peaks may be included multiple times (e.g. as TN+1 of one gene then -1 of following gene)
    
    # find correct peaks file
    my ($matched_peaks,%peaks_hash);
    
    my ($peaks_filename, $peaks_dir, $peaks_suffix) = fileparse($peaks_file,@suffixes);
    die "Error: $peaks_filename should be either a .txt or .sgr file\n" unless ($peaks_suffix);
        
    if ($peaks_filename =~ /$condition.*/) {
        
        print "Found and processing peaks from $peaks_filename\n";
        
        print ($info_out "Peaks file used: $peaks_filename\n");
        
        open (my $peaks_in, '<', "$peaks_file") || die "Unable to open $peaks_file: $!\n";
        
        PEAK: while (<$peaks_in>) {
            
            chomp;
            my @line = split('\t');
            my ($chrn,$pos) = @line;
            
            # skip header if present
            next PEAK if ($chrn =~ /chrn/i);
            
            # get dyad values in window
            for ( map {$_ * $bin} ($pos-$nuc_window)/$bin..($pos+$nuc_window)/$bin ) {
                
                push( @{$peaks_hash{$chrn}{$pos}[0]}, $sgr_hash{$chrn}{$_} || 0 );
                
            }
            
            # check if peak in elim hash
            for my $change (@elim_cat) {
                if (exists $elim_hash{$chrn}{$pos}{$change}) {
                    
                    # record any elimination category
                    push(@{$peaks_hash{$chrn}{$pos}[1]},$change);
                }
            }	
        }
    
    } else {
        
        die "Error: Peaks file doesn't appear to match .sgr - check file names\n";
    
    }
    
    ################################ Peak Mapping ################################
    
    print "Looking up dyad distributions\n";
    
    my (%genic_hash,%upstream_hash,%downstream_hash);
    
    # loop through all peaks and determine where they lie relative to gene start
    for my $chrn (sort keys %peaks_hash) {
        
        my ($prev_chrn,$prev_pos,$last_gene);
        
        POSITION: for my $pos (sort {$a <=> $b} keys %{$peaks_hash{$chrn}}) {
            
            my ($genic_flag);
            
            # skip first nuc of each chromosome
            if (! defined $prev_chrn || $chrn ne $prev_chrn ) {
                $prev_pos = $pos;
                $prev_chrn = $chrn;
                delete($peaks_hash{$chrn}{$pos});
                next POSITION;
            }
            
            # if bin is associated with a gene
            if (exists $gene_ref{$chrn}{$pos}) {
                                                
                # loop through every ID that bin is associated with
                for my $id (@{$gene_ref{$chrn}{$pos}}) {
                    
                    # set orientation relative to termini
                    my $mid_gene = ($coord_hash{$id}[2]-$coord_hash{$id}[1])/2; ### CHANGE INDEX WHEN COORD ALTERED
                    my $rel_pos = $pos-$coord_hash{$id}[1]; # position relative to start
                    my $orientation = $coord_hash{$id}[3];
                    
                    # mark current peak as genic and record orientation relative to termini
                    $peaks_hash{$chrn}{$pos}[5] = "genic";
                    $peaks_hash{$chrn}{$pos}[6] = "F" if ($rel_pos <= $mid_gene);
                    $peaks_hash{$chrn}{$pos}[6] = "R" if ($rel_pos > $mid_gene);
                    
                    # if want to elim overlapping & prev peak too close, skip this one
                    next POSITION if ($overlap && ($pos-$prev_pos) <= $overlap );
                    
                    # push/unshift onto array of nucleosomes mapping to this gene
                    push(@{$genic_hash{$id}},$pos) if ($orientation eq "F");
                    unshift(@{$genic_hash{$id}},$pos) if ($orientation eq "R");
                    
                    # flag peak as genic and keep track of id
                    $genic_flag++;
                    $last_gene = $id;
                    
                }
            }
                
            # for intergenic peaks...
            unless ($genic_flag) {
                
                my $found; 
                
                # if less than specified limit from previous gene, add to upstream/downstream list for this gene
                if ($last_gene && ($pos - $coord_hash{$last_gene}[2] < $extend) && $chrn eq $coord_hash{$last_gene}[0]) {  ### CHANGE INDEX WHEN COORD ALTERED
                        
                    # look up gene info
                    my $id = $last_gene;
                    my ($orientation) = ($coord_hash{$id}[3]); ### CHANGE INDEX WHEN COORD ALTERED
                    
                    # mark current peak as genic and record orientation
                    $peaks_hash{$chrn}{$pos}[5] = "intergenic";
                    $peaks_hash{$chrn}{$pos}[6] = "F";
                                        
                    # if want to elim overlapping & prev peak too close, skip this one
                    next POSITION if ($overlap && ($pos-$prev_pos) <= $overlap );
                    
                    # push/unshift onto array of nucleosomes mapping to this gene
                    push(@{$downstream_hash{$id}},$pos) if ($orientation eq "F");
                    unshift(@{$upstream_hash{$id}},$pos) if ($orientation eq "R");
                    
                    # flag as found
                    $found++;
                }
                
                # see if in range of next gene along
                WALK: for (my $index = $pos; $index <= ($pos+$extend); $index += $bin) {
                    
                    if (exists $gene_ref{$chrn}{$index}) {
                        
                        # record gene id
                        my $id = $gene_ref{$chrn}{$index}[0];
                        
                        # look up gene info
                        my ($orientation) = ($coord_hash{$id}[3]); ### CHANGE INDEX WHEN COORD ALTERED
                        
                        # mark current peak as genic and record orientation
                        $peaks_hash{$chrn}{$pos}[5] = "intergenic";
                        $peaks_hash{$chrn}{$pos}[6] = "R";
                        
                        # if want to elim overlapping & prev peak too close, skip this one
                        next POSITION if ($overlap && ($pos-$prev_pos) <= $overlap );
                        
                        # push/unshift onto array of nucleosomes mapping to this gene
                        push(@{$upstream_hash{$id}},$pos) if ($orientation eq "F");
                        unshift(@{$downstream_hash{$id}},$pos) if ($orientation eq "R");
                        
                        $found++;
                        last WALK;
                    }
                }
                
                # For any peaks not in range of a gene, just add as intergenic and give arbitrary orientation
                if (!$found) { 
                    $peaks_hash{$chrn}{$pos}[5] = "intergenic";
                    $peaks_hash{$chrn}{$pos}[6] = "F";
                }
                
            }
            $prev_pos = $pos;
            $prev_chrn = $chrn;
        }	
    }
    
    
    ################################ Peak Position Categories ################################	

    # peaks numbered 0..n depending on number plotted					
        
    print "Determining nucleosome position categories\n";
    
    for my $gene (keys %genic_hash) {
                        
        my $chrn = $coord_hash{$gene}[0];
        my $limit = round(@{$genic_hash{$gene}}/2); # set maximum amount of nucleosomes to count from each end
        my ($orientation) = ($coord_hash{$gene}[3]); ### CHANGE INDEX WHEN COORD ALTERED
        
        # gene start
        my ($index,$nuc_order) = (1,($up_limit+1)); # index number of +1 nucleosome will be equal to limit for upstream nuc+1
        until ( $index > $start_limit || $index > $limit) {
                        
            my $pos = shift(@{$genic_hash{$gene}});
            
            push(@{$peaks_hash{$chrn}{$pos}[2]},$nuc_order); 
            push(@{$peaks_hash{$chrn}{$pos}[4]},$gene); 
            
            push(@{$peaks_hash{$chrn}{$pos}[3]},"F") if ($orientation eq "F");
            push(@{$peaks_hash{$chrn}{$pos}[3]},"R") if ($orientation eq "R");
            
            $index++;
            $nuc_order++;
        }
            
        # gene end 
        my $pos = pop(@{$genic_hash{$gene}});
        
        ($index,$nuc_order) = (-1,($up_limit+$start_limit+1+$end_limit)); # have to add up desired number of nuc pos to get TN index number
        until ( $index < -$end_limit || $index < -$limit || ! defined $pos ) {
            
            push(@{$peaks_hash{$chrn}{$pos}[2]},$nuc_order);
            push(@{$peaks_hash{$chrn}{$pos}[4]},$gene);
            
            push(@{$peaks_hash{$chrn}{$pos}[3]},"F") if ($orientation eq "F");
            push(@{$peaks_hash{$chrn}{$pos}[3]},"R") if ($orientation eq "R");

            $index--;
            $nuc_order--;
            $pos = pop(@{$genic_hash{$gene}});
        
        }
            
        # if any peaks left add them to mid nuc average
        undef $pos;
        $pos = shift(@{$genic_hash{$gene}});
        
        $nuc_order = ($up_limit+$start_limit+1);
        until (! defined $pos) {
            
            push(@{$peaks_hash{$chrn}{$pos}[2]},$nuc_order);
            push(@{$peaks_hash{$chrn}{$pos}[4]},$gene);
            
            push(@{$peaks_hash{$chrn}{$pos}[3]},"F") if ($orientation eq "F");
            push(@{$peaks_hash{$chrn}{$pos}[3]},"R") if ($orientation eq "R");
            
            $pos = shift(@{$genic_hash{$gene}});
        }
    }
        
    # upstream:
    for my $gene (keys %upstream_hash) {
        
        my $chrn = $coord_hash{$gene}[0];
        my $nuc_order = $up_limit;
        my $orientation = $coord_hash{$gene}[3]; ### CHANGE INDEX WHEN COORD ALTERED
        
        my $pos = pop(@{$upstream_hash{$gene}});
        until (! defined $pos) {
            
            push(@{$peaks_hash{$chrn}{$pos}[2]},$nuc_order);
            push(@{$peaks_hash{$chrn}{$pos}[4]},$gene);
            
            push(@{$peaks_hash{$chrn}{$pos}[3]},"F") if ($orientation eq "F");
            push(@{$peaks_hash{$chrn}{$pos}[3]},"R") if ($orientation eq "R");
            
            # next pos
            $pos = pop(@{$upstream_hash{$gene}});
            $nuc_order-- unless ($nuc_order == 0); # will only decrease to 0 then keep adding to upstream average
        }
    }
        
    #downstream
    for my $gene (keys %downstream_hash) {
        
        my $chrn = $coord_hash{$gene}[0];
        my $nuc_order = ($up_limit+$start_limit+2+$end_limit);
        my $orientation = $coord_hash{$gene}[3]; ### CHANGE INDEX WHEN COORD ALTERED	
                    
        my $pos = shift(@{$downstream_hash{$gene}});
        until (! defined $pos) {
            
            push(@{$peaks_hash{$chrn}{$pos}[2]},$nuc_order);
            push(@{$peaks_hash{$chrn}{$pos}[4]},$gene);
            
            push(@{$peaks_hash{$chrn}{$pos}[3]},"F") if ($orientation eq "F");
            push(@{$peaks_hash{$chrn}{$pos}[3]},"R") if ($orientation eq "R");
            
            # next pos
            $pos = shift(@{$downstream_hash{$gene}});
            $nuc_order++ unless ($nuc_order >= ($up_limit+$start_limit+2+$end_limit+$down_limit));
        }
    }
    # pppfffftttt! that was confusing, moving on. 
    
    # print categorisations to check
    my $temp_file = "$condition\_nucleosome_position_categories.sgr";
    open (my $temp_out, '>', "$outdir/$temp_file") || die "Unable to open $temp_file: $!\n";
    
    for my $chrn (sort keys %peaks_hash) {
        for my $pos (sort {$a<=>$b} keys %{$peaks_hash{$chrn}}) {
            my $temp_pos = $pos;
            for (0..$#{$peaks_hash{$chrn}{$pos}[2]}) {
                print ($temp_out "$chrn\t$temp_pos\t$peaks_hash{$chrn}{$pos}[2][$_]\n");
                $temp_pos++;
            }
        }
    }
    
    ################################ Summing Dyads ################################	
    
    print "Summing dyads distributions around peaks\n";
    my (%plot_hash);
    
    # get average dyad distributions for each nuc set (total + elim)
    for my $chrn (keys %peaks_hash) {
        
        for my $pos (keys %{$peaks_hash{$chrn}}) {
            
            # for genic/intergenic average
            if ($peaks_hash{$chrn}{$pos}[6] eq "F") {
                
                my $cat = $peaks_hash{$chrn}{$pos}[5]; # get category (genic/intergenic)
                $plot_hash{"Total"}{$cat}{"dyads"}[$_] += $peaks_hash{$chrn}{$pos}[0][$_] for (0..$#{$peaks_hash{$chrn}{$pos}[0]}); # add dyads to sum
                $plot_hash{"Total"}{$cat}{"count"} ++;
                                    
            } elsif ($peaks_hash{$chrn}{$pos}[6] eq "R") {
                
                my $cat = $peaks_hash{$chrn}{$pos}[5]; # get category (genic/intergenic)
                my @yarra = reverse(@{$peaks_hash{$chrn}{$pos}[0]}); # reverse the dyad array
                $plot_hash{"Total"}{$cat}{"dyads"}[$_] += $yarra[$_] for (0..$#yarra); # add dyads to sum
                $plot_hash{"Total"}{$cat}{"count"} ++;
            
            }
            
            # for specific positions
            for my $nuc_position (@{$peaks_hash{$chrn}{$pos}[2]}) {
                
                if ($peaks_hash{$chrn}{$pos}[3] eq "F") {
                    
                    $plot_hash{"Total"}{$nuc_position}{"dyads"}[$_] += $peaks_hash{$chrn}{$pos}[0][$_] for (0..$#{$peaks_hash{$chrn}{$pos}[0]}); # add dyads to sum
                    $plot_hash{"Total"}{$nuc_position}{"count"} ++;
                                        
                } elsif ($peaks_hash{$chrn}{$pos}[6] eq "R") {
                    
                    my @yarra = reverse(@{$peaks_hash{$chrn}{$pos}[0]}); # reverse the dyad array
                    $plot_hash{"Total"}{$nuc_position}{"dyads"}[$_] += $yarra[$_] for (0..$#yarra); # add dyads to sum
                    $plot_hash{"Total"}{$nuc_position}{"count"} ++;
                
                }
            }
            
            # repeat whole thing again skipping any nuc in elimination sets
            ELIM: for my $elim (@elim_cat) {
                
                # check if current elim category in list
                if ( any { lc($_) eq $elim } @{$peaks_hash{$chrn}{$pos}[1]} ) {
                    
                    next ELIM;
                    
                } else {
                    
                    # for genic/intergenic average
                    if ($peaks_hash{$chrn}{$pos}[6] eq "F") {
                    
                        my $cat = $peaks_hash{$chrn}{$pos}[5]; # get category (genic/intergenic)
                        $plot_hash{$elim}{$cat}{"dyads"}[$_] += $peaks_hash{$chrn}{$pos}[0][$_] for (0..$#{$peaks_hash{$chrn}{$pos}[0]}); # add dyads to sum
                        $plot_hash{$elim}{$cat}{"count"} ++;
                                            
                    } elsif ($peaks_hash{$chrn}{$pos}[6] eq "R") {
                        
                        my $cat = $peaks_hash{$chrn}{$pos}[5]; # get category (genic/intergenic)
                        my @yarra = reverse(@{$peaks_hash{$chrn}{$pos}[0]}); # reverse the dyad array
                        $plot_hash{$elim}{$cat}{"dyads"}[$_] += $yarra[$_] for (0..$#yarra); # add dyads to sum
                        $plot_hash{$elim}{$cat}{"count"} ++;
                        
                    }
                    
                    # for specific positions
                    for my $nuc_position (@{$peaks_hash{$chrn}{$pos}[2]}) {
                        
                        if ($peaks_hash{$chrn}{$pos}[3] eq "F") {
                            
                            $plot_hash{$elim}{$nuc_position}{"dyads"}[$_] += $peaks_hash{$chrn}{$pos}[0][$_] for (0..$#{$peaks_hash{$chrn}{$pos}[0]}); 
                            $plot_hash{$elim}{$nuc_position}{"count"} ++;
                                                
                        } elsif ($peaks_hash{$chrn}{$pos}[6] eq "R") {
                            
                            my @yarra = reverse(@{$peaks_hash{$chrn}{$pos}[0]});
                            $plot_hash{$elim}{$nuc_position}{"dyads"}[$_] += $yarra[$_] for (0..$#yarra);
                            $plot_hash{$elim}{$nuc_position}{"count"} ++;
                        
                        }
                    }
                }
            }
        }
    }
    
    #print Dumper(\%plot_hash);
    
    ################################ Error Calculation ###############################
    
    print "Calculating errors\n";
    
    # loop through each set and calc mean
    unshift (@elim_cat,"Total");
    for my $elim (@elim_cat) {
        
        for my $nuc_pos (keys %{$plot_hash{$elim}}) {
            
            for (0..$#{$plot_hash{$elim}{$nuc_pos}{"dyads"}}) {
                
                $plot_hash{$elim}{$nuc_pos}{"dyads"}[$_] = $plot_hash{$elim}{$nuc_pos}{"dyads"}[$_]/$plot_hash{$elim}{$nuc_pos}{"count"};
            }
        }
    }
    
    # loop through each peak position and calc deviation (this seems long winded - theres probably a less stupid way to do this!)
    for my $chrn (keys %peaks_hash) {
        for my $pos (keys %{$peaks_hash{$chrn}}) {
            
            # copy dyad values
            my @dyads = @{$peaks_hash{$chrn}{$pos}[0]};
            
            # genic/intergenic
            my $nuc_pos = $peaks_hash{$chrn}{$pos}[5];
            @dyads = reverse (@dyads) if ($peaks_hash{$chrn}{$pos}[6] eq "R");
            
            CAT: for my $elim (@elim_cat) {
                
                # skip if current change category flagged for current peak
                next CAT if ($elim ~~ @{$peaks_hash{$chrn}{$pos}[1]});
                
                # otherwise calc squared deviation from the mean and add to appropriate plot set					
                $plot_hash{$elim}{$nuc_pos}{"error"}[$_] += ($plot_hash{$elim}{$nuc_pos}{"dyads"}[$_] - $dyads[$_])**2 for (0..$#dyads);
            
            }
            
            # repeat for specific nuc positions
            for (0..$#{$peaks_hash{$chrn}{$pos}[2]}) {
                my $nuc_pos = $peaks_hash{$chrn}{$pos}[2][$_];
                my @dyads = reverse (@dyads) if ($peaks_hash{$chrn}{$pos}[3][$_] eq "R");
                
                CAT2: for my $elim (@elim_cat) {
                    
                    # skip if current change category flagged for current peak
                    next CAT2 if ( any {lc($_) eq $elim} @{$peaks_hash{$chrn}{$pos}[1]} );
                    
                    # otherwise calc squared deviation from the mean and add to appropriate plot set					
                    $plot_hash{$elim}{$nuc_pos}{"error"}[$_] += ($plot_hash{$elim}{$nuc_pos}{"dyads"}[$_] - $dyads[$_])**2 for (0..$#dyads);
                
                }
            }
        }
    }
    
    # then loop through different plotting sets and calc standard deviation 
    for my $elim (keys %plot_hash) {
        for my $nuc_pos (keys %{$plot_hash{$elim}}) {
            for (0..$#{$plot_hash{$elim}{$nuc_pos}{"error"}}) {
                
                # divide sum of squares by n-1 and take square root
                $plot_hash{$elim}{$nuc_pos}{"error"}[$_] = ($plot_hash{$elim}{$nuc_pos}{"error"}[$_]/($plot_hash{$elim}{$nuc_pos}{"count"}-1)) ** 0.5;
                
                # divide SD by square root of sample size (=SE)
                $plot_hash{$elim}{$nuc_pos}{"error"}[$_] = ( $plot_hash{$elim}{$nuc_pos}{"error"}[$_]/( $plot_hash{$elim}{$nuc_pos}{"count"} ** 0.5 ) );
                
                # multiply by critical value to get CI
                $plot_hash{$elim}{$nuc_pos}{"error"}[$_] = nearest( 0.00001, $plot_hash{$elim}{$nuc_pos}{"error"}[$_] * $crit_val );
            }
        }
    }
    
    ################################ Printing ################################
    
    print "Printing output files\n";
    
    ## print average nuc
    
    # Open & header
    my $av_nuc_file = $condition."_average_nuc_profiles_$date.txt";
    my $nuc_file = $condition."_nuc_position_profiles_$date.txt";
    open (my $av_nuc_out, '>', "$outdir/$av_nuc_file") || die "Unable to open $av_nuc_file: $!\n";
    open (my $nuc_out, '>', "$outdir/$nuc_file") || die "Unable to open $nuc_file: $!\n";
    
    my $header = "Nuc_pos\tBin\t".join("\t",@elim_cat)."\n";
    my @av_cats = ("genic","intergenic");
    
    ## Print average profiles
    
    print ($av_nuc_out "Dyad Distribution\n");
    print ($av_nuc_out $header);
    
    # Dyad values
    for my $cat (@av_cats) {
        
        print ($av_nuc_out "$cat:");
        
        for (my $pos = -$nuc_window, my $index = 0; $pos <= $nuc_window; $pos += $bin, $index++) {
                
            print ($av_nuc_out "\t$pos");
            
            # print dyads for each set
            for my $elim (@elim_cat) {
                
                my $val = $plot_hash{$elim}{$cat}{"dyads"}[$index];
                print ($av_nuc_out "\t$val");
                
            }
            
            print ($av_nuc_out "\n");
        }
        
        print ($av_nuc_out "\n");
    }
    
    # Errors
    print ($av_nuc_out "\n$CI Confidence Intervals:\n");
    print ($av_nuc_out $header);
    for my $cat (@av_cats) {
        
        print ($av_nuc_out "$cat:");
        
        for (my $pos = -$nuc_window, my $index = 0; $pos <= $nuc_window; $pos += $bin, $index++) {
                
            print ($av_nuc_out "\t$pos");
            
            for my $set (@elim_cat) {
                
                my $val = $plot_hash{$set}{$cat}{"error"}[$index];
                print ($av_nuc_out "\t$val");
                
            }
            
            print ($av_nuc_out "\n");
        }
        
        print ($av_nuc_out "\n");
    }
    
    # print number of nucleosomes in each set
    print ($info_out "\nNucleosome numbers for genic/intergenic averages:\n");
    $header = "Nuc_pos\t".join("\t",@elim_cat)."\n";
    print ($info_out $header);
    for my $cat (@av_cats) {
        print ($info_out "$cat:");
        for my $elim (@elim_cat) {
            print ($info_out "\t$plot_hash{$elim}{$cat}{count}");
        }
        print ($info_out "\n");
    }
    
    # new header for +/- 95% CI
    $header = "Nuc_pos\tBin\t";
    my @header;
    for my $elim (@elim_cat) {
        push(@header,("$elim+CI","$elim","$elim-CI"));
    }
    $header = $header.join("\t",@header)."\n";
    
    # Dyad values +/- 95% CI (just for easier plotting)
    print ($av_nuc_out "\nDyad values +/- $CI\% CI:\n");
    print ($av_nuc_out $header);
    for my $cat (@av_cats) {
        
        print ($av_nuc_out "$cat:");
        
        for (my $pos = -$nuc_window, my $index = 0; $pos <= $nuc_window; $pos += $bin, $index++) {
                
            print ($av_nuc_out "\t$pos");
            
            for my $set (@elim_cat) {
                
                my $val = $plot_hash{$set}{$cat}{"dyads"}[$index];
                my $lower = $val - $plot_hash{$set}{$cat}{"error"}[$index];
                my $upper = $val + $plot_hash{$set}{$cat}{"error"}[$index];
                print ($av_nuc_out "\t$upper\t$val\t$lower");
                
            }
            
            print ($av_nuc_out "\n");
        }
        
        print ($av_nuc_out "\n");
    }
    
    
    ## Print specific nucleosome positions
    
    # determine nucleosome position reference and print
    my ($ref,@nuc_pos_ref,$string) = (-1);
    for (1..$up_limit) {
        unshift(@nuc_pos_ref,$ref) if ($up_limit > 0);
        $ref--;
    }
    for (1..$start_limit) {
        push(@nuc_pos_ref,"+$_") if ($start_limit > 0);
    }
    $ref=-$end_limit+1;
    for (1..$end_limit) {
        $string = "TN$ref";
        $string = "TN" if ($ref == 0);
        push(@nuc_pos_ref,"$string");
        $ref++;
    } 
    for (1..$down_limit) {
        $string = "TN+$_";
        push(@nuc_pos_ref,"$string") if ($down_limit > 0);
    }
    unshift (@nuc_pos_ref, "Upstream average");
    push (@nuc_pos_ref, "Downstream average");
    $nuc_pos_ref[$_] = "$_ = $nuc_pos_ref[$_]" for (0..$#nuc_pos_ref);
    $string = join(", ",@nuc_pos_ref);
    print ($info_out	"\nNucleosome order for specified positions: $string\n\n");
    
    # get numeric positions
    my @nuc_pos;
    for my $nuc_pos (keys %{$plot_hash{"Total"}}) {
        push(@nuc_pos,$nuc_pos) if ($nuc_pos =~ /\d+/);
    }
    
    print ($nuc_out $header);
    
    for my $nuc_pos (sort {$a <=> $b} @nuc_pos) {
            
        print ($nuc_out "$nuc_pos");
        
        for (my $pos = -$nuc_window, my $index = 0; $pos <= $nuc_window; $pos += $bin, $index++) {
        
            print ($nuc_out "\t$pos");
            
            for my $set (@elim_cat) {
                
                my $val = $plot_hash{$set}{$nuc_pos}{"dyads"}[$index];
                my $lower = $val - $plot_hash{$set}{$nuc_pos}{"error"}[$index];
                my $upper = $val + $plot_hash{$set}{$nuc_pos}{"error"}[$index];
                print ($nuc_out "\t$upper\t$val\t$lower");						
            }
            print ($nuc_out "\n");
        }
        
        print ($nuc_out "\n");
        
    }	
    
    # print number of nucleosomes in each set
    print ($info_out "Nucleosome numbers for specified positions:\n");
    $header = "Nuc_pos\t".join("\t",@elim_cat)."\n";
    print ($info_out $header);
    for my $cat (sort {$a <=> $b} @nuc_pos) {
        print ($info_out "$cat");
        for my $elim (@elim_cat) {
            my $count = $plot_hash{$elim}{$cat}{"count"};
            print ($info_out "\t$count");
        }
        print ($info_out "\n");
    }
}

print "End ",`date`,"\n";