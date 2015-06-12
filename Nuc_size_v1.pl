#!/usr/bin/perl
# Written: Mark Robinson, April 2015
# Last updated:
# USAGE: Nuc_size_vn.pl
#
# Sums fragment sizes in given size range within a window either side of identified peaks. 
#
# Input files required: peak locations (.sgr) & an aligned SAM file of reads from MNase-seq. 
#
# Output: Average 3D & 2D peak profile (.txt), adds peak fragment size to peaks file (.sgr/.txt),
#		  + c3 file of frag distribution around each peak for clustering.
#		  + standard dyad histogram (with 3-bin smooth) for checking mapping
# 
################################ Load Modules ################################

use strict;
use warnings;
use Math::Round qw(:all);
use Data::Dumper;
use List::Util qw(sum min max);
use Cwd qw();
use Getopt::Long;
use File::Basename;

################################ Settings ################################

# Default settings
my $outdir = Cwd::cwd(); # current directory
my $pos_window = 50; # bp either side of the peak position to include
my $min_frag = 100; # smallest fragment size to include in mapping
my $max_frag = 200; # largest fragment size to include in mapping
my $frag_bin = 5; # binning for fragment sizes
my $pos_bin = 10; # binning for positions (should match peaks file)
my $av_bin = ''; # optional higher resolution binning for 2d average, set to '' if want to save memory
my $PE = ''; # whether data is paired ends or single (set to 1 for PE, or '' for SE)
my $help;

# getting user specified options
GetOptions ("wind=i" => \$pos_window,
            "max=i" => \$max_frag,
            "min=i" => \$min_frag,
            "fbin=i" => \$frag_bin,
            "pbin=i" => \$pos_bin,
            "abin=i" => \$av_bin,
            "PE" => \$PE,
            "help" => \$help,
            "out=s" => \$outdir)
or die("Error in command line arguments\n");

my $help_string = "\nUsage: $0 -options SAM_file.sam peaks_file.txt\nOptional settings:\n";
$help_string.="--wind|-w = bp either side of peak to count as peak region (default = 50)\n";
$help_string.="--min|-mi = minimum fragment size to include (default = 100)\n";
$help_string.="--max|-ma = maximum fragment size to include (default = 200)\n";
$help_string.="--fbin|-f = binning for fragment sizes (default = 5)\n";
$help_string.="--pbin|-p = binning for genomic position (default = 10)\n";
$help_string.="--abin|-a = optional binning for the average 2D/3D profile so can set seperate binning for c3 and average if desired (default = off)\n";
$help_string.="--PE|-p = flag paired end reads (default = off)\n";
$help_string.="--out|-o = output directory (default = current working directory)\n\n";

print $help_string if ($help);

# get files from command line
my $SAM_file = shift(@ARGV) or die "Inusufficient arguments supplied\n$help_string";
my $peaks_file = shift(@ARGV) or die "Inusufficient arguments supplied\n$help_string\n";

################################ Main Program ################################

print "start:\t", `date`."\n";
print "Analysing fragment size variability in fragment window $min_frag to $max_frag bp\n";
print "Analysing for sum of fragments +/- ".$pos_window/$pos_bin." bins either side of called peaks in $pos_bin bp bins\n";

my $condition;

my ($filename, $dir, $suffix) = fileparse($SAM_file,".sam");
die "Error: $filename doesn't appear to be a SAM file\n" unless ($suffix);


if ( $filename =~ /^([A-Za-z\d]{1,10}_[A-Za-z\d]{1,10}_[A-Za-z\d]{1,10})/ ) {
    
    $condition = $1;
    
    print "Found and processing $filename \n";
    
    open(my $sam_in, '<', "$SAM_file") || die "Unable to open $SAM_file: $!";
    
    my (%samhash,%avprofilehash,$read_count);
    
    LINE: while(<$sam_in>) {
        
        chomp;

        # detect header lines
        my $curr_line = $_;
        if ($curr_line =~ /^@/) {
            
            # if reference sequence line
            if ( $curr_line =~ /\@SQ.+SN:(chr[\dRMFBIV]+).+LN:(\d+)/ ) {
                
                # get sequence name and length
                my $chrn = $1;
                my $length = nearest($pos_bin,$2);
                
                # loop through all chromosome bins
                for ( map {$_*$pos_bin} 0..($length/$pos_bin) ) {
                    
                    my $pos = $_;
                    
                    # loop through all fragment sizes possible
                    for ( map {$_ * $frag_bin} ($min_frag/$frag_bin)..($max_frag/$frag_bin) ) {
                        
                        my $frag_size = $_;
                        $samhash{$chrn}{$pos}{$frag_size} = 0;
                    }
                }
                
            }
            next LINE;
        }
        
        my @line = split('\t');
        
        # Check for chromosome alignment and skip if unaligned
        my $chrn;
        if ($line[2] =~ /chr[\dRBMFVI]+/) { 
            $chrn = $line[2];
        } elsif ($line[2] =~ /(NC_\d+)/) { # handle refseq references in stead of chromosome names using a sub
            my $refseq = $1;
            $chrn = &refseqname($refseq);
        } 
        
        next LINE unless ($chrn);
        
        # Determine fragment size bin according to user specified bins
        my $frag_size;
        if ($PE) {
            $frag_size = nearest($frag_bin,$line[8]);
        } else {
            my $read =  $line[10];
            $frag_size = nearest($frag_bin,length($read));
        }
        
        # Check fragment length is within window
        if ($frag_size >= $min_frag && $frag_size <= $max_frag) {
            
            # Determine dyad position
            my $dyad_pos = nearest($pos_bin,($line[3] + ($frag_size * 0.5)));
            
            # Increment the count for caluclated fragment size bin at dyad position
            $samhash{$chrn}{$dyad_pos}{$frag_size}++;
            $samhash{$chrn}{$dyad_pos}{"Total"}++;
            
            # If higher resolution wanted also increment for this binning
            if ($av_bin) {
                my $av_fragsize;
                if ($PE) {
                    $av_fragsize = nearest($av_bin,$line[8]);
                } else {
                    $av_fragsize = nearest($av_bin,length($line[10]));
                }
                $avprofilehash{$chrn}{$dyad_pos}{$av_fragsize}++;
            }
            
            # Increment count of reads within window
            $read_count++; 
        }
    }
    
    print "Total reads within specified fragment range: $read_count\n";
    
    # open sgr output file
    my $sgrname = $condition."_$min_frag"."-"."$max_frag"."bp.sgr";
    open(my $sgrout,'>',"$outdir/$sgrname") || die "Unable to open output file: $sgrname\n";			
    
    # 3 bin smoothing average print sgr
    for my $chrn (sort keys %samhash) {
        for my $pos (sort {$a <=> $b} keys %{$samhash{$chrn}}) {
            # smooth sgr sum with 3 bin average and print
            my $sum_val = nearest(0.0001,( (($samhash{$chrn}{$pos-$pos_bin}{"Total"}||0)+($samhash{$chrn}{$pos}{"Total"}||0)+($samhash{$chrn}{$pos+$pos_bin}{"Total"}||0))/3 ));
            print ($sgrout "$chrn\t$pos\t$sum_val\n");
        }
    }
    
    # open peaks files
    my @suffixes = (".txt",".sgr");
    my ($peaks_filename, $peaks_dir, $peaks_suffix) = fileparse($peaks_file,@suffixes);
    die "Error: $peaks_filename doesn't appear to be a .txt or .sgr file\n" unless ($peaks_suffix);
    
    if ( $peaks_filename =~ /^$condition\_([A-Za-z\d\-]{1,10}).+_(t\d+)/ ) {
        
        my $size_range = $1;
        my $peaks_thresh = $2;
        
        # open output files
        print "Found and processing $peaks_filename \n";
        open(my $peaks_in, '<', "$peaks_file") || die "Unable to open $peaks_file: $!";

        my $c3_outfile = "$condition\_$size_range\_$peaks_thresh\_frag_size_C3.txt";
        open(my $c3out, '>', "$outdir/$c3_outfile")|| die "Unable to open $c3_outfile: $!";
        
        my $peaks_outfile = "$peaks_filename\_size.txt";
        open(my $peaks_out, '>', "$outdir/$peaks_outfile")|| die "Unable to open $peaks_outfile: $!";
        
        # print header for C3 file
        my @header = map{$frag_bin*$_} ($min_frag/$frag_bin)..($max_frag/$frag_bin);
        print ($c3out "ID_string\t".join("\t", @header)."\n");
        
        # initiate some variables 
        my (%sum_average,@hi_res_sum,$total_frag_sum,%peakshash); 
        my $peaks_count = 0;
        
        # loop through peaks file
        PEAK: while(<$peaks_in>) { 
            
            chomp;
            my @line = split('\t');
            
            my $chrn;
            if ($line[0] =~ /^chr[\dMRFBVXI]+/ ) {
                $chrn = $line[0];
            } elsif ($peaks_count == 0) { #if header line, print to new peaks file
                print ($peaks_out join("\t",@line)."\tSize\n");
                next PEAK;
            }
            
            my $pos = $line[1];
            
            # get fragment sizes surrounding each peak
            my (%nuc_hash,$frag_sum); # reset sum for individual nucleosome region
            for (my ($i,$k) = ($pos-$pos_window,-$pos_window);$i <= ($pos + $pos_window);$i+=$pos_bin,$k+=$pos_bin) { #foreach pos around peak
                
                # Sum fragments for this nucleosome
                for ( my $fragsize = $min_frag; $fragsize <= $max_frag; $fragsize+=$frag_bin ) {
                    $nuc_hash{$fragsize} += $samhash{$chrn}{$i}{$fragsize} ? $samhash{$chrn}{$i}{$fragsize} : 0;
                    $frag_sum += $samhash{$chrn}{$i}{$fragsize} ? $samhash{$chrn}{$i}{$fragsize} : 0;
                }
                
                # Sum fragments for average
                for (my ($fragsize,$index) = ($min_frag,0);$fragsize<=$max_frag;$fragsize+=$frag_bin,$index++) {
                    $sum_average{$k}[$index] += $samhash{$chrn}{$i}{$fragsize} ? $samhash{$chrn}{$i}{$fragsize} : 0;
                    $total_frag_sum += $samhash{$chrn}{$i}{$fragsize} ? $samhash{$chrn}{$i}{$fragsize} : 0;
                    $sum_average{"total"}[$index] += $samhash{$chrn}{$i}{$fragsize} ? $samhash{$chrn}{$i}{$fragsize} : 0;
                }
                
                # Sum fragments for average at higher res if wanted
                if ($av_bin) { 
                    for (my ($fragsize,$index) = ($min_frag,0);$fragsize<=$max_frag;$fragsize+=$av_bin,$index++) {
                        $hi_res_sum[$index] += $avprofilehash{$chrn}{$i}{$fragsize} ? $avprofilehash{$chrn}{$i}{$fragsize} : 0;
                    }
                }
            }
            
            my ($max_size,@nuc_sum);
            
            # skip if no reads mapping to this pos (shouldn't happen but just incase!)
            if (!defined $frag_sum || $frag_sum == 0) {
                
                next PEAK;
            
            } else {
                
                # normalise this nucleosomes values by total frag count
                for my $fragsize (sort {$a <=> $b} keys %nuc_hash) {
                    
                    $nuc_hash{$fragsize} = $nuc_hash{$fragsize}/$frag_sum;
                    
                }
                
                # 3 bin smooth of frag profile for this peak						
                my ($temp_count,$max_val) = (0,0);
                
                for my $fragsize (sort {$a <=> $b} keys %nuc_hash) {
                    
                    my $temp_sum = ( ($nuc_hash{$fragsize-$frag_bin}//0) + ($nuc_hash{$fragsize}//0) + ($nuc_hash{$fragsize+$frag_bin}//0) );
                    $nuc_sum[$temp_count] = nearest( 0.00001, ($temp_sum/3) );
                    
                    # keep track of max value
                    ($max_size,$max_val) = ($fragsize,$nuc_hash{$fragsize}) if ($nuc_hash{$fragsize} > $max_val);
                    
                    $temp_count++;
                    
                }
                
                $peaks_count++;
                
            }
            
            # print each nucleosome frag profile to C3 file with an ID made up of chromosome and pos
            print ($c3out "$chrn"."-"."$pos\t".join("\t",@nuc_sum)."\n");
            
            # print peaks out again with maximum fragment size added to last column
            print ($peaks_out join("\t",@line)."\t$max_size\n");
        }
        
        print "$peaks_count peaks found and processed from $filename\n";
        
        # print 2d average profile
        if ($av_bin) {
            
            @header = map{$av_bin*$_} ($min_frag/$av_bin)..($max_frag/$av_bin);
            
            my $outfile = $condition."-".$peaks_thresh."-".$min_frag."-".$max_frag."bp_2D_fragment_profile.txt";
            open (my $av_out, '>', "$outdir/$outfile")|| die "Unable to open $outfile: $!";
            
            print($av_out "Fragment Size\tFreq\n");
        
            for (0..$#hi_res_sum) {
                my $val = nearest(0.000001,($hi_res_sum[$_]/$total_frag_sum));
                print ($av_out "$header[$_]\t$val\n");
            }
        } 
        
        # print 3d profile
        @header = map{$frag_bin*$_} ($min_frag/$frag_bin)..($max_frag/$frag_bin);
        
        my $outfile = $condition."-".$peaks_thresh."-".$min_frag."-".$max_frag."bp_3D_fragment_profile.txt";
        open (my $profile_out, '>', "$outdir/$outfile")|| die "Unable to open $outfile: $!";
        
        print ($profile_out "3D average fragment profile for $condition with fragments between $min_frag and $max_frag bp\n");
        print ($profile_out "SAM file used: $filename\n");
        print ($profile_out "Peaks file used: $peaks_filename\n");
        print ($profile_out "Total reads mapped to peak regions: $read_count\n");
        print ($profile_out "Total number of peaks with mapped reads: $peaks_count\n\n");
        print ($profile_out "Position\t".join("\t",@header)."\n");
        
        my @positions = map{$_*$pos_bin} (-$pos_window/$pos_bin)..($pos_window/$pos_bin);
        for (0..$#positions) {
            my $pos = $positions[$_];
            print($profile_out "$pos");
            for (0..$#{$sum_average{$pos}} ) {
                my $val = nearest(0.000001,($sum_average{$pos}[$_]/$total_frag_sum));
                print($profile_out "\t$val");
            }
            print($profile_out "\n");
        }
    }
}

print "End: ".`date`."\n";

sub refseqname {
	# converts refseq id to chromosome name
	my $temp;
	$temp = "chrI" if ($_[0] eq "NC_001133");
	$temp = "chrII" if ($_[0] eq "NC_001134");
	$temp = "chrIII" if ($_[0] eq "NC_001135");
	$temp = "chrIV" if ($_[0] eq "NC_001136");
	$temp = "chrV" if ($_[0] eq "NC_001137");
	$temp = "chrVI" if ($_[0] eq "NC_001138");
	$temp = "chrVII" if ($_[0] eq "NC_001139");
	$temp = "chrVIII" if ($_[0] eq "NC_001140");
	$temp = "chrIX" if ($_[0] eq "NC_001141");
	$temp = "chrX" if ($_[0] eq "NC_001142");
	$temp = "chrXI" if ($_[0] eq "NC_001143");
	$temp = "chrXII" if ($_[0] eq "NC_001144");
	$temp = "chrXIII" if ($_[0] eq "NC_001145");
	$temp = "chrXIV" if ($_[0] eq "NC_001146");
	$temp = "chrXV" if ($_[0] eq "NC_001147");
	$temp = "chrXVI" if ($_[0] eq "NC_001148");
	$temp //= "";
	return $temp;
}
