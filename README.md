# Marks_scripts

## Introduction 

General approach was to identify, characterise & plot peaks individually while maintaining some level of positioning information. This eliminates some of the confounding issues from plotting genomic sites and allows seperation of the four main nucleosome parameters I was interested in: positioning relative to genomic sites (position), uniformity of positioning accross the population (distribution), peak height (occupancy) and size of protected region (size).

## Obtaining & using scripts

I converted all of these scripts to command line tools and tried to keep them as flexible as possible to work with different organisms and accept customisable settings. I *think* that if you just clone the repository and add it to your path you should be good to go?

One rather in-flexible part is the filenaming - uses regexp to check correct file input so your SAM file need to start with a description with the following format: 
    `[A-Za-z\d]{1,10}_[A-Za-z\d]{1,10}_[A-Za-z\d]{1,10}`
    e.g. WT_rep1_ExpA_blahbahblah.sam 
Nuc_mapper.pl will add the particle size range - don't delete this.

Using options - I used the perl Getopt::Long for command line options, which is basically like usual command line flags with either - or -- but don't group them (i.e. `-f -b` or `--flag --b` good, `-fb` bad!)

## Nuc_mapper.pl 

Simplified version of Nicks mapping method - just takes SAM file and converts straight to an .sgr histogram with a fragment size range of your choosing. Includes option to output SAMparser style files if needed for other Kent lab scripts. Will just plot any aligned reads in the SAM - so make sure your SAM only contains the alignments you want to contribute to the final dyad histogram! 

Output: Dyad histogram file (.sgr).

Usage: perl Nuc_mapper_vn.pl --options SAMfile.sam

Arguments:

    --bin|-b = binning desired (default = 10bp)
    --min|-mi = minimum fragment length to include (after rounding to closest bin) (default = 120bp)
    --max|-ma = maximum fragment length to include (after rounding to closest bin) (default = 180bp)
    --SAMparser|-s = flag SAMparser style output on (default = off)
    --PE|-p = flag paired end reads (default = off)
    --out|-o = output directory (default = current working directory)

Example: `Nuc_mapper_v1.pl --PE -mi 100 -ma 190 --bin 5 /path/WT_rep1_0H_aligned.sam`

## Nuc_caller.pl

Basically just Nick's peak marker with a few changes: requires two bins either side of central bin to have lower dyad frequency, has a maximum threshold as well as minimum, prevents dual peaks being called within n bins of each other. Output peak height is 3 bin sum around peak.

Output: Peak positions and height (.sgr).

Usage: perl Nuc_caller_vn.pl --options dyad_histogram.sgr

Arguments:

    --thresh|-t = minimum dyad frequency for calling peaks (default = 10)
    --max|-m = maximum dyad frequency of peaks (default = 200)
    --scale|-s = scaling factor (default = 1.0)
    --dual|-d = minimum number of bins between called peaks (default = 2)
    --out|-o = output directory (default = current working directory)

Example: `perl Nuc_caller_v1.pl -t 20 -m 1000 /path/WT_rep1_0H_120-180bp_histogram_10.sgr`

## Nuc_dist.pl

Takes cumulative sum of dyads accross a specified window surrounding each peak and measures gradient as a measure of dyad distribution. Also allows elimination of peaks within n bins of each other and those with low read depth (to deal with artifacts from smoothing).

    I would reccommend using un-smoothed data so that distribution hasn't been altered.

    NB: Input peaks file can be .txt or .sgr to allow use of modified peaks files.

Output: A new peaks file with same format as input but with additional column containing gradient, and without excluded peaks (_grad.txt). Also outputs .sgr of overlapping peak locations for checking in IGB etc.

Usage: perl Nuc_dist_vn.pl --options dyad_histogram.sgr peaks_file.sgr

Arguments: 

    --bin|-b = binning (default = 10)
    --pwind|-p = peak window - bins either side of peak to sum dyads over (default = 5)
    --gwind|-g = gradient window - bins around peak over which to measure gradient from summed dyads (default = 2)
    --ewind|-e = exclusion window - minimum distance between peaks (default = 60) 
    --min|-m = minimum number sum of dyads within peak window to keep peak (default = 10)
    --out|-o = output directory (default = current working directory)

Example: `perl Nuc_dist_v1.pl --pwind 10 -g 5 -e 50 /path/WT_rep1_0H_120-180bp_histogram_10.sgr /path/WT_rep1_0H_120-180bp_peaks_t10.sgr`

## Nuc_