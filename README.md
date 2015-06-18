# Marks_scripts

## Introduction 

General approach was to identify, characterise & plot peaks individually while maintaining some level of positioning information. This eliminates some of the confounding issues from plotting genomic sites and allows seperation of the four main nucleosome parameters I was interested in: positioning relative to genomic sites (position), uniformity of positioning accross the population (distribution), peak height (occupancy) and size of protected region (size).

## Obtaining & using scripts

I converted all of these scripts to command line tools and tried to keep them as flexible as possible to work with different organisms and accept customisable settings. I *think* that if you just clone the repository and add it to your path you should be good to go?

One rather in-flexible part is the filenaming - uses regexp to check correct file input so your SAM file need to start with a description with the following format: 
    `[A-Za-z\d]{1,10}_[A-Za-z\d]{1,10}_[A-Za-z\d]{1,10}`
    (e.g. WT_rep1_ExpA_blahbahblah.sam). Nuc_mapper.pl will add the particle size range - don't delete this.

Using options - I used the perl Getopt::Long for command line options, which is basically like usual command line flags with either - or -- but don't group them (i.e. `-f -b` or `--flag --b` good, `-fb` bad!). All options have a default value which I tend to use for Dicty but it often won't make any sense for your data, make sure to check what each option does!'

## Nuc_mapper.pl 

Simplified version of Nicks mapping method - just takes SAM file and converts straight to an .sgr histogram with a fragment size range of your choosing. Includes option to output SAMparser style files if needed for other Kent lab scripts. Will just plot any aligned reads in the SAM - so make sure your SAM only contains the alignments you want to contribute to the final dyad histogram! 

**Output:** Dyad histogram file (.sgr).

**Usage:** perl Nuc_mapper_vn.pl --options SAMfile.sam

**Arguments:**

    --bin|-b = binning desired (default = 10bp)
    --min|-mi = minimum fragment length to include (after rounding to closest bin) (default = 120bp)
    --max|-ma = maximum fragment length to include (after rounding to closest bin) (default = 180bp)
    --SAMparser|-s = flag SAMparser style output on (default = off)
    --PE|-p = flag paired end reads (default = off)
    --out|-o = output directory (default = current working directory)

**Example:** `Nuc_mapper_v1.pl --PE -mi 100 -ma 190 --bin 5 /path/WT_rep1_0H_aligned.sam`

## Nuc_caller.pl

Basically just Nick's peak marker with a few changes: requires two bins either side of central bin to have lower dyad frequency, has a maximum threshold as well as minimum, prevents dual peaks being called within n bins of each other. Output peak height is 3 bin sum around peak. 

v2 edit: Added ability to measure distribution of peaks and treat disorganised regions as a single peak in order to avoid calling lots of false peaks in noisy regions. 

NB: I usually smooth and then normalise all sgrs to same read-depth before doing peak calling but can be run on un-smoothed and use scaling factors if prefered.

**Output:** Peak positions and height (.sgr).

**Usage:** perl Nuc_caller_vn.pl --options dyad_histogram.sgr

**Arguments:**

    --thresh|-t = minimum dyad frequency for calling peaks (default = 10)
    --max|-m = maximum dyad frequency of peaks (default = 200)
    --scale|-s = scaling factor (default = 1.0)
    --dual|-d = minimum number of bins between called peaks (default = 2)
    --grad|-g = minimum gradient threshold - below which peaks treated as a "disorganised region" rather than a single peak, only checks gradient if set (default = off)
    --window|-w = window over which to calc gradient (default = 5)
    --out|-o = output directory (default = current working directory)

**Example:** `perl Nuc_caller_v1.pl -t 20 -m 1000 /path/WT_rep1_0H_120-180bp_histogram_10.sgr`

## Nuc_dist.pl

Takes cumulative sum of dyads accross a specified window surrounding each peak and measures gradient as a measure of dyad distribution. Also allows elimination of peaks within n bins of each other and those with low read depth (to deal with artifacts from smoothing).

I would reccommend using un-smoothed data so that distribution hasn't been altered.

NB: Input peaks file can be .txt or .sgr to allow use of modified peaks files.

**Output:** 

1. A new peaks file with same format as input but with additional column containing gradient, and without excluded peaks (_grad.txt). 

2. An .sgr of overlapping peak locations for checking in IGB etc.

**Usage:** perl Nuc_dist_vn.pl --options dyad_histogram.sgr peaks_file.sgr

**Arguments:** 

    --bin|-b = binning (default = 10)
    --pwind|-p = peak window - bins either side of peak to sum dyads over (default = 5)
    --gwind|-g = gradient window - bins around peak over which to measure gradient from summed dyads (default = 2)
    --ewind|-e = exclusion window - minimum distance between peaks (default = 60) 
    --min|-m = minimum number sum of dyads within peak window to keep peak (default = 10)
    --out|-o = output directory (default = current working directory)

**Example:** `perl Nuc_dist_v1.pl --pwind 10 -g 5 -e 50 /path/WT_rep1_0H_120-180bp_histogram_10.sgr /path/WT_rep1_0H_120-180bp_peaks_t10.sgr`

## Nuc_size.pl

Creates a histogram of fragment sizes contributing to a given peak. Will define a peak area (peak position +/- n bins) and sum all fragments within specified size range in that region. Can be used to check fragment sizes outside of that used for the dyad histogram from which peaks were called.

NB: pretty slow since has to re-map from SAM file to get different fragment sizes

**Output:** 

1. Adds max/mode fragment size to peaks file.

2. C3 type file for clustering with a fragment profile for each individual peak. 

3. A new .sgr file of mapped fragment sizes - to check in IGB or similar.

4. 2D and 3D profiles (2D = just fragment sizes, 3D =  position and frag size).

**Usage:** perl Nuc_size_vn.pl --options SAMfile.sam peaksfile.sgr 

**Arguments:**

    --wind|-w = bp either side of peak to count as peak region (default = 50)
    --min|-mi = minimum fragment size to include (default = 100)
    --max|-ma = maximum fragment size to include (default = 200)
    --fbin|-f = binning for fragment sizes (default = 5)
    --pbin|-p = binning for genomic position (default = 10)
    --abin|-a = optional binning for the average 2D/3D profile so can set seperate binning for c3 and average if desired (default = off)
    --PE|-p = flag paired end reads (default = off)
    --out|-o = output directory (default = current working directory)

**Example:** `perl Nuc_size_v1.pl --a 1 --wind 75 WT_rep1_0H_aligned.sam WT_rep1_0H_120-180bp_peaks_t10.sgr`

## Nuc_changes

Will compare two peaks files against one another and measure change in each provided peak parameter for any peaks occupying same genomic position between conditions (peak pos +/- n bp). Calls all changes between common peaks over given thresholds as nucleosome changes.

Will measure whatever is present in peaks file - e.g. if just use unmodified peaks file from Nuc_caller.pl/PeakMarker.pl will just measure positioning changes and height changes but not gradient or size.

NB: Default settings were set as ~5x SD of the changes observed between Dicty WT conditions - I'd reccommend running this with WT reps first and determining appropriate thresholds.

**Output:** 

1. One .sgr file for each change present in provided peaks file
        -Units for changes: pos and size = bp, grad, height and occ = log2 fold change
        -Occupancy change = peaks with height change AND no gradient change

2. Info file

**Usage:** perl Nuc_changes_vn.pl reference_peaks_file.txt test_peaks_file.txt

**Arguments:**

    --bin|-b = binning (default = 10)
    --common|-c = bp window either side of peak to look for common peaks between conditions (default = 40)
    --pos|-p = bp position change threshold for calling peaks with altered positioning (default = 10)
    --size|-s = bp size change threshold for calling peaks with altered fragment size (default = 20)
    --grad|-g = fold gradient change threshold for calling peaks with altered distribution (default = 1.1)
    --height|-h = fold height change threshold for calling peaks with altered height (default = 1.3)
    --out|-o = output directory (default = current working directory)

**Example:** `perl Nuc_changes_v1.pl -p 25 -s 40 WT_rep1_ExpA_120-180bp_peaks_dist_size.txt Mutant_rep1_ExpA_120-180bp_peaks_dist_size.txt`

## Nuc_plotter

Maps identified peaks back to genes, and determines nucleosome order, before plotting the average dyad distributions for each specified position. Uses the gene termini as boundaries and allows you to plot as many nucleosomes as you like upstream/downstream of these boundaries (e.g. you could plot the first 3 nucleosomes downstream of the TSS, one upstream and none around TTS). 

Also allows you to input peaks for elimination - the output will include the profile of total peaks then the profiles after eliminating each input elimination set. 

**Input:**

-Peaks for elimination should be in format from Nuc_changes.pl - i.e. basically just a peaks file but with `*_changename_changes.sgr` naming and a header with change name.

-Coordinates file - you basically want to parse a GFF type file to get whatever genome locations you want but needs to have following format: "chrn\tstart\tend\tID\tstrandedness\n"

**Rules:** This plotting method is a little complicated so to clarify exactly whats happening without you having to go through this ridiculous script...
    
-Uses half the number of nucleosomes mapped to that gene as a limit - e.g. if you set the number of nucleosomes plotted at the start of the gene to be 3 - a gene with 6+ nucleosomes will plot the first three, whereas a gene with 4 nucleosomes will only plot 2.
    
-Any nucleosomes not plotted as a defined position are added to the average (there are three averages - upstream av, mid av & downstream av) (e.g. if gene has 5 genes and set plotting of 3 from each end, will just plot 2 and add the middle nucleosome to mid average)
    
-Intergenic nucleosomes can be included multiple times - e.g. if a nucleosome is at the -1 of a downstream gene and TN+1 position of an upstream gene, it is included in the average profile of each position
    
-Peaks are plotted in a gene-relative orientation
   
-Will include all nucleosomes up to a specified distance from each gene as intergenic nucleosomes belonging to that gene, but won't go into coding regions.

**Output:**

1. Average genic and intergenic nucleosome profiles along with confidence intervals

2. Average nucleosome position profile including each specified position

3. An .sgr of determined peak position categories if you want to check you agree with called positions

**Usage:** perl Nuc_plotter_vn.pl dyad_histogram_file.sgr peaks_file.sgr gene_coord_file.txt [elim_files.sgr (optional)]
NB: for elimination files, just drag all change.sgr files wanted into the terminal window after other arguements or list (each file seperated by a space)

**Arguments:**
    
    --bin|-b = binning (default 10)
    --window|-w = bp either side of peak to plot (default = 50)
    --nrl|-n = minimum distace between peaks (default = off)
    --add|-a = bp distance to extend into intergenic regions flanking gene for finding intergenic peaks belonging to a given gene (default = 1500)
    --up|-u = number of nucleosomes to plot upstream of TSS/start boundary (default = 1)
    --start|-s = number of nucleosomes to plot downstream of TSS/start boundary (default = 3)
    --end|-en = number of nucleosomes to plot upstream of TTS/end boundary (default = 2)
    --down|-d = number of nucleosomes to plot downstream of TTS/end boundary (default = 1)
    --ci|-c = % confidence interval to calculate - (I only included crit values for 90%, 95% and 99% so choose one of those (or add it!)) (default = 95)
    --out|-o = output directory (default = current working directory)

**Example:** `perl Nuc_plotter_v1.pl Mutant_rep1_expA_120-180bp_10.sgr Mutant_rep1_expA_120-180bp_peaks_t10.sgr gene_coordinates_2009.txt mutant_vs_WT_dist_changes.sgr mutant_vs_WT_occ_changes.sgr`
