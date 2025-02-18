#!/usr/bin/perl -w
###
# Process output from calcDxy.R to get windowed Dxy
# AUTHOR: Katharine Korunes
# MODIFIED: Michael Connelly, Oct. 2024 (w/ help from ChatGPT)
# USAGE: perl angsd_dxy_step3_processoutputoverwindows_chromosomes.pl Dxy_persite.txt chromosome_lengths.txt window_size combined_angsd_DxySummary.txt
###
use strict;

my $win = "$ARGV[2]"; # Window size
my $chrom_file = "$ARGV[1]"; # File with chromosome names and lengths
my $dxy_file = "$ARGV[0]"; # Input file with Dxy per site
my $combined_output = "$ARGV[3]"; # Combined output file

# Open the combined output file
open (COMBINED_OUT, ">$combined_output") or die "file not found $!";

# Print header to the combined output file
print COMBINED_OUT "Chrom\tWindowStart\tWindowEnd\tAvgOverAllSites\tAvgOverSitesWithDxyEstimate\n";

# Read chromosome lengths from the input file
open (CHROMFILE, $chrom_file) or die "file not found $!\n";
my %chrom_lengths;
while (<CHROMFILE>) {
    chomp;
    my ($chrom, $len) = split(/\s+/); # Assuming two columns: chromosome name and length
    $chrom_lengths{$chrom} = $len;
}
close (CHROMFILE);

# Process each chromosome
foreach my $chrom (keys %chrom_lengths) {
    my $totalLen = $chrom_lengths{$chrom}; # Get chromosome length
    
    my %snps;
    open (CHR, $dxy_file) or die "file not found $!\n";
    while (<CHR>){
        my $line = $_;
        next if ($line =~ m/^chromo/); # Skip header
        my @fields = split(/\s+/, $line);
        my $chr_in_file = $fields[0]; # Assuming chromosome name is the first column
        next unless ($chr_in_file eq $chrom); # Only process lines for the current chromosome
        my $position = $fields[1];
        my $dxy = $fields[2];
        $snps{$position} = $dxy;
    }
    close (CHR);

    # Now get avg from each window and write to combined file
    for (my $i = 1; $i <= ($totalLen - $win); $i += $win){
        my $start = $i;
        my $end = ($i + ($win - 1));
        my $sum = 0;
        my $withData = 0;
        for (my $x = $start; $x <= $end; $x++){
            if (exists $snps{$x}){
                my $val = $snps{$x};
                $withData++;
                $sum += $val;
            }
        }
        my $avg = ($sum / $win);
        my $avgSitesWithData = 0;
        if ($withData > 0){
            $avgSitesWithData = ($sum / $withData);
        }
        print COMBINED_OUT "$chrom\t$start\t$end\t$avg\t$avgSitesWithData\n";
    }
}

close COMBINED_OUT;
exit;
