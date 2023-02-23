use strict;
use warnings;
use lib "/users/lvelten/lvelten/Software/perl/lib/perl5/x86_64-linux-thread-multi/";
use String::Approx 'amatch';
use List::Util qw(min max sum);

#What this script should do
#0. Read in the filtered assignments from the previous step
#1. Check that the fwd and the rev. BC read ar the same
#2. Check if it matches any of the assigned BCs perfectly.
#3. (optional) if not, use kmers to rapidly find a match and then amatch to confirm - currently count these to see how big the problem is.
     #the problem here is that even with an editing distance of 1 there can be several matches etc. require an extact match, but use the larger barcode list.
#4. Remember the UMI that this BCs has and the number of reads for that UMI
#5. Merge similar UMIs based (remove more lowly covered UMIs)
#6. Report statistics - how covered is each UMI (sequencing saturation)
#7. Output number of UMIs per BC, assignment

my ($umi, $fwd, $rev, $assignments, $outfile, $statfile, $nonmapped) = @ARGV;

my @MINREADUMI = (1..10); #report UMIs with at least this many reads
my $MAXREADUMI = 10000000; #set to a high value in order to not filter
my $downsample = 0; #set to 1 if downsampling of reads should be performed
my $dsto = 0.5; #fraction of reads to downsample

###########################################
###  Read in  filtered assignments.    ####
###########################################

my %BC_CRS;

if ($assignments =~ /gz$/) {
      open(ASSI, "zcat $assignments|");
  } else {
      open(ASSI, "<$assignments");
  }

my $line = <ASSI>; #skip header
while($line = <ASSI>) {
  chomp $line;
  my @elements = split /,/ , $line;
  $BC_CRS{$elements[0]} = $elements[1];
}
close(ASSI);

###########################################
###          Go through reads          ####
###########################################

sub rc {
  my $revcomp = reverse $_[0];
  $revcomp =~ tr/ATGCatgc/TACGtacg/;
  return $revcomp;
}

open(UMI, "zcat $umi|"); my $umiline;
open(FWD, "zcat $fwd|"); my $fwdline;
open(REV, "zcat $rev|"); my $revline;
open(NONMAP, "|gzip -c > $nonmapped ");
my $i = 0;
my $counter = 0;
my $no_perfect_match = 0;
my $perfect = 0;
my $readname;

my %BC_UMI; #hash that maps BC to UMIs
#initialize this hash with known BC, then there will be no chaos down below.
#eahc element in this hash is a hash of UMIs (value: number of reads for that BC - UMI association)

foreach my $bc (keys %BC_CRS) {
  $BC_UMI{$bc} = {};
}

while ($umiline = <UMI>) {
  chomp $umiline;
  $fwdline = <FWD>; chomp $fwdline;
  $revline = <REV>; chomp $revline;
  
  $i = 0 if ($i == 4); 
  $i++;
  
  if ($i == 1) {
  $umiline =~ s/\s.+//;
  $fwdline =~ s/\s.+//;
  $revline =~ s/\s.+//;
  $readname = $umiline;
  
  die "FASTQ files not matching\n" unless ($readname eq $fwdline and $readname eq $revline);

  }
  
  next unless ($i == 2); #skip read name, quality string etc. Only care about the actual read
  $counter++;
  warn "$counter\n" if ($counter % 1000000 == 0);
  #downsampling
  next if (rand() > $dsto and $downsample);
  #check if fwd BC read = rev BC read
  $fwdline = substr $fwdline, 0, 15;
  $revline = substr rc($revline), 1,15;
  
  if ($fwdline ne $revline) { #if they map to different CRS: weird
    if (defined($BC_CRS{$fwdline}) and defined($BC_CRS{$revline})) {
      print NONMAP "$fwdline,$revline,weird\n";
    } elsif (defined($BC_CRS{$fwdline})) {
      $BC_UMI{$fwdline}->{$umiline}++;
      $perfect++;
    } elsif (defined($BC_CRS{$revline})) {
      $BC_UMI{$revline}->{$umiline}++;
      $perfect++;
    } else {
      $no_perfect_match++;
    }
  } else { #if the fwd and the rev. read are the same
    if (defined($BC_CRS{$fwdline})) {
      $BC_UMI{$fwdline}->{$umiline}++;
      $perfect++;
    } else {
      #allow for one mismatch
      $no_perfect_match++;
      print NONMAP "$fwdline,$revline,normal\n";
    }
  }
  
}


close(UMI);
close(FWD);
close(REV);
close(NONMAP);

warn "Perfect: $perfect , Not: $no_perfect_match\n";

###########################################
###            Merge UMIs              ####
###########################################
warn "Now merging UMIs...";
my %BC_UMI_fixed; 
while ( my ($barcode, $umis) = each (%BC_UMI)) {
  $BC_UMI_fixed{$barcode} = {} if (%{$umis});
}

my $total = scalar(keys %BC_UMI_fixed);
$counter =0;
while ( my ($barcode, $umis) = each (%BC_UMI)) {
  #sort UMIs by number of reads
  next unless %{$umis};
  $counter++;
   warn "$counter of $total\n" if ($counter % 1000000 == 0);
  my @sorted = sort {$umis->{$a} <=> $umis->{$b} } keys %{$umis}; #sort umis by number of reads
  while (my $this_umi = shift @sorted) { 
    #starting with the UMI backed by the smallest number of reads
    
     my @MATCHES;
     if (@sorted) {
        @MATCHES = amatch($this_umi,  ["S1","I0","D0"], @sorted); #find approx. matches with at least one substitution in all barcodes with more reads
     }
     
    if (@MATCHES) {
      @MATCHES = sort { $umis->{$b} <=> $umis->{$a}} @MATCHES; #if any are found, sort matches by barcode counts
      if (defined($BC_UMI_fixed{$barcode}->{$MATCHES[0]})) {
        $BC_UMI_fixed{$barcode}->{ $MATCHES[0] }  += $umis->{$this_umi};
      } else {
        $BC_UMI_fixed{$barcode}->{ $MATCHES[0] } = $umis->{$this_umi};
      }
      #warn "Replace $this_umi by $MATCHES[0]\n";
    } else {
      if (defined($BC_UMI_fixed{$barcode}->{$this_umi})) {
          #if for a certain barcode a barcode with lower counts has been mapped, add the statistics
          $BC_UMI_fixed{$barcode}->{$this_umi} += $umis->{$this_umi}; #add read counts
          } else {
          $BC_UMI_fixed{$barcode}->{$this_umi} = $umis->{$this_umi}; #add read counts
        }
    }
  }
  
}


###########################################
###           print output to csv      ####
###########################################
open(OUT, "|gzip -c >$outfile");
print OUT "BARCODE,CRS,UMI,READS";
foreach my $thresh (@MINREADUMI) {
  print OUT ",UMI$thresh";
}
print OUT "\n";

my %reads_per_umi; #a histogram: How often are there x reads per UMI?
while ( my ($barcode, $umis) = each (%BC_UMI_fixed)) {
  my @umi_reads = values %{$umis};
  my $numi = scalar(@umi_reads);
  my $nreads = 0;
  my @numinread = (0) x @MINREADUMI;
  foreach my $nn (@umi_reads) {
    $reads_per_umi{$nn}++;
    $nreads += $nn;
    for (my $ii = 0; $ii <= $#MINREADUMI; $ii++) {
      $numinread[$ii]++ if ($nn >= $MINREADUMI[$ii] and $nn <= $MAXREADUMI);
    }
    #$numioneread++ if ($nn == 1);
    #$numimanyread++ if ($nn >= $MINREADUMI);
    #$numifiveread++ if ($nn >= 5);
    #$numitenread++ if ($nn >= 10);
  }
  print OUT "$barcode,$BC_CRS{$barcode},$numi,$nreads";
  foreach my $val (@numinread) {
    print OUT ",$val";
  }
  print OUT "\n";
}

close(OUT);

open(STAT, "|sort -n>$statfile");
while (my ($r, $c) = each (%reads_per_umi)) {
  print STAT "$r,$c\n";
}
close(STAT);