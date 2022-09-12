use strict;
use warnings;
use lib "/users/lvelten/lvelten/Software/perl/lib/perl5/x86_64-linux-thread-multi/";
use String::Approx 'amatch';
use List::Util qw(min max sum);

my ($bc, $bam, $outfile) = @ARGV;



#################################################
####      Mapping barcodes and alignments    ####
#################################################

open(BC, "zcat $bc|"); my $bcline;
open(BAM, "samtools view $bam| "); my $bamline;

my %BC_CRS_raw; #hash mapping barcode to CRS
my ($readname);
my $lastread = "";
my $i = 0;
my $counter = 0;
while ($bcline = <BC>) {

  $i = 0 if ($i == 4); 
  $i++;
  
  chomp $bcline;
  
  $readname = $bcline if $i == 1;
  next unless ($i == 2); #skip read name, quality string etc. Only care about the actual read
  $counter++;
  $readname =~ s/^@(\S+)\s.+/$1/;
  
  my @bamsplit; #split bam entries
  my($fwdmap, $revmap, $match);
  my $j = 0;
  my $multiple_alignments = 0;
  #### read in 2 lines from the bam file (fwd and rev)
  while ($bamline = <BAM>) {
    @bamsplit = split /\t/, $bamline;
    if ($bamsplit[0] eq $lastread) { #if the read name in bam is behind the readname in the barcode fastq file, skip a few lines in the bam file
      #warn "$counter: Mutliple alignments for $lastread, using first\n";
      next;
    }
    $j++;
    #make sure that there are no reads in the barcode fastq file that do not appear in the bam file
    die "bam and fastq not matching: $bamsplit[0] vs $readname\n" unless ($bamsplit[0] eq $readname); 
    #this works apparently with bwa mem output - all reads are there, no reads are left out from the bam file.
    if ($j==1) {
      $fwdmap = $bamsplit[2];
      $match = $bamsplit[5]; #parse out the CIGAR string (e.g. 120M4D20M = 120 bases match the reference, deletion of 4, 20 match the reference)
    } else {
      $revmap = $bamsplit[2];
      $match .= $bamsplit[5];
      $lastread = $readname;
      last;
    }
  }
  
  if ($fwdmap ne $revmap) {
    warn "$counter: $readname: fwd and rev read do not map to same CRS\n"; #this sometimes happens, in this case issue warning and skip. could look into what is happening here.
    next;
  }
  
  ###of the fwd and rev file how many matches are there
  $match =~ s/\*/0M/g;
  my @numbers = split /[A-Z]/, $match;
  my @chars = split /\d+/, $match;
  shift @chars;
  my $score;
  for ($j=0;$j<=$#chars;$j++) {
    $score += $numbers[$j] if ($chars[$j] eq "M");
  }
  


  #the data structure this gets stored into is:
  # HASH of BARCODES -> HASH of CRS (mappinmg to that barcode) -> ARRAY with statistics
  if (defined($BC_CRS_raw{$bcline})) {
    if (defined($BC_CRS_raw{$bcline}->{$fwdmap})) {
      #if the barcode has been seen with the same mapping, count it as one read and add mapping score to list
      push @{$BC_CRS_raw{$bcline}->{$fwdmap}->[0]}, $score;
      $BC_CRS_raw{$bcline}->{$fwdmap}->[1]++;
    } else {
        #if the barcode has been seen with a different mapping, add an entry in the hash. They will later be merged based on the read count.
        $BC_CRS_raw{$bcline}->{$fwdmap} = [ [$score], 1  ];
      }
  } else {
    #if this barcode has not been seen, create an anonymous hash
    $BC_CRS_raw{$bcline} = { $fwdmap => [ [$score], 1 ] } ; #here the data structure is set up.
  }
  
}

close(BC);
close(BAM);

#############################################
####       Clean up multiple assignments#####
#############################################

#There are some barcodes that get assigned to multiple CRS. actually this happens a lot. In these cases count how often this happens but assign the barcode to the
#CRS that is supported with more reads, buty report the number of reads supporting a different assignment. This does not remove anything 

my %BC_CRS; #this is a more simple datastructure: a hash mapping each barcode to an array with the following entries
#[assigned CRS, [mapping scores], reads supporting this assignment, reads supporting other assignments]

while ( my ($barcode, $crs) = each (%BC_CRS_raw)) {
  #sort the CRS by number of reads
  my @sorted = sort {$crs->{$b}->[1] <=> $crs->{$a}->[1] } keys %{$crs}; #remember each value in BC_CRS_raw s a (reference to) a hash with CRS names as keys
  my $current = shift @sorted;
  $BC_CRS{$barcode} = [$current, $crs->{$current}->[0], $crs->{$current}->[1], 0 ]; #main entry: the crs with most reads supporting this assignment
  while ($current = shift @sorted) {
    $BC_CRS{$barcode}->[3] += $crs->{$current}->[1]; #in field 3, count how many reads upport other assignments.
  }
  
}

#############################################
####        Run barcode correction      #####
#############################################

#if two barcodes map to the same CRS and are similar, add the less abundant barcode to the more abundant one

my %CRS_BC;
#for this, create a hash that maps each CRS to a list of barcodes
while ( my ($barcode, $crs) = each (%BC_CRS)) {
  unless (defined($CRS_BC{ $crs->[0] })) {
    $CRS_BC{ $crs->[0] } = [$barcode];
  } else {
    push @{$CRS_BC{ $crs->[0] }} , $barcode;
  }
  
}

my %BC_CRS_fixed;
while (my ($crs, $barcodes) = each(%CRS_BC)) {
  my @BARCODES = @{$barcodes}; #all barcodes mapping to a given CRS
  @BARCODES = sort { $BC_CRS{$a}->[2] <=> $BC_CRS{$b}->[2] } @BARCODES;  #sort barcodes by number of raeds
  while (my $this_bc = shift @BARCODES) {
    #starting with the barcode with the least reads
    
    my @MATCHES;
    if (@BARCODES) {
      @MATCHES = amatch($this_bc,  ["S1","I0","D0"], @BARCODES); #find approx. matches with at least one substitution in all barcodes with more reads
    }
    
    if (@MATCHES) {
      @MATCHES = sort { $BC_CRS{$b}->[2] <=> $BC_CRS{$a}->[2]} @MATCHES; #if any are found, sort matches by barcode counts
       if (defined($BC_CRS_fixed{$MATCHES[0]})) {
          #if for a certain barcode a barcode with lower counts has been mapped, add the statistics
          $BC_CRS_fixed{$MATCHES[0]}->[2] += $BC_CRS{$this_bc}->[2]; #add read counts
          $BC_CRS_fixed{$MATCHES[0]}->[3] += $BC_CRS{$this_bc}->[3]; #add deviant read counts
          push @{$BC_CRS_fixed{$MATCHES[0]}->[1]} , @{$BC_CRS{$this_bc}->[1]}; #append mapping statistics
          } else {
          $BC_CRS_fixed{$MATCHES[0]} = $BC_CRS{$this_bc}; #if not match, simply copy it over
        }
      #warn "Replace $this_bc ( $BC_CRS{$this_bc}->[2] ) by $MATCHES[0] ($BC_CRS{$MATCHES[0]}->[2] )\n";
    } else {
      if (defined($BC_CRS_fixed{$this_bc})) {
          #if for a certain barcode a barcode with lower counts has been mapped, add the statistics
          $BC_CRS_fixed{$this_bc}->[2] += $BC_CRS{$this_bc}->[2]; #add read counts
          $BC_CRS_fixed{$this_bc}->[3] += $BC_CRS{$this_bc}->[3]; #add deviant read counts
          push @{$BC_CRS_fixed{$this_bc}->[1]} , @{$BC_CRS{$this_bc}->[1]}; #append mapping statistics
          } else {
          $BC_CRS_fixed{$this_bc} = $BC_CRS{$this_bc}; #if not match, simply copy it over
        }
    }
  }
  
}


#############################################
####         Output results to csv       ####
#############################################


open(OUT, ">$outfile");
print OUT "BARCODE,CRS,READS,DEVIANTREADS,MEANMATCHES,MINMATCHES,MAXMATCHES\n";
while ( my ($barcode, $crs) = each (%BC_CRS_fixed)) {
  my $line = join(",", ($barcode , $crs->[0], $crs->[2],$crs->[3], sum(@{$crs->[1]}) / scalar(@{$crs->[1]}), min(@{$crs->[1]}),max(@{$crs->[1]})    ));
  print OUT "$line\n";
  
}
close(OUT);