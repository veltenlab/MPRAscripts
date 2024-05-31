use strict;
use warnings;
use String::Approx 'amatch';
use List::Util qw(min max sum);

my ($bc, $bam, $outfile) = @ARGV;

my $time = CORE::localtime;
warn "Started the job at: $time\n";

#################################################
####      Mapping barcodes and alignments    ####
#################################################

open(BC, "zcat $bc|");
open(BAM, "samtools view $bam| "); my $bamline;


sub get_next_barcode {
  my $i = 0;
  my $readname;
  my $barcode;
  my $bcline;
  for (my $i=1;$i<=4;$i++) {
    return undef unless ($bcline = <BC>);
    chomp $bcline;
    $readname = $bcline if $i == 1;
    $barcode = $bcline if $i == 2;
  }
  
  $readname =~ s/^@(\S+)\s.+/$1/;  
  my @OUT = ($readname, $barcode);
  return \@OUT;
}

sub cigar_to_score {
    ###of the fwd and rev file how many matches are there
  my $match = shift @_;
  $match =~ s/\*/0M/g;
  my @numbers = split /[A-Z]/, $match;
  my @chars = split /\d+/, $match;
  shift @chars;
  my $score;
  for (my $j=0;$j<=$#chars;$j++) {
    $score += $numbers[$j] if ($chars[$j] eq "M");
  }
  return $score;
}

  my %BC_CRS_raw; #hash mapping barcode to CRS
  my $counter = 0;
  my @bamsplit; #split bam entries
  my %maps_score;
  my %maps_count_fwd;
  my %maps_count_rev;
  my $multiple_alignments = 0;
  ### read in all lines in the BAM file that correspond to that read name
  my $nbc = get_next_barcode();
  my ($readname, $bcline) = @{$nbc};
  my $score;
  my $selectedmap;
  
  while ($bamline = <BAM>) {
    
    
    @bamsplit = split /\t/, $bamline;
    my @bits = reverse(split(//, sprintf("%b", $bamsplit[1]))); #parse out flags
    
    if ($bamsplit[0] ne $readname) { #if the bam file is ahead
      #wrap up
      #select only alignments where both fwd and rev read match

      my @considered = grep {$maps_count_fwd{$_} == 1 and $maps_count_rev{$_} == 1} keys %maps_score;
      
      
      if (scalar(@considered) == 1) {
        $selectedmap = $considered[0];
        $score = $maps_score{$considered[0]};
      } elsif(scalar(@considered)  > 1) {
        my @sorted = sort {$maps_score{$b} <=> $maps_score{$a} } @considered;
        $selectedmap = $sorted[0];
        $score = $maps_score{$sorted[0]};
        #warn "$counter: $readname: multimapping\n"; # CRS to " . scalar(@considered) . " candidates, best: $sorted[0] , second: $sorted[1]\n";
        
      }
      
      if (@considered and $selectedmap ne "*") {
        if (defined($BC_CRS_raw{$bcline})) {
        if (defined($BC_CRS_raw{$bcline}->{$selectedmap})) {
          #if the barcode has been seen with the same mapping, count it as one read and add mapping score to list
          push @{$BC_CRS_raw{$bcline}->{$selectedmap}->[0]}, $score;
          $BC_CRS_raw{$bcline}->{$selectedmap}->[1]++;
        } else {
          #if the barcode has been seen with a different mapping, add an entry in the hash. They will later be merged based on the read count.
          $BC_CRS_raw{$bcline}->{$selectedmap} = [ [$score], 1  ];
        }
      } else {
        #if this barcode has not been seen, create an anonymous hash
        $BC_CRS_raw{$bcline} = { $selectedmap => [ [$score], 1 ] } ; #here the data structure is set up.
      }
      
      } #else {
        #warn "$counter: $readname: not mapping\n";
      #}
      
      
      #Start next round get next barcode
      %maps_score = ();
      %maps_count_fwd = ();
      %maps_count_rev = ();
      $nbc = get_next_barcode();
      last unless $nbc;
      ($readname, $bcline) = @{$nbc};
      $counter++;
      die "Something is wrong" unless $bamsplit[0] eq $readname;
    }
    
    
    if (defined($maps_score{$bamsplit[2]})) {
      $maps_score{$bamsplit[2]} += cigar_to_score($bamsplit[5]);
      if ($bits[4]) {
        $maps_count_fwd{$bamsplit[2]}++;
      } else {
        $maps_count_rev{$bamsplit[2]}++;
      }
    } else {
      $maps_score{$bamsplit[2]} = cigar_to_score($bamsplit[5]);
      $maps_count_fwd{$bamsplit[2]} = $bits[4];
      $maps_count_rev{$bamsplit[2]} = 1 - $bits[4];
    }
    
  }
  
  #the data structure this gets stored into is:
  # HASH of BARCODES -> HASH of CRS (mappinmg to that barcode) -> ARRAY with statistics


close(BC);
close(BAM);

$time = CORE::localtime;
warn "$time: Completed reading BAM file\n";

#############################################
####       Clean up multiple assignments#####
#############################################

#There are some barcodes that get assigned to multiple CRS. actually this happens a lot. In these cases count how often this happens but assign the barcode to the
#CRS that is supported with more reads, buty report the number of reads supporting a different assignment. This does not remove anything 

my %BC_CRS; #this is a more simple datastructure: a hash mapping each barcode to an array with the following entries
#[assigned CRS, [mapping scores], reads supporting this assignment, reads supporting other assignments]
my $total_reads;

while ( my ($barcode, $crs) = each (%BC_CRS_raw)) {
  #sort the CRS by number of reads
  $total_reads += sum(map {$_->[1]} values %{$crs});
  my @sorted = sort {$crs->{$b}->[1] <=> $crs->{$a}->[1] } keys %{$crs}; #remember each value in BC_CRS_raw s a (reference to) a hash with CRS names as keys
  my $current = shift @sorted;
  $BC_CRS{$barcode} = [$current, $crs->{$current}->[0], $crs->{$current}->[1], 0 ]; #main entry: the crs with most reads supporting this assignment
  while ($current = shift @sorted) {
    $BC_CRS{$barcode}->[3] += $crs->{$current}->[1]; #in field 3, count how many reads upport other assignments.
  }
  
}

undef %BC_CRS_raw;

$time = CORE::localtime;
warn "$time: Completed creating BC_CRS\n";


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

$time = CORE::localtime;
warn "$time: Completed creating CRS_BC\n";

my %BC_CRS_fixed;
my $nstep;
my $ncrs = scalar keys %CRS_BC;
my $removed_oneread = 0;
my $total_mapped_reads =0;

while (my ($crs, $barcodes) = each(%CRS_BC)) {

  ######
  $nstep++;
  my @BARCODES = @{$barcodes}; #all barcodes mapping to a given CRS
#print information
  my $totreads = sum(map { $BC_CRS{$_}->[2] } @BARCODES);
  $total_mapped_reads+= $totreads;
  $time = CORE::localtime;
  #warn "$time: Correcting barcodes for CRS $crs with " .  scalar(@BARCODES) . " barcodes and $totreads reads, $nstep of $ncrs\n";

#filter
  if (scalar(@BARCODES) > 1000) { #for all CRS with more than 1000 non error corrected barcodes
  $removed_oneread += scalar( grep {$BC_CRS{$_}->[2] == 1} @BARCODES );
  @BARCODES = grep {$BC_CRS{$_}->[2] > 1} @BARCODES; #simply through out all barcodes with 1 read but count how many these are
  #warn "for $crs removed barcodes with less than 2 reads\n";
  }

#SortBC
  @BARCODES = sort { $BC_CRS{$a}->[2] <=> $BC_CRS{$b}->[2] } @BARCODES;  #sort barcodes by number of reads
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

undef %CRS_BC;
undef %BC_CRS;

$time = CORE::localtime;
warn "$time: Completed BC correction\n";


#############################################
####         Output results to csv       ####
#############################################


open(OUT, "| gzip -c > $outfile") or die "Can't write to file '$outfile' [$!]\n";
print OUT "BARCODE,CRS,READS,DEVIANTREADS,MEANMATCHES,MINMATCHES,MAXMATCHES\n";
while ( my ($barcode, $crs) = each (%BC_CRS_fixed)) {
  my $line = join(",", ($barcode , $crs->[0], $crs->[2],$crs->[3], sum(@{$crs->[1]}) / scalar(@{$crs->[1]}), min(@{$crs->[1]}),max(@{$crs->[1]})    ));
  print OUT "$line\n";
}
close(OUT);

$time = CORE::localtime;

warn "$total_reads reads seen in BC_CRS_raw.\n";
warn "$total_mapped_reads reads seen in CRS_BC.\n";
warn "$removed_oneread reads removed - barcode covered only once.\n";
warn "Finished the job at: $time\n";
