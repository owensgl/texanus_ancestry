#/bin/perl
use strict;
use warnings;

#This takes a DA file (Derived-Ancestral), a sample label file, and a population file. The population file should have groups 1, 2 and 3. It outputs for each site the ABBA, BABA, and Fd values for each site for windowed plotting.

my $samplefile = $ARGV[0]; #List of samples with pop label
my $groupfile = $ARGV[1]; #List of populations that represent groups.


my %pop;
open SAMPLE, $samplefile;
while(<SAMPLE>){
  chomp;
  my @a = split(/\t/,$_);
  $pop{$a[0]} = $a[1];
}
close SAMPLE;
my %group;
open GROUP, $groupfile;
while(<GROUP>){
  chomp;
  my @a = split(/\t/,$_);
  $group{$a[0]} = $a[1];
}

my %sample;
my @group_1;
my @group_2;
my @group_3;
my %hash;
my %cols;
my $counter = 1;

print "chr\tpos\tabba\tbaba\tp3_abba\tp3_baba";
while(<STDIN>){
  my %derived_count;
  my %total_count;
  chomp;
  my @a = split(/\t/,$_);
  $counter++;
  if ($counter % 10000 == 0){ print STDERR "Processing $a[0] $a[1]\n";}
  if ($. == 1){
    foreach my $i (2..$#a){
      $sample{$i} = $a[$i];
      unless($pop{$a[$i]}){next;}
      unless($group{$pop{$a[$i]}}){next;}
      if ($group{$pop{$a[$i]}} == 1){
        $cols{$i}= 1;
      }elsif ($group{$pop{$a[$i]}} == 2){
        $cols{$i} = 2;
      }elsif ($group{$pop{$a[$i]}} == 3){
        $cols{$i} = 3;
      }
    }
  }else{
    foreach my $i (2..$#a){
      if ($cols{$i}){
        if ($a[$i] eq "A"){
          $total_count{$cols{$i}}+=2;
        }elsif ($a[$i] eq "H"){
	  $total_count{$cols{$i}}+=2;
	  $derived_count{$cols{$i}}+=1;
	}elsif ($a[$i] eq "B"){
  	  $total_count{$cols{$i}}+=2;
          $derived_count{$cols{$i}}+=2;
	}
      }
    }
    my %freq;
    foreach my $i (1..3){
      unless($total_count{$i}){goto SKIPLINE;}
      unless($derived_count{$i}){
        $derived_count{$i} = 0;
      }
      $freq{$i} = $derived_count{$i}/$total_count{$i};
    }
    my $abba = (1- $freq{1}) * $freq{2} * $freq{3};
    my $baba = $freq{1} * (1 - $freq{2}) * $freq{3};
    my $p3_abba;
    my $p3_baba;
    if ($freq{3} > $freq{2}){
      $p3_abba = (1- $freq{1}) * $freq{3} * $freq{3};
      $p3_baba = $freq{1} * (1 - $freq{3}) * $freq{3};
    }else{
      $p3_abba = (1- $freq{1}) * $freq{2} * $freq{2};
      $p3_baba = $freq{1} * (1 - $freq{2}) * $freq{2};
    }
    if ($abba == 0 and $baba == 0){next;}
    print "\n$a[0]\t$a[1]\t$abba\t$baba\t$p3_abba\t$p3_baba";
    SKIPLINE:
  }
}
