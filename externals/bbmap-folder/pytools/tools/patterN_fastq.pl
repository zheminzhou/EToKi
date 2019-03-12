#!/usr/bin/env perl
# counts N bases occuring in same spot on a read to identify stutters in illumina runs
# by packing bases into a bitstring and using as hash key. 
# TODO: these have been observed occuring tilewise most often but also across lanes so
#       counting by tile rather than by all input is probably useful.


use strict;
use warnings;
use Getopt::Long;

my $min_count  = 1;
my $min_pct    = 1 / 120;  # based on 120 tiles to a lane
my $length     = 0;
my $i          = 0;
my %hash       = ();
my $tile_parse = 0;
my $analog     = 0;
my ( $in_file, $out_file, $help );

GetOptions( 
  'count=i' => \$min_count,
  'pct=f'   => \$min_pct,
  'in=s'    => \$in_file,
  'help'    => \$help,
  'tile'    => \$tile_parse,
  'analog'  => \$analog,
  'out=s'   => \$out_file,
) or die $!;

if( $help || ! defined $in_file ) {
  print STDERR "$0 -count MIN_OCCURENCES -pct MIN_PCT_READS_WITH_PATTERN -in s_1_1_sequence.fastq [ -out out.stats ][ -tile ][ -analog ]\n";
  print STDERR "\n";
  print STDERR "  -count  INT   min. occurences of a pattern before it's reported\t[ $min_count ]\n";
  print STDERR "  -pct    FLOAT min. % of reads with a particular pattern before it's reported\t[ 1 / 120 ]\n";
  print STDERR "  -tile         evaluate pct per tile\n";
  print STDERR "  -analog       output like -N--- instead of the hipper 01000\n";
  print STDERR "\n";
  exit;
}

open( my $fh, $in_file ) or die "ERROR: failed to open $in_file for reading $!";

my $out_fh;
if ( $out_file ) {
  open( $out_fh, '>', $out_file ) or die $!;
}
else {
  $out_fh = *STDOUT;
}

my %tile_counts = ();

while ( <$fh> ) {              # readname 
 my $tile = $tile_parse ? ( split /:/,$_,4 )[2] : '';
 $i++;
 if( $tile_parse ) { $tile_counts{ $tile }++; }
 $_=<$fh>;                     # sequence
 if ( /N/ ) {
   chomp;
   tr/ACGTN/00001/;
   if ( ! $length ) {
     $length = length;
   }
   if ( $tile_parse ) {
     $hash{ $tile }{ pack( "b$length", $_ ) }++;
   }
   else {
     $hash{ pack( "b$length", $_ ) }++;
   }
 }
 <$fh>;                        # seperator of seq/qual
 <$fh>;                        # quality
}
my $sum_of_bad_past_threshold_pct = 0;
my $sum_of_bad_past_threshold     = 0;
my %tile_sums = ();

if ( $tile_parse ) {       #key is a tile, then maps as sub keys.
  for my $key ( sort { $a <=> $b } keys %hash ) {
    for my $pattern ( sort { $hash{ $key }{ $a } <=> $hash{ $key }{ $b } } keys %{ $hash{ $key } } ) { 
      if ( $hash{ $key }->{$pattern} >= $min_count && $hash{ $key }->{$pattern} / $tile_counts{ $key } * 100 >= $min_pct  ) { 
        my $print_pattern = unpack( "b$length", $pattern );
        if ( $analog ) { $print_pattern =~ tr/01/-N/; }
        print $out_fh "$key\t$hash{ $key }{ $pattern }\t$print_pattern\n";
        $sum_of_bad_past_threshold_pct += $hash{ $key }->{$pattern} / $i * 100;
        $sum_of_bad_past_threshold     += $hash{ $key }->{$pattern};
        $tile_sums{ $key }         += $hash{ $key }->{$pattern} / $tile_counts{ $key } * 100;
      }
    }
  }
}
else {
  for my $key ( keys %hash ) {
    if ( $hash{ $key } >= $min_count && $hash{ $key } / $i * 100 >= $min_pct ) { 
      my $print_pattern = unpack( "b$length", $key );
      if ( $analog ) { $print_pattern =~ tr/01/-N/; }
      print $out_fh $hash{$key} , "\t$print_pattern\n";
      $sum_of_bad_past_threshold_pct += $hash{ $key } / $i * 100;
      $sum_of_bad_past_threshold += $hash{ $key };
    }
  }
}

# summary stuff

if ( $sum_of_bad_past_threshold_pct ) {
  print $out_fh "sum pct patterNs past $min_pct == $sum_of_bad_past_threshold_pct ( $sum_of_bad_past_threshold / $i * 100 )\n";
}

if ( $tile_parse && $sum_of_bad_past_threshold_pct ) {
  for my $tile ( sort { $tile_sums{ $a } <=> $tile_sums{ $b } } keys %tile_sums ) {
    print $out_fh "tile $tile pct patterNs past $min_pct == $tile_sums{ $tile }\n";
  }
}
