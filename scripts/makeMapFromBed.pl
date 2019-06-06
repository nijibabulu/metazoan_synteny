#!/usr/bin/perl -w
use strict;
if (@ARGV!=2) {
 die(" Usage: BED PREFIX\n\n".
     "BED should be a standard bed12 file\n".
      "PREFIX should be some identifier to add to the beginning fo IDs\n".
      "  if they are not already unique to the study. Use '.' for no prefix.");
}

my $pre=$ARGV[1];
my $pre2="$ARGV[1]\_";
if ($ARGV[1] eq ".") {
 $pre="";
 $pre2="";
}
open(I,"<$ARGV[0]");
my %chr=();
while (<I>) {
	chomp;
	my @tmp = split /\t/;
	if (/^#/) { next }
   $chr{$tmp[3]} = [$tmp[0],$tmp[5],$tmp[1],$tmp[2]] ;
}
close I;
print STDERR " ".keys(%chr)." genes found!\n";

for my $x (keys %chr) {
 my @ar= @{$chr{$x}};
 print "$pre\t$pre2$x\t$ar[0]\t$ar[1]\t$ar[2]\t$ar[3]\n";
}
