#!/usr/bin/perl -w
use strict;

open(I,"<$ARGV[0]");
while (<I>) {
	chomp;
	if (/^>(.*)/) { print ">$ARGV[1]\_$1\n" } else { print "$_\n" }
}
close I;
