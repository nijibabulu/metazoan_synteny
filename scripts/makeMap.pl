#!/usr/bin/perl -w
use strict;
if (@ARGV!=4) {
 die("Usage: $0 GFF_FILE PREFIX FEATURE ATTR\n".
     "\nPREFIX should be added to the beginning of gene names\n".
     "  when the gene names are not already unique. A '.' indicates\n".
     "  no prefix should be added.\n".
     "FEATURE is a feature of the GFF file, e.g. CDS or mRNA\n".
     "ATTR is the gff attribute of the GFF file in which genes are\n".
     "  organized, e.g. gene_id or transcript_id\n")
}

my $pre=$ARGV[1];
my $pre2="$ARGV[1]\_";
if ($ARGV[1] eq ".") {
 $pre="";
 $pre2="";
}
my $key=$ARGV[2];
my $name=$ARGV[3];

open(I,"<$ARGV[0]");
my %chr=();
while (<I>) {
	chomp;
	my @tmp = split /\t/;
	if (/^#/ || /^$/) { next }
	if ($tmp[2] eq $key) {
		my $id="";
		if ($tmp[$#tmp]=~/$name=([^\;]*)/) {
			$id=$1;
		}
		if ($tmp[$#tmp]=~/$name \"([^\;\"]*)\"/) {
                        $id=$1;
                }
		if ($id eq "") { next }
		push @{$chr{$id}}, [ ($tmp[0],$tmp[6],$tmp[3],$tmp[4]) ];
	}
}
close I;
print STDERR " ".keys(%chr)." genes found!\n";

for my $x (keys %chr) {
 my @ar= @{$chr{$x}};
 @ar=sort {${$a}[2]<=>${$b}[2]} @ar;
 print "$pre\t$pre2$x\t$ar[0][0]\t$ar[0][1]\t$ar[0][2]\t$ar[$#ar][3]\n";
}
