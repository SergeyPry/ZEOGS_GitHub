#!c:/Perl/bin/perl.exe -w

use strict;

open READ, "<staged_anatomy.txt" or die $!;

undef $/;
my @stage_terms = split(/\n\n/, <READ>);

my $count = 1;


foreach(@stage_terms){
	my @terms = split(/\n/, $_);
	open WRITE, ">stage$count.txt" or die $!;
	
	foreach(@terms[1..$#terms]){
		my $term = trim($_);
		print WRITE "$term\n";
	}
	$count++;
	
	close WRITE;

}

sub trim {
	my $string = shift;
	$string =~ s/^\s+//;
	$string =~ s/\s+$//;
	return $string;
}