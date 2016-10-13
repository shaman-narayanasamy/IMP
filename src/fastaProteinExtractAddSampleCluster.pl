#!/usr/bin/perl

use Bio::DB::Fasta;

my $fastaFile = shift;
my $queryFile = shift; #0: gene, 1: sample, 2: cluster

my $db = Bio::DB::Fasta->new( $fastaFile );

open (IN, $queryFile);
while (<IN>){
    chomp;
    my @fields = split /\t/;
    my $seq = $fields[0];
    unless ($seq eq "gene") {
        my $sequence = $db->seq($seq);
        if  (!defined( $sequence )) {
            print STDERR "Sequence $seq not found. \n";
            next;
        }
        my $lib = $fields[1];
        my $cluster = $fields[2];
        print ">$lib","_","$cluster","_","$seq\n", "$sequence\n";
    }
}
