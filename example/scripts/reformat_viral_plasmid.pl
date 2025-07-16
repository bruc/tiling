#!/usr/bin/perl

use strict;
use warnings;

use Bio::Frescobi::Sequtil;
use Bio::Frescobi::Genutil;

use Getopt::Long;

my $itr_file = "flip.fa";
my $minlength = 100;
my $file = "";

GetOptions("file=s" => \$file,
	   "itr_file=s" => \$itr_file,
	   "minlength=i" => \$minlength);

my $usage = <<EOF;
reformat_viral_plasmid.pl -file=<file>
                        [ -itr_file=<file>     ]
                        [ -minlength=<integer> ]
EOF

if (scalar(@ARGV) != 0) {
    print STDERR $usage;
    exit(1);
}

open (SEQ, "<$file") ||
    die "Unable to open $file: $!\n";
$_ = <SEQ>;
my ($name, $annotation, $seq) = read_1_fasta_sequence(\*SEQ, $_, 0);
close SEQ;

my $cmd = "blastall -i $file -d $itr_file -p blastn -m 9 -F F -G 1";
open (PIPE, "$cmd |" ) ||
    die "Unable to open pipe command $cmd: $!\n";

my @itr_hits;

while (<PIPE>) {
    next if m/^#/;
    chomp;
    my ($query_id, $subject_id, $pct_identity, $length, $mismatches,
	$gaps, $q_start, $q_end, $s_start, $s_end, $e, $bits) = split(/\t/);
    if ($q_start > $q_end) {
	print STDERR "Bad blastn result (q_start > q_end): $_\n";
	exit(1);
    }
    if ($length >= $minlength) {
	push (@itr_hits,
	      [ $q_start, $q_end ]);
    }
}

@itr_hits = sort {$a->[0] <=> $b->[0]} @itr_hits;

if (scalar(@itr_hits) != 2) {
    print STDERR "$file does not have exactly two ITR hits.\n";
    print STDERR "No splitting will be done and no ITR added.\n";
    if ($name eq "pHM-00224") {
	$name = "RepCap $name";
    }
    elsif ($name eq "pHMI-00023" or $name eq "pOX-00023") {
	$name = "Helper $name";
    }
    print ">$name $annotation\n";
    print format_seq($seq);
}
else {
    my $payload_len = $itr_hits[1]->[0] - $itr_hits[0]->[1] - 1;
    my $payload = substr($seq, $itr_hits[0]->[1], $payload_len);
    my $backbone = substr($seq, $itr_hits[1]->[1]) . substr($seq, 0, $itr_hits[0]->[0] - 1);
    my $itr_seq = `cat $itr_file`;

    print $itr_seq;
    print ">Payload $name $annotation\n";
    print format_seq($payload);
    print ">Backbone $name $annotation\n";
    print format_seq($backbone);
}
