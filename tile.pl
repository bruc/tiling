#!/usr/bin/perl

# Tiling prototype.

# Author: Robert Bruccoleri, Congenomics, LLC.

use warnings;
use strict;
use English '-no_match_vars';
use POSIX;

use Data::Dumper;
use Getopt::Long;
use Bio::Frescobi::Genutil;
use Bio::Frescobi::CgPg;
use Bio::Frescobi::Sequtil;
use Bio::Frescobi::BigQuery;
use Bio::Frescobi::Dbseq;

my $gap_limit = 30;
my $max_iter = 100000;
my $dbname = 'itrext';
my $driver = 'Pg';
my $library = "";
my $consolidate = 0;
my $db = "";
my @arg_other_dbs;
my $limit = 0;
my $min_gap_seq = 6;
my $min_gap_search = 11;
my $min_poly = 50;
my $by_strand = 0;
my $end_gap = 0;
my $write_seq_limit = 0;
my $debug = 0;
my $raw = 1;
my $fastq = "";
my $by_component = 0;
my $add_gaps = 0;
my $quant_homopolymers = 1;
my $details = 0;
my $add_unknowns = 0;

GetOptions("gap_limit=i" => \$gap_limit,
	   "max_iter=i" => \$max_iter,
	   "dbname=s" => \$dbname,
	   "driver=s" => \$driver,
	   "library=s" => \$library,
	   "consolidate!" => \$consolidate,
	   "db=s" => \$db,
	   "other_db=s" => \@arg_other_dbs,
	   "limit=i" => \$limit,
	   "min_gap_seq=i" => \$min_gap_seq,
	   "min_gap_search=i" => \$min_gap_search,
	   "min_poly=i" => \$min_poly,
	   "by_strand!" => \$by_strand,
	   "end_gap=i" => \$end_gap,
	   "write_seq_limit=i" => \$write_seq_limit,
	   "raw!" => \$raw,
	   "debug!" => \$debug,
	   "fastq=s" => \$fastq,
	   "by_component!" => \$by_component,
	   "add_gaps!" => \$add_gaps,
	   "quant_homopolymers!" => \$quant_homopolymers,
	   "details!" => \$details,
	   "add_unknowns!" => \$add_unknowns);

my $usage = <<EOF;
tile.pl -library=<string>
        -db=<string>
      [ -gap_limit=<integer>      ]
      [ -max_iter=<integer>       ]
      [ -dbname=<string>          ]
      [ -driver=<string>          ]
      [ -[no]consolidate          ]
      [ -other_db=<string>...     ]
      [ -min_gap_seq=<integer>    ]
      [ -min_gap_search=<integer> ]
      [ -limit=<integer>          ]
      [ -[no]by_strand            ]
      [ -end_gap=<integer>        ]
      [ -write_seq_limit=<integer>]
      [ -[no]raw                  ]
      [ -fastq=<file>             ]
      [ -[no]by_component         ]
      [ -[no]add_gaps             ]
      [ -[no]quant_homopolymers   ]
      [ -[no]details              ]
EOF
    ;

my %other_dbs;

if ($library eq "") {
    print STDERR "No library specified.\n";
    print STDERR $usage;
    exit(1);
}

if ($db eq "") {
    print STDERR "No template database specified.\n";
    print STDERR $usage;
    exit(1);
}

if (scalar(@arg_other_dbs) == 0) {
    @arg_other_dbs = ("Helper", "RepCap");
}

foreach my $other_db (@arg_other_dbs) {
    $other_dbs{$other_db} = 1;
}

if ($by_component and $fastq eq "") {
    print STDERR "The by_component option requires FASTQ input.\n";
    exit(1);
}

if (scalar(@ARGV) != 0) {
    print STDERR $usage;
    exit(1);
}

if ($min_poly < 5 or $min_poly >= 32766) {
    print STDERR "min_poly ($min_poly) must be 5 or greater and less than 32766.\n";
    exit(1);
}

$min_poly -= 1;

my $pg = Bio::Frescobi::CgPg::new(dbname => $dbname,
				  driver => $driver,
				  cgi => 0);

if ($debug) {
    $pg->default_echo(1);
}

open (STDOUT, ">${library}.tiles") ||
    die "Unable to open for writing ${library}.tiles: $!\n";

open (GAPSEQ, ">${library}.gap.seq") ||
    die "Unable to open for writing ${library}.gap.seq: $!\n";

open (UNMATCHED, ">${library}.unmatched.seq") ||
    die "Unable to open for writing ${library}.unmatched.seq: $!\n";

open (SUMMARY, ">${library}.summary") ||
    die "Unable to open for writing ${library}.summary: $!\n";

if ($by_component) {
    open (SPLIT_READS, ">${library}.split.reads") ||
	die "Unable to open for writing ${library}.split.reads: $!\n";
}

if ($details) {
    open (HSPS, ">${library}.hsps") ||
	die "Unable to open ${library}.hsps: $!\n";
    print HSPS join("\t",
		    "Seqid",
		    "Status",
		    "Seq Length",
		    "Reference DB",
		    "Match Name",
		    "Match Length",
		    "Fraction Identical",
		    "Query Start",
		    "Query Stop",
		    "Subject Strand",
		    "Subject Start",
		    "Subject Stop",
		    "Used in Tile?",
		    "Number Identical",
		    "Gaps"), "\n";
    open (TILEDETAILS, ">${library}.tile.details") ||
	die "Unable to open ${library}.tile.details: $!\n";
    print TILEDETAILS join("\t",
			   "Seqid",
			   "Occurrence Count",
			   "Status",
			   "Tile Code"), "\n";
    open (NAME2SEQID, ">${library}.name2seqid") ||
	die "Unable to open ${library}.name2seqid: $!\n";
    print NAME2SEQID join("\t",
			 "Seqid",
			 "Name"), "\n";
}

my %seqid_count = $pg->get_hash_for_field("select seqid, count(*) " .
					  "  from raw_seqs " .
					  " where library = " . quotify($library) .
					  " group by seqid");
if (scalar(keys %seqid_count) == 0) {
    print STDERR "library $library appears to have no sequences.\n";
    exit(1);
}

my %names;
if ($by_strand) {
    my $string_agg_function;
    if ($driver eq 'Pg') {
	$string_agg_function = 'string_agg';
    }
    else {
	$string_agg_function = 'group_concat';
    }
    %names = $pg->get_hash_for_field("select seqid, ${string_agg_function}(name, ':') " .
				     "  from raw_seqs " .
				     " where library = " . quotify($library) .
				     " group by seqid")
}

my %fastq_data;
if ($fastq ne "") {
    open (FASTQ, "<$fastq") ||
	die "Unable to open $fastq: $!\n";
    while (my ($name, $annotation, $seq, $qual) =
	   read_1_fastq_sequence(\*FASTQ, 'b')) {
	if (exists($fastq_data{$name})) {
	    die "Duplicate name $name found in $fastq\n";
	}
	$fastq_data{$name} = [ $annotation, $seq, $qual ];
    }
}

my %feat_counts;
my %code_counts;
my %code_seqids;
my %unidentified_seqs;
my %tile_failure_seqs;
my %component_files;

my $count = 0;
my $big_gap_count = 0;
my $error_gap_count = 0;
my $recursion_error_count = 0;
my $unknown_tiling_error_count = 0;
my $tiles_with_unknown_count = 0;
my $complete_tiling_count = 0;
my $partial_tiling_count = 0;
my $split_references_count = 0;
my $total_unknown_count = 0;

if ($pg->get_single_value("select count(*) from hits where lower(key) = 'unknown'") > 0 and
    $add_unknowns) {
    print STDERR "Unable to add unknowns to the tiling patterns. At least one references has a name of 'unknown'\n";
    exit(1);
}

my $db_clause = " and hits.common_db_name in ( " .
    quotify($db) .
    ", " .
    join(",", map {quotify($_)} keys %other_dbs) .
    ")";
my $limit_clause = $limit > 0 ? "limit $limit" : "";
my $round_frac;
if ($driver eq 'Pg') {
    $round_frac = 'round(hsps.frac_identical::numeric, 4)';
}
else {
    $round_frac = 'round(hsps.frac_identical, 4)';
}

my @hsps;
my @other_hsps;
my %seqids_with_hsps;
my %codes_for_zmw;

my $query = Bio::Frescobi::BigQuery->new($pg, 
					 "select hits.seqid, " . # 0
					 "       r.library, " .
					 "       d.length, " . # 2
					 "       hits.common_db_name as db, " .
					 "       hsps.bits, " . # 4
					 "       hsps.expect, " .
					 "       hsps.length, " . # 6
					 "       hsps.gaps, " .
					 "       ${round_frac}, " . # 8
					 "       hsps.query_start as q_start, " .
					 "       hsps.query_end as q_end, " . # 10
					 "       hsps.query_strand_positive as qs, " .
					 "       hsps.sbjct_start as s_start, " . # 12
					 "       hsps.sbjct_end as s_end, " .
					 "       hsps.sbjct_strand_positive as ss, " . # 14
					 "       hits.key, " .
					 "       hsps.num_identical as nid " .
					 "  from hits, " .
					 "       hsps, " .
					 "       seq_data d, " .
					 "       (select distinct seqid, library from raw_seqs where library = '$library' $limit_clause) r " .
					 " where hits.hit_id = hsps.hit_id " .
					 "   and hits.seqid = d.seqid " .
					 "   and hits.seqid = r.seqid " .
					 "       $db_clause " .
					 " order by hits.seqid");
my $prev_seqid = "";
my $seqid = "";
while (my $row = $query->next) {
    my ($library, $seq_length, $db, $bits,
	$e, $length, $gaps, $frac_id, $q_start,
	$q_stop, $qs, $s_start, $s_stop, $ss, $key,
	$n_id);
    ($seqid, $library, $seq_length, $db, $bits,
     $e, $length, $gaps, $frac_id, $q_start,
     $q_stop, $qs, $s_start, $s_stop, $ss, $key,
     $n_id) = @{$row};
    # DBI and DBD::Pg can handle Booleans differently. In one case,
    # "true" is "t" and the other, "true" is 1. The query strand must
    # always be "true", so handle both cases, but otherwise, it's a
    # fatal error.
    if ($qs eq 't') { }
    elsif ($qs eq "1") {
	$qs = logical2bool($qs);
	$ss = logical2bool($ss);
    }
    else {
	print STDERR "Unexpected query strand: ", join(" ", @{$row}), "\n";
	exit(1);
    }
    if ($ss eq 't') {
	$ss = '+';
    }
    elsif ($ss eq 'f') {
	$ss = '-';
    }
    if ($key eq "ITR-FLOP") {
	# Switch it to flip.
	$key = "ITR-FLIP";
	my $new_stop = 165 - $s_start + 1;
	my $new_start = 165 - $s_stop + 1;
	if ($ss eq '+') {
	    $ss = '-';
	}
	elsif ($ss eq '-') {
	    $ss = '+';
	}
	else {
	    die "Unexpected ss: $ss\n";
	}
	$s_start = $new_start;
	$s_stop = $new_stop;
    }
    if ($seqid ne $prev_seqid) {
	if (scalar(@hsps) != 0 or scalar(@other_hsps) != 0) {
	    &process_seqid($prev_seqid);
	}
	$prev_seqid = $seqid;
	@hsps = ();
	@other_hsps = ();
    }
    my $hsp = { "seq_length" => $seq_length,
		"db" => $db,
		"key" => sprintf("%s[%d-%d]", $key, $s_start, $s_stop),
		"length" => $length,
		"frac_id" => $frac_id,
		"q_start" => $q_start,
		"q_stop" => $q_stop,
		"s_start" => $s_start,
		"s_stop" => $s_stop,
		"ss" => $ss,
		"n_id" => $n_id,
		"gaps" => $gaps,
		"used" => 0};
    if (exists($other_dbs{$db})) {
	push (@other_hsps, $hsp);
    }
    else {
	push (@hsps, $hsp);
    }
}

if (scalar(@hsps) != 0 or scalar(@other_hsps) != 0) {
    die "Assertion violated: seqid ne prev_seqid: seqid = $seqid  prev_seqid = $prev_seqid\n"
	if $prev_seqid ne $seqid;
    &process_seqid($prev_seqid);
}

foreach my $seqid (keys %seqid_count) {
    if (not exists($seqids_with_hsps{$seqid})) {
	print "No matches for $seqid.\n";
	$total_unknown_count += $seqid_count{$seqid};
	if ($details) {
	    print TILEDETAILS join("\t",
				   $seqid,
				   $seqid_count{$seqid},
				   "Completely unknown",
				   ""), "\n";
	}
    }
}

sub process_seqid {
    my $seqid = shift;
    $seqids_with_hsps{$seqid} = 1;
    printf "Processing $seqid with occurrence count of %d\n", $seqid_count{$seqid};

    my $status;
    my ($seq) = retrieve_sequence_data($pg, $seqid);
    my @filtered_hsps;
    my $coverage;

    # if there are no HSP's from the primary reference, then just look at
    # the other dbs. If there are HSP's from the primary reference, then
    # we need to look at the primaries first, and if needed, add the others.
    
    if (scalar(@hsps) == 0) {
	@hsps = @other_hsps;
	@other_hsps = ();
	&add_homopolymers($seq, $seqid);
	@filtered_hsps = &filter_hsps($seqid,
				      "other",
				      @hsps);
	$coverage = &compute_coverage(@filtered_hsps);
    }
    else {
	&add_homopolymers($seq, $seqid);
	@filtered_hsps = &filter_hsps($seqid,
				      "primary",
				      @hsps);
	$coverage = &compute_coverage(@filtered_hsps);

	if ($coverage =~ m/X\s{${min_gap_search},}X/) { # Look for an internal gap in coverage
	    # add additional HSP's from the other DB's.
	    my $offset = 0;
	    my $added = 0;
	    while ((my $start = index($coverage, ' ', $offset)) > -1) {
		my $stop = index($coverage, 'X', $start + 1);
		if ($stop == -1) {
		    $stop = length($coverage) - 1;
		}
		last if $stop >= length($coverage) - 1;
		$offset = $stop + 1;
		my $len = $stop - $start + 1;
		if ($len >= $min_gap_search) {
		    for (my $ihsp = 0; $ihsp < scalar(@other_hsps); $ihsp++) {
			my $hsp = $other_hsps[$ihsp];
			next if not defined($hsp);
			if (range_overlap($start + 1,
					  $stop + 1,
					  $hsp->{"q_start"},
					  $hsp->{"q_stop"})) {
			    push(@hsps, $hsp);
			    $added = 1;
			    $other_hsps[$ihsp] = undef;
			}
		    }
		}
	    }
	    if ($added) {
		@filtered_hsps = &filter_hsps($seqid,
					      "other",
					      @hsps);
		$coverage = &compute_coverage(@filtered_hsps);
	    }
	}
    }
    if ($coverage =~ m/X\s{${min_gap_seq},}X/ or # Look for an internal gap in coverage
	$coverage =~ m/^\s{${min_gap_seq},}X/ or # Look for 5'  gap in coverage)
	$coverage =~ m/X\s{${min_gap_seq},}$/ ) { # Look for 3'  gap in coverage)
	&emit_gap_seq($seqid, $coverage, $seq);
    }
    if ($coverage =~ m/X\s{${gap_limit},}X/ or # Look for internal gap in coverage
	$coverage =~ m/^\s{${gap_limit},}X/ or # Look for missing coverages at ends.
	$coverage =~ m/X\s{${gap_limit},}$/ ) { 
	print "Seqid $seqid has an internal gap or an uncovered end that is too large.\n";
	$big_gap_count += $seqid_count{$seqid};
	$tile_failure_seqs{$seqid} = 1;
	$status = "an internal gap or an uncovered end that is too large.";
	&print_hsps(\*STDOUT,
		    $seqid,
		    $status,
		    @filtered_hsps);
	if ($debug) {
	    print "All HSPS:\n";
	    &print_hsps(\*STDOUT, $seqid, $status, @hsps);
	}
	if (not $add_unknowns) {
	    if ($details) {
		&print_hsps(\*HSPS, $seqid, $status, @hsps);
		print TILEDETAILS join("\t",
				       $seqid,
				       $seqid_count{$seqid},
				       $status,
				       ""), "\n";
	    }
	    return;
	}
	else {
	    my $offset = 0;
	    while ((my $start = index($coverage, ' ', $offset)) > -1) {
		my $stop = index($coverage, 'X', $start + 1);
		if ($stop == -1) {
		    $stop = length($coverage) - 1;
		}
		$offset = $stop + 1;
		my $len = $stop - $start + 1;
		if ($len >= $min_gap_seq) {
		    my $hsp = { "seq_length" => length($seq),
				    "db" => "unknown",
				    "key" => sprintf("unknown[%d]", $len),
				    "length" => $len,
				    "frac_id" => 1.0,
				    "q_start" => $start + 1,
				    "q_stop" => $stop,
				    "s_start" => 1,
				    "s_stop" => $len,
				    "ss" => '+',
				    "n_id" => $len,
				    "gaps" => 0,
				    "used" => 0};
		    push (@hsps, $hsp);
		}
	    }
	    @filtered_hsps = &filter_hsps($seqid,
					  "unknowns",
					  @hsps);
	    $tiles_with_unknown_count += $seqid_count{$seqid};
	}
    }
    my $root = $filtered_hsps[0];
    for (my $i = 0; $i < scalar(@filtered_hsps); $i++) {
	my $q_start = $filtered_hsps[$i]->{"q_start"};
	my $q_stop = $filtered_hsps[$i]->{"q_stop"};
	my $q_len = $q_stop - $q_start + 1;
	$filtered_hsps[$i]->{"next"} = [];
	for (my $j = $i+1; $j < scalar(@filtered_hsps); $j++) {
	    my $gap = abs($q_stop + 1 - $filtered_hsps[$j]->{"q_start"});
	    if ($gap < $gap_limit) {
		push (@{$filtered_hsps[$i]->{"next"}}, $j);
	    }
	}
    }
    my $iters = 0;
    my $obj = &find_best_path(\@filtered_hsps, \$iters, 0 );
    if ($obj->{"error"}) {
	if ($obj->{"error"} == 1) {
	    print "Gap too large for $seqid\n";
	    $status = "Gap too large";
	    $error_gap_count += $seqid_count{$seqid};
	}
	elsif ($obj->{"error"} == 2) {
	    print "Recursion error for $seqid\n";
	    $status = "Recursion error";
	    $recursion_error_count += $seqid_count{$seqid};
	}
	else {
	    printf "Unknown error (%d) for $seqid\n", $obj->{"error"};
	    $unknown_tiling_error_count += $seqid_count{$seqid};
	    $status = "Unknown error";
	}
	$tile_failure_seqs{$seqid} = 1;
	&print_hsps(\*STDOUT,
		    $seqid,
		    $status,
		    @filtered_hsps);
	if ($debug) {
	    print "All HSPS:\n";
	    &print_hsps(\*STDOUT,
			$seqid,
			$status,
			@hsps);
	}
	if ($details) {
	    &print_hsps(\*HSPS,
			$seqid,
			$status,
			@filtered_hsps);
	    print TILEDETAILS join("\t",
				   $seqid,
				   $seqid_count{$seqid},
				   $status,
				   ""), "\n";
	}
    }
    else {
	printf "Score for $seqid is %d\n", $obj->{"score"};
	my @components;
	my @path = @{$obj->{"path"}};
	my @read_pieces;
	foreach my $ind (@path) {
	    my $hsp = $filtered_hsps[$ind];
	    printf "  [%d-%d (%s)] %5.2f%% %d gaps %s\n",
	    $hsp->{"q_start"},
	    $hsp->{"q_stop"},
	    $hsp->{"ss"},
	    $hsp->{"frac_id"}*100,
	    $hsp->{"gaps"},
	    $hsp->{"key"};
	    $feat_counts{$hsp->{"key"}}++;
	    $hsp->{"used"} = 1;
	    push(@components, sprintf("%s(%s)", $hsp->{"key"}, $hsp->{"ss"}));
	    push(@read_pieces, [ $hsp->{"key"}, $hsp->{"q_start"}, $hsp->{"q_stop"} ]);
	}
	if ($add_gaps) {
	    # This block of code adds the portions of the reads
	    # outside of the matched components. At the ends, we add
	    # the unmatched parts.
	    my $last_piece_ind = $#read_pieces;
	    $read_pieces[0]->[1] = 1;
	    $read_pieces[$last_piece_ind]->[2] = length($seq);
	    for (my $i = 0; $i < $last_piece_ind; $i++) {
		my $i_end_piece = $read_pieces[$i]->[2];
		my $i_start_next_piece = $read_pieces[$i+1]->[1];
		# If the end of one component overlaps with the
		# beginning of the next one, then we'll leave it
		# alone. We only need to handle gaps.  Also, if one of
		# the components around a gap is an ITR and one is
		# not, then assign the gap to the non-ITR piece. If both
		# components are ITR's, then split the gap in the middle.
		if ($i_end_piece + 1 < $i_start_next_piece) {
		    my $gap_size = $i_start_next_piece - $i_end_piece - 1;
		    my $split_point;
		    if ($read_pieces[$i]->[0] =~ m/itr/i and
			$read_pieces[$i+1]->[0] =~ m/itr/i) {
			$split_point = "middle";
			if ($gap_size > 0) {
			    print STDERR "Double ITR found for $seqid.\n";
			}
		    }
		    elsif ($read_pieces[$i]->[0] =~ m/itr/i) {
			$split_point = "left";
		    }
		    elsif ($read_pieces[$i+1]->[0] =~ m/itr/i) {
			$split_point = "right";
		    }
		    else {
			$split_point = "middle";
		    }
		    my $i_new_end;
		    if ($split_point eq "left") {
			$i_new_end = $i_end_piece;
		    }
		    elsif ($split_point eq "middle") {
			$i_new_end = int($gap_size/2) + $i_end_piece;
		    }
		    else {
			$i_new_end = $gap_size + $i_end_piece;
		    }
		    $read_pieces[$i]->[2] = $i_new_end;
		    $read_pieces[$i+1]->[1] = $i_new_end + 1;
		}
	    }
	}
	my $last_ind = $#path;
	if ($consolidate) {
	    for (my $i = 0; $i < scalar(@components); $i++) {
	    # Apply bespoke rules to reduce complexity.
		$components[$i] =~ s/ITR-FLIP/ITR/;
		$components[$i] =~ s/\([tf]\)//g;
		my ($start, $stop) = ($components[$i] =~ m/\[(\d+)-(\d+)\]/);
		if ($components[$i] =~ /ITR/) {
		    if (($start >= 1 and $start <= 6 and
			 $stop >= 140 and $stop <= 165) or
			($stop >= 160 and $start <= 25)) {
			$components[$i] =~ s/\[\d+-\d+\]//;
		    }
		}
	    }
	}
	if (not $quant_homopolymers or $consolidate) {
	    for (my $i = 0; $i < scalar(@components); $i++) {
		if ($components[$i] =~ m/poly(.)\[\d+\](\([tf]\))/) {
		    $components[$i] = "poly" . $1;
		}
	    }
	}
	my $seq_code = join(" ", @components);
	if ($by_component) {
	    &output_components($seqid, \@read_pieces, \@components);
	}
	if (not exists($seqid_count{$seqid})) {
	    die "Unable to find seqid count for $seqid.\n";
	}
	$code_counts{$seq_code} += $seqid_count{$seqid};
	push (@{$code_seqids{$seq_code}}, $seqid);
	print "Sequence code: ", $seq_code, "\n";
	my $start_diff = $filtered_hsps[0]->{"q_start"} - 1;
	my $end_diff = $filtered_hsps[$path[$last_ind]]->{"seq_length"} - $filtered_hsps[$path[$last_ind]]->{"q_stop"};
	if ( $start_diff <= $end_gap and $end_diff <= $end_gap ) {
	    $status = "Tiled";
	    print "  Complete\n\n";
	    $complete_tiling_count += $seqid_count{$seqid};
	}
	else {
	    $status = "Partial tiled with end gaps";
	    printf "  Partial: seq_length = %d\n\n", $filtered_hsps[$last_ind]->{"seq_length"};
	    $partial_tiling_count += $seqid_count{$seqid};
	    &print_hsps(\*STDOUT,
			$seqid,
			$status,
			@filtered_hsps);
	    if ($debug) {
		print "All HSPS:\n";
		&print_hsps(\*STDOUT,
			    $seqid,
			    $status,
			    @hsps);
	    }
	}
	if ($details) {
	    &print_hsps(\*HSPS,
			$seqid,
			$status,
			@filtered_hsps);
	    print TILEDETAILS join("\t",
				   $seqid,
				   $seqid_count{$seqid},
				   $status,
				   $seq_code), "\n";
	}
	if ($by_strand) {
	    foreach my $name (split(/:/, $names{$seqid})) {
		my $zmw = (split(/\//, $name))[1];
		if (not exists($codes_for_zmw{$zmw})) {
		    $codes_for_zmw{$zmw} = $seq_code;
		}
		elsif ($codes_for_zmw{$zmw} eq $seq_code) {
		    $codes_for_zmw{$zmw} .= " x 2";
		}
		elsif ($seq_code lt $codes_for_zmw{$zmw}) {
		    $codes_for_zmw{$zmw} = "$seq_code U $codes_for_zmw{$zmw}";
		}
		else {
		    $codes_for_zmw{$zmw} = "$codes_for_zmw{$zmw} U $seq_code";
		}
	    }
	}
    }
}

if ($by_strand) {
    my %zmw_code_counts;
    my $lib_total = 0;
    foreach my $zmw (keys %codes_for_zmw) {
	$zmw_code_counts{$codes_for_zmw{$zmw}}++;
	$lib_total++;
    }
    open (ZCOUNTS, ">${library}.tile.zmw.counts") ||
	die "Unable to open for writing ${library}.tile.zmw.counts: $!\n";
    foreach my $seq_code (sort {$zmw_code_counts{$b} <=> $zmw_code_counts{$a} ||
				    $a cmp $b}  keys %zmw_code_counts) {
	printf ZCOUNTS "%7d %9.7f %s\n", $zmw_code_counts{$seq_code}, $zmw_code_counts{$seq_code} / $lib_total, $seq_code;
    }
}
    
if (scalar(keys %feat_counts) > 0) {
    open (FEATS, ">${library}.feats") ||
	die "Unable to open for writing ${library}.feats: $!\n";
    printf FEATS "%20s %s\n", "Feature", "Counts";
    foreach my $feat (sort {$feat_counts{$b} <=> $feat_counts{$a}} keys %feat_counts) {
	printf FEATS "%20s %d\n", $feat, $feat_counts{$feat};
    }
}

my $count_total = 0;
foreach my $count (values %code_counts) {
    $count_total += $count;
}
print SUMMARY "Total sequence code count = $count_total\n";

# The next block of code does two things -- counts total number of sequences
# in the library and emits any sequences that had no Blast hits at all or which were not tiled.

my $lib_total = 0;
foreach my $seqid (keys %seqid_count) {
    $lib_total += $seqid_count{$seqid};
    if (not exists($seqids_with_hsps{$seqid}) or exists($tile_failure_seqs{$seqid})) {
	&emit_unmatched_seq($seqid);
    }
}

print SUMMARY "Total sequence library count = $lib_total\n";
my $n_unaccounted = $lib_total - $count_total + $tiles_with_unknown_count;
if ($add_unknowns) {
    printf SUMMARY "Number of sequences where 'unknown' was added to the tiling pattern is %d\n", $tiles_with_unknown_count;
}
printf SUMMARY "Unaccounted sequences number %d which is %7.4f %%\n", $n_unaccounted, $n_unaccounted/$lib_total * 100;
printf SUMMARY "Number of big gap sequences is %d\n", $big_gap_count;
printf SUMMARY "Number of oversize gap errors is %d\n", $error_gap_count;
printf SUMMARY "Number of recursion errors in tiling is %d\n", $recursion_error_count;
printf SUMMARY "Number of unknown errors in tiling is %d\n", $unknown_tiling_error_count;
printf SUMMARY "Number of complete tilings found is %d\n", $complete_tiling_count;
printf SUMMARY "Number of partial tilings found is %d\n", $partial_tiling_count;
if ($by_component) {
    printf SUMMARY "Number of split references reads is %d\n", $split_references_count;
}
printf SUMMARY "Number of completely unknown sequences is %d\n", $total_unknown_count;

my $iseq = 0;
if (scalar(keys %code_counts) > 0) {
    $write_seq_limit = min($write_seq_limit, scalar(keys %code_counts));
    open (COUNTS, ">${library}.tile.counts") ||
	die "Unable to open for writing ${library}.tile.counts: $!\n";
    foreach my $seq_code (sort {$code_counts{$b} <=> $code_counts{$a} ||
				    $a cmp $b}  keys %code_counts) {
	my $freq = $code_counts{$seq_code} / $lib_total;
	printf COUNTS "%7d %9.7f %s\n", $code_counts{$seq_code}, $freq, $seq_code;
	if (++$iseq <= $write_seq_limit) {
	    my $dir = "${library}.seqs/" . &make_dir_for_iseq($iseq);
	    if (not -d $dir) {
		system_with_check("mkdir -p $dir", 1)
	    }
	    my $ext = $fastq eq "" ? "fasta" : "fastq";
	    open (SEQOUT, ">$dir/${library}_${iseq}.${ext}") ||
		die "Unable to open $dir/${library}_${iseq}.${ext} for writing";
	    my $i = 0;
	    foreach my $seqid (@{$code_seqids{$seq_code}}) {
		my ($seq) = retrieve_sequence_data($pg, $seqid);
		if ($raw or $fastq ne "") {
		    foreach my $name ($pg->get_array_for_field("select name " .
							       "  from raw_seqs " .
							       " where seqid = " . quotify($seqid) .
							       "   and library = " . quotify($library))) {
			if ($fastq eq "") {
			    print SEQOUT ">$name\n";
			    print SEQOUT format_seq($seq);
			}
			else {
			    if (not exists($fastq_data{$name})) {
				die "Unable to find $name in $fastq\n";
			    }
			    my ($annotation, $seq, $qual) = @{$fastq_data{$name}};
			    print SEQOUT format_fastq($name, $annotation, $seq, $qual);
			}
		    }
		}
		else {
		    if (++$i == 1) {
			printf SEQOUT ">$seqid number of reads = %d Code count = %d Frequency = %9.7f %s\n",
			    $seqid_count{$seqid},
			    $code_counts{$seq_code},
			    $freq,
			    $seq_code;
		    }
		    else {
			printf SEQOUT ">$seqid number of reads = %d\n",
			    $seqid_count{$seqid};
		    }
		    print SEQOUT format_seq($seq);
		}
	    }
	    close (SEQOUT);
	}
    }
}

if (scalar(keys %unidentified_seqs) > 0) {
    open (GAPCOUNTS, ">${library}.gap.counts") ||
	die "Unable to open for writing ${library}.gap.counts: $!\n";
    foreach my $seq (sort {$unidentified_seqs{$b} <=> $unidentified_seqs{$a} ||
			       $a cmp $b}  keys %unidentified_seqs) {
	printf GAPCOUNTS "%7d %s\n", $unidentified_seqs{$seq}, $seq;
    }
}

if ($details) {
    my $rows = $pg->get_all_rows("select seqid, name " .
				 "  from raw_seqs " .
				 " where library = " . quotify($library) .
				 " order by seqid, name");
    foreach my $row (@{$rows}) {
        print NAME2SEQID join("\t", @{$row}), "\n";
    }
}

sub add_homopolymers {
    my $seq = shift;
    my $seqid = shift;

    my $done = 0;
    my $found = 0;
    until ($done) {
	if ($seq =~ m/([acgt])\1{$min_poly,}/ip) {
	    my $base = uc($1);
	    my $start = length(${^PREMATCH}) + 1;
	    my $len = length(${^MATCH});
	    my $stop = $start + $len - 1;
	    substr($seq, $start - 1, $len) = "X" x $len;
	    # Check to see if there are additional runs after this one.
	    my $done_extension = 0;
	    until ($done_extension) {
		if ($seq =~ m/([acgt])\1{$min_poly,}/ip) {
		    my $new_start = length(${^PREMATCH}) + 1;
		    my $new_base = uc($1);
		    if ($stop + 2 == $new_start) {
			$len = $len + 1 + length(${^MATCH});
			$stop = $start + $len - 1;
			substr($seq, $start - 1, $len) = "X" x $len;
		    }
		    else {
			$done_extension = 1;
		    }
		}
		else {
		    $done_extension = 1;
		}
	    }
#	    print STDERR "In $seqid: start = $start  len = $len  stop = $stop\n";
	    $found = 1;
	    my $hsp = { "seq_length" => length($seq),
			"db" => "poly",
			"key" => sprintf("poly%s[%d]", $base, $len),
			"length" => $len,
			"frac_id" => 1.0,
			"q_start" => $start,
			"q_stop" => $stop,
			"s_start" => 1,
			"s_stop" => $len,
			"ss" => '+',
			"n_id" => $len,
			"gaps" => 0,
			"used" => 0};
	    push (@hsps, $hsp);
	}
	else {
	    $done = 1;
	}
    }
}

sub filter_hsps {
    my $seqid = shift;
    my $operation = shift;
    my @sorted_hsps = sort {$a->{"q_start"} <=> $b->{"q_start"} ||
				$a->{"q_stop"} <=> $b->{"q_stop"} } @_;
    my @delete_me = (0) x scalar(@sorted_hsps);
    # I'm not certain if we can optimize this code by exiting loops when an array
    # element is marked for deletion, so I'm not going to try. The efficiency loss
    # is negligible.
    for (my $i = 0; $i < scalar(@sorted_hsps) - 1; $i++) {
	next if $delete_me[$i];
	for (my $j = $i + 1; $j < scalar(@sorted_hsps); $j++) {
	    next if $delete_me[$j];
	    # Check for exact overlaps and then keep the better one.
	    if ($sorted_hsps[$i]->{"q_start"} == $sorted_hsps[$j]->{"q_start"} and
		$sorted_hsps[$i]->{"q_stop"} == $sorted_hsps[$j]->{"q_stop"}) {
		# Best score is 0 and increasing score means worse HSP.
		my $i_score = $sorted_hsps[$i]->{"length"} - $sorted_hsps[$i]->{"n_id"} + $sorted_hsps[$i]->{"gaps"};
		my $j_score = $sorted_hsps[$j]->{"length"} - $sorted_hsps[$j]->{"n_id"} + $sorted_hsps[$j]->{"gaps"};
		my $i_key = $sorted_hsps[$i]->{"db"} . $sorted_hsps[$i]->{"key"};
		my $j_key = $sorted_hsps[$j]->{"db"} . $sorted_hsps[$j]->{"key"};
		if ($i_score > $j_score) {
		    $delete_me[$i] = 1;
		}
		elsif ($i_score < $j_score) {
		    $delete_me[$j] = 1;
		}
		elsif ($i_key lt $j_key) {
		    $delete_me[$i] = 1;
		}
		elsif ($j_key lt $i_key) {
		    $delete_me[$j] = 1;
		}
		elsif ($sorted_hsps[$i]->{"ss"} ne $sorted_hsps[$j]->{"ss"}) {
		    # Save the positive strand.
		    if ($sorted_hsps[$i]->{"ss"} eq '-') {
			$delete_me[$i] = 1;
		    }
		    else {
			$delete_me[$j] = 1;
		    }
		}
		else {
		    $delete_me[$i] = 1;
		    print STDERR "Unable to pick from two exactly overlapping HSPS:\n";
		    print STDERR "i = $i  j = $j\n";
		    &print_hsps(\*STDERR,
				$seqid,
				"filtering",
				@sorted_hsps);
		}
		next;
	    }
	    # Check for complete overlaps and remove the smaller entity.
	    if (range_overlap($sorted_hsps[$i]->{"q_start"},
			      $sorted_hsps[$i]->{"q_stop"},
			      $sorted_hsps[$j]->{"q_start"},
			      $sorted_hsps[$j]->{"q_stop"})) {
		if ($sorted_hsps[$i]->{"q_start"} >= $sorted_hsps[$j]->{"q_start"} and
		    $sorted_hsps[$i]->{"q_stop"} <= $sorted_hsps[$j]->{"q_stop"}) {
		    $delete_me[$i] = 1;
		    # printf "Deleting $i; $i = [%d, %d]  $j = [%d, %d]\n",
		    # 	$sorted_hsps[$i]->{"q_start"},
		    # 	$sorted_hsps[$i]->{"q_stop"},
		    # 	$sorted_hsps[$j]->{"q_start"},
		    # 	$sorted_hsps[$j]->{"q_stop"}
		}
		elsif ($sorted_hsps[$j]->{"q_start"} >= $sorted_hsps[$i]->{"q_start"} and
		       $sorted_hsps[$j]->{"q_stop"} <= $sorted_hsps[$i]->{"q_stop"}) {
		    $delete_me[$j] = 1;
		    # printf "Deleting $j; $i = [%d, %d]  $j = [%d, %d]\n",
		    # 	$sorted_hsps[$i]->{"q_start"},
		    # 	$sorted_hsps[$i]->{"q_stop"},
		    # 	$sorted_hsps[$j]->{"q_start"},
		    # 	$sorted_hsps[$j]->{"q_stop"}
		}
	    }
	}
    }
    my @filtered_hsps;
    for (my $i = 0; $i < scalar(@sorted_hsps); $i++) {
	if (not $delete_me[$i]) {
	    push (@filtered_hsps,
		  $sorted_hsps[$i]);
	}
    }
    return @filtered_hsps;
}

sub compute_coverage {
    my @filtered_hsps = @_;
    my $root = $filtered_hsps[0];
    my $coverage = " " x $root->{"seq_length"};
    for (my $i = 0; $i < scalar(@filtered_hsps); $i++) {
	my $q_start = $filtered_hsps[$i]->{"q_start"};
	my $q_stop = $filtered_hsps[$i]->{"q_stop"};
	my $q_len = $q_stop - $q_start + 1;
	substr($coverage,
	       $q_start - 1,
	       $q_len) = "X" x $q_len;
    }
    return $coverage;
}

sub find_best_path {
    my $hspsp = shift;
    my $itersp = shift;
    my $ind = shift;
    my $ret = { "path" => [ ],
		"score" => 2000000000,
		"error" => 0 };
		
    # Error codes: 1 = gap too big. 2 = recursion error.

    my $last_ind = scalar(@{$hspsp}) - 1;
    if (++$$itersp > $max_iter) {
	print "Max_iter = $max_iter exceeded.\n";
	$ret->{"error"} = 2;
    }
    else {
	my $hsp = $hspsp->[$ind];
	if (scalar(@{$hsp->{"next"}}) == 0) {
	    $ret->{"path"} = [ $ind ];
	    # Terminal score is the number of mismatches + number of gaps + gap at end.
	    $ret->{"score"} = abs($hsp->{"length"} - $hsp->{"n_id"}) + $hsp->{"gaps"} + $hsp->{"seq_length"} - $hsp->{"q_stop"}
	}
	else {
	    my $best_score = 2000000000;
	    my $best_ind;
	    my $best_obj;
	    foreach my $nexti (@{$hspsp->[$ind]->{"next"}}) {
		my $add_score = abs($hsp->{"q_stop"} + 1 - $hspsp->[$nexti]->{"q_start"});
		next if $add_score >= $gap_limit;
		$add_score += abs($hspsp->[$nexti]->{"length"} -
				  $hspsp->[$nexti]->{"n_id"}) +
				  $hspsp->[$nexti]->{"gaps"};
		my $obj = &find_best_path($hspsp,
					  $itersp,
					  $nexti);
		if ($obj->{"error"}) {
		    $ret->{"error"} = $obj->{"error"};
		    return $ret;
		}
		my $new_score = $obj->{"score"} + $add_score;
		if ($new_score < $best_score) {
		    $best_score = $new_score;
		    $best_ind = $nexti;
		    $best_obj = $obj;
		}
	    }
	    if ($best_score == 2000000000) {
		$ret->{"error"} = 1;
	    }
	    else {
		$ret = $best_obj;
		unshift(@{$ret->{"path"}}, $ind);
		$ret->{"score"} = $best_score;
	    }
	}
    }
    return $ret;
}


sub emit_gap_seq {
    my $seqid = shift;
    my $coverage = shift;
    my $seq = shift;
    
    my $pos = 0;
    my $igap = 0;
    while ((my $ind = index($coverage, ' ', $pos)) >= 0) {
	$igap++;
	my $start = $ind + 1;
	my $last = index($coverage, 'X', $start);
	if ($last < 0) {
	    $last = length($coverage);
	}
	my $len = $last - $ind;
	if ($len >= $min_gap_seq) {
	    my $gap_seq = substr($seq, $ind, $last - $ind);
	    $unidentified_seqs{$gap_seq} += $seqid_count{$seqid};
	    printf GAPSEQ ">${seqid}_%d_%d\n", $start, $last + 1;
# 	    print GAPSEQ format_seq($gap_seq);
	    print GAPSEQ $gap_seq, "\n";
	}
	$pos = $last;
	last if $last == length($coverage);
    }
}

sub emit_unmatched_seq {
    my $seqid = shift;
    
    my ($seq) = retrieve_sequence_data($pg, $seqid);
    print UNMATCHED ">${seqid}\n";
    print UNMATCHED format_seq($seq);
}

sub print_hsps {
    my $fh = shift;
    my $seqid = shift;
    my $status = shift;
    my @sorted_hsps = sort {$a->{"q_start"} <=> $b->{"q_start"} ||
				$a->{"q_stop"} <=> $b->{"q_stop"} } @_;
    foreach my $hsp (@sorted_hsps) {
	print $fh join("\t",
		       $seqid,
		       $status,    
		       $hsp->{"seq_length"},
		       $hsp->{"db"},
		       $hsp->{"key"},
		       $hsp->{"length"},
		       $hsp->{"frac_id"},
		       $hsp->{"q_start"},
		       $hsp->{"q_stop"},
		       $hsp->{"ss"},
		       $hsp->{"s_start"},
		       $hsp->{"s_stop"},
		       $hsp->{"used"},
		       $hsp->{"n_id"},
		       $hsp->{"gaps"}), "\n";
    }
}

sub make_dir_for_iseq {
    # Construct a multilevel directory name for the passed sequence number
    # assuming there are $write_seq_limit total.
    my $iseq = shift;
    $iseq--; 			# We need to have sequence numbers start at 0;
    die "Programming error: write_seq_limit is not positive.\n"
	if $write_seq_limit < 1;
    my @levels;
    my $nlevels = POSIX::ceil(log($write_seq_limit) / log(100));
    for (my $i = 1; $i <= $nlevels; $i++) {
	my $piece = sprintf("%02d", $iseq % 100);
	unshift(@levels, $piece);
	$iseq = int($iseq / 100);
    }
    if ($nlevels > 0) {
	pop(@levels);
    }
    return join("/", @levels);
}

sub open_file {
    my $file = shift;
    my $fh;

    open ($fh, "$file") ||
	die "Unable to open $file: $!\n";
    return $fh;
}

sub output_components {
    my $seqid = shift;
    my $read_pieces_p = shift;
    my $components_p = shift;
    my @read_pieces = @{$read_pieces_p};
    my @components = @{$components_p};

    if (scalar(@read_pieces) != scalar(@components)) {
	die "Programming error: mismatch of read_pieces and components sizes\n";
    }
    my $prev_reference = "";
    foreach my $piece (@read_pieces) {
	my $reference = $piece->[0];
	$reference =~ s/\[\d+-\d+\]$//;
	if ($prev_reference ne "") {
	    if ($reference eq $prev_reference) {
		if ($debug) {
		    print "$seqid had split references.\n";
		}
		$split_references_count += $seqid_count{$seqid};
		printf SPLIT_READS ">$seqid count = %d\n", $seqid_count{$seqid};
		my ($seq) = retrieve_sequence_data($pg, $seqid);
		print SPLIT_READS format_seq($seq);
		return;
	    }
	}
	$prev_reference = $reference;
    }
    foreach my $piece (@read_pieces) {
	my ($reference, $start, $stop) = @{$piece};
	next if $reference =~ m/^poly/;
	$reference =~ s/^([^]]+)\[.*/$1/;
	if (not exists($component_files{$reference})) {
	    my $ref_file = "${library}.${reference}.fastq";
	    $component_files{$reference} = open_file(">$ref_file");
	}
	my $fh = $component_files{$reference};
	foreach my $name ($pg->get_array_for_field("select name " .
						   "  from raw_seqs " .
						   " where seqid = " . quotify($seqid) .
						   "   and library = " . quotify($library))) {
	    if (not exists($fastq_data{$name})) {
		die "Unable to find $name in $fastq\n";
	    }
	    my ($annotation, $seq, $qual) = @{$fastq_data{$name}};
	    $seq = substr($seq, $start - 1, $stop - $start + 1);
	    $qual = substr($qual, $start - 1, $stop - $start + 1);
	    $name .= sprintf("_%d_%d", $start, $stop);
	    print $fh format_fastq($name, $annotation, $seq, $qual);
	}
    }
}
