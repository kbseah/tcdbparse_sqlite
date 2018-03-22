#!/usr/bin/env perl

=head1 NAME

tcdbparse_sqlite.pl - Utility to manage data for transporter analysis with the
TCDB annotations in an SQLite database

=head1 SYNOPSIS

B<perl tcdbparse_sqlite.pl> -db dbfile -task setup -tcfam tc_family_names.tsv -uptake uptake_transporters.tsv

B<perl tcdbparse_sqlite.pl> -db dbfile -task read -blast blast_files.tsv -tmhmm tmhmm_files.tsv

B<perl tcdbparse_sqlite.pl> -db dbfile -task report -out output_suffix

B<perl tcdbparse_sqlite.pl> -help

B<perl tcdbparse_sqlite.pl> -man

=head1 DESCRIPTION

tcdbparse_sqlite is an utility to set up and populate an SQLite database for
transporter analysis of microbial genomes from Blast results vs. the Transporter
Classification Database (TCDB) and from TMHMM results.

The intended use is for reporting counts of proteins per genome that have hits
to a shortlist of uptake transporters, however, the SQLite database can also be
used for other types of analysis, and makes it easier to manage the data.

The required Blast and TMHMM output files can be generated with the script
B<do_blastp_tmhmm.pl>.

=cut

use strict;
use warnings;
use 5.010;

use DBI;
use Getopt::Long;
use Pod::Usage;
#use Data::Dumper; # For diagnosing data problems

my $task;
my $dbfile;
# Params for task setup
my ($tcfam,$tcsubfam,$uptakefam); 
# Params for task read
my ($tmhmm_tsv, $blastout_tsv);
# Params for task report
my ($blast_minqcovs,$blast_minpident) = (70, 30);
my $out = 'tcdbparse_uptake_transporter';

## INPUT OPTIONS

pod2usage(verbose=>0,exitval=>2) if !@ARGV;

GetOptions ("task=s"=> \$task,
            "tcfam=s" => \$tcfam,
            "tcsubfam=s" => \$tcsubfam,
            "uptake=s" => \$uptakefam,
            "db=s" =>\$dbfile,
            "tmhmm=s" => \$tmhmm_tsv,
            "blast=s" => \$blastout_tsv,
            "qcov=i" => \$blast_minqcovs,
            "pid=i" => \$blast_minpident,
            "out=s" => \$out,
            'help' => sub { pod2usage(verbose=>1,exitval=>1) },
            'man' => sub { pod2usage(verbose=>2,exitval=>0) },
            ) or pod2usage(2);

=head1 ARGUMENTS

=head2 Common arguments

=over 8

=item -task I<setup|read|report>

There are three functions, specified with this option:
B<setup> (create database with tables, and read in list of TC family names and
the shortlist of uptake transporter families), B<read> (read in Blast and TMHMM
result files), and B<report> (output tables of summary statistics for further
analysis).

=item -db F<DBFILE>

Name of SQLite database file

=item -help

Short help message

=item -man

This manual page

=back

=head2 Arguemnts for task "setup"

=over 8

=item -tcfam F<tc_family_names.tsv>

Tab-separated table of TC family identifiers (column 1) and names (column 2).
Lines beginning with # are treated as comments and ignored. Other columns ignored.

=item -uptake F<uptake_transporters.tsv>

Tab-separated table of uptake-related transporter TC family identifiers (column
1) and names (column 2). Lines beginning with # are treated as comments and
ignored. Other columns ignored.

=back

=head2 Arguments for task "read"

=over 8

=item -blast F<blast_output_files.tsv>

Tab-separated table of Blastp output files (tabular format 6 'std qcovs') of
amino acid sequences (one file per genome) vs. TCDB sequences. Column 1 are the
genome names/identifiers, and column 2 are the blastp output file names.

=item -tmhmm F<tmhmm_output_files.tsv>

Tab-separated table of TMHMM output files of amino acid sequences (one file per
genome) vs. TCDB sequences. Column 1 are the genome names/identifiers, and
column 2 are the TMHMM output filenames.

=back

=head2 Arguments for task "report"

=over 8

=item -out I<STRING>

Prefix for output files.

=item -qcov I<NUMERIC>

Minimum percentage coverage of Blastp hit vs. TCDB to be included in results.

=item -pid I<NUMERIC>

Minimum percentage identity of Blastp hit vs. TCDB to be included in results.

=back

=cut 

## MAIN ########################################################################

die ("Please specify database file to parameter -db") if !$dbfile;

## SETUP DATABASE
if ($task eq 'setup') {
    # Check for input files
    die ('Please supply table of TC family names') if !$tcfam;
    die ('Please supply shortlist of uptake transporters') if !$uptakefam;
    
    # Hash TC family and subfamily names from TSV file shortlists
    my $tchref = tbl2hash ($tcfam);
    my $uptakehref = tbl2hash ($uptakefam);
    
    # Set up database
    my $setup_return = setup_DB ($dbfile,
                                 $tchref,
                                 $uptakehref);
    
}
## READ BLAST AND TMHMM RESULTS INTO DB
elsif ($task eq 'read') {
    # Check for input files
    die ('Please supply table of Blastp output files') if !$blastout_tsv;
    die ('Please suppyl table of tmhmm output files') if !$tmhmm_tsv;
    
    # Connect to SQLite database
    say STDERR "Connecting to SQLite database file";
    my $dsn = "dbi:SQLite:dbname=$dbfile";
    my $user = "";
    my $password = "";
    my $dbh = DBI->connect ($dsn, $user, $password, {  # Database handle
        PrintError => 0,
        RaiseError => 1,
        AutoCommit => 0,
        FetchHashKeyName => 'NAME_lc'
    });
    
    # Hash lists of Blast and TMHMM results
    my $blastfileshref = tbl2hash ($blastout_tsv);
    my $tmhmmfileshref = tbl2hash ($tmhmm_tsv);
    
    # For each blast file and TMHMM file, read results
    foreach my $genome (keys %$blastfileshref) {
        say STDERR "Reading Blast and TMHMM results for genome $genome from file $blastfileshref->{$genome}";
        my $blastparse_href = read_blastout6 ($blastfileshref->{$genome});
        read_tmhmm_add_to_hash($tmhmmfileshref->{$genome},$blastparse_href);
        #print Dumper($blastparse_href);
        insert_blast_tmhmm_hash($dbh,$blastparse_href,$genome);
        $dbh->commit or die $dbh->errstr;
    }
    
    # Disconnect DB file
    say STDERR "Disconnecting from database";
    $dbh->disconnect;
}
## REPORT COUNTS FROM DB
elsif ($task eq 'report') {
    # Check for input parameters

    # Connect to SQLite database
    say STDERR "Connecting to SQLite database file";
    my $dsn = "dbi:SQLite:dbname=$dbfile";
    my $user = "";
    my $password = "";
    my $dbh = DBI->connect ($dsn, $user, $password, {  # Database handle
        PrintError => 0,
        RaiseError => 1,
        AutoCommit => 0,
        FetchHashKeyName => 'NAME_lc'
    });
    
    # Write report for all hits
    open(my $report_all_fh, ">", "$out\.all.tsv") or die ("$!");
    report_uptake_all($dbh,$report_all_fh);
    close ($report_all_fh);
    
    # Write report summarizing all hits above cutoffs, total and in shortlist
    open(my $report_summary_fh, ">", "$out\.summary.tsv") or die ("$!");
    report_uptake_summary($dbh,$report_summary_fh);
    close ($report_summary_fh);

    # Write report summarizing uptake-TC families per genome
    open(my $report_summary2_fh, ">", "$out\.uptakepergenome.tsv") or die ("$!");
    report_uptake_pergenome($dbh,$report_summary2_fh);
    close ($report_summary2_fh);
    
    # Disconnect DB file
    say STDERR "Disconnecting from database";
    $dbh->disconnect;
    
} else {
    say STDERR '\"task\" should be one of: setup, read, report';
}


## SUBROUTINES #################################################################

sub tbl2hash {
    # Convert TSV table to hash
    # First column is key, second col is value
    # Other columns ignored
    # Lines beginning with # char are comments and ignored
    # Returns a hash ref
    my ($file) = @_;
    my $href;
    open(my $fhin, "<", $file) or die ("Cannot open file $file: $!");
    while (my $line = <$fhin>) {
        chomp $line;
        next if $line =~ m/^#/; # Skip header or comment lines
        my @splitline = split /\t/, $line;
        $href->{$splitline[0]} = $splitline[1];
    }
    close ($fhin);
    return ($href);
}

sub setup_DB {
    # Set up sqlite database and populate TC family name tables
    my ($dbfile,
        $tchref,
        $uptakehref,
        ) = @_;
    # Initialize SQLite database
    say STDERR "Initializing SQLite database file";
    my $dsn = "dbi:SQLite:dbname=$dbfile";
    my $user = "";
    my $password = "";
    my $dbh = DBI->connect ($dsn, $user, $password, {  # Database handle
        PrintError => 0,
        RaiseError => 1,
        AutoCommit => 0, # Set to zero and explicitly commit, otherwise slow
        FetchHashKeyName => 'NAME_lc'
    });
    say STDERR 'Creating table main';
    my $create_main_sql = << 'END_SQL';
CREATE TABLE main (
    id              TEXT PRIMARY KEY,
    genome          TEXT,
    tcdbacc         TEXT,
    tcfull          TEXT,
    tcfam           TEXT,
    tcsubfam        TEXT,
    tcsubsubfam     TEXT,
    pidmult         NUMERIC,
    alnlen          INTEGER,
    pident          NUMERIC,
    qcovs           INTEGER,
    evalue          TEXT,
    bitscore        NUMERIC,
    tmdomains       INTEGER
)
END_SQL

    $dbh->do($create_main_sql);
    $dbh->commit or die $dbh->errstr;
    
    say STDERR 'Creating table tcfam_names';
    my $create_tcfam_sql = << 'END_SQL';
CREATE TABLE tcfam_names (
    tcnum           TEXT PRIMARY KEY,
    name            TEXT
)
END_SQL

    $dbh->do($create_tcfam_sql);
    say STDERR 'Inserting values into table tcfam_names';
    my $counter = 0;
    foreach my $tc (keys %$tchref) {
        $counter++;
        my $sql_cmd = "INSERT INTO tcfam_names (tcnum, name) VALUES (?, ?)";
        $dbh -> do($sql_cmd,
                   undef,
                   $tc, $tchref->{$tc});
        if ($counter % 50 == 0) {
            say STDERR "... inserted $counter entries ..."
        }
    }
    $dbh->commit or die $dbh->errstr;

    
    say STDERR 'Creating table uptakefam_names';
    my $create_uptakefam_sql = << 'END_SQL';
CREATE TABLE uptakefam_names (
    tcnum           TEXT PRIMARY KEY,
    name            TEXT
)
END_SQL

    $dbh->do($create_uptakefam_sql);
    say STDERR 'Inserting values into table uptakefam_names';
    $counter = 0;
    foreach my $tc (keys %$uptakehref) {
        $counter ++;
        my $sql_cmd = "INSERT INTO uptakefam_names (tcnum, name) VALUES (?, ?)";
        $dbh -> do($sql_cmd,
                   undef,
                   $tc, $uptakehref->{$tc});
        if ($counter % 50 == 0) {
            say STDERR "... inserted $counter entries ... ";
        }
    }
    $dbh->commit or die $dbh->errstr;

    # Disconnect DB file
    say STDERR "Disconnecting from database";
    $dbh->disconnect;
}

# Read Blastp results (output format 6 std qcovs)
sub read_blastout6 {
    my ($file) = @_;
    my %hash; # Store intermediate results
    open(IN, "<", $file) or die ("Cannot read file $file: $!");
    my $counter = 0;
    while (<IN>) {
        chomp;
        $counter ++;
        my @theline = split "\t", $_;
        my @tcheader = split /\|/, "$theline[1]";
        my $tc_acc = $tcheader[2]; # TC accession of subject hit
        my @tcfull_split = split /\./, $tcheader[3];
        my $tcfam = join ".", @tcfull_split[0..2];
        my $tcsubfam = join ".", @tcfull_split[0..3];
        my $tcsubsubfam = join ".", @tcfull_split[0..4];
        # TC hit
        $hash{$theline[0]}{"tcdbacc"} = $tc_acc;
        # TC number 
        $hash{$theline[0]}{"tcfull"} = $tcheader[3];
        # TC family
        $hash{$theline[0]}{"tcfam"} = $tcfam;
        # TC subfamily
        $hash{$theline[0]}{"tcsubfam"} = $tcsubfam;
        # TC subsubfamily
        $hash{$theline[0]}{"tcsubsubfam"} = $tcsubsubfam;
        # Percent identity (pident) multiplied by alignment length
        # In case there are multiple alignments on same subject, want to take the mean
        $hash{$theline[0]}{"pidmult"} += $theline[2] * $theline[3];
        # Cumulative alignment length - in case more than one alignment
        $hash{$theline[0]}{"alnlen"} += $theline[3];
        # Coverage of alignment on subject
        $hash{$theline[0]}{"qcovs"} += $theline[12];
        # evalue
        $hash{$theline[0]}{'evalue'} = $theline[10];
        # Strip leading whitespace from bitscore
        $theline[11] =~ s/^\s+//;
        $hash{$theline[0]}{'bitscore'} = $theline[11];
        if ($counter % 100 == 0) {
            say STDERR "... processed $counter Blast results ...";
        }
        
    }
    close(IN);
    foreach my $query (keys %hash) {
        # Get percentage identity of alignments and coverage of alignments on reference
        $hash{$query}{'pident'} = $hash{$query}{"pidmult"} / $hash{$query}{"alnlen"};
    }
    return (\%hash);
}

# Insert Blast results from hash into database table main
sub insert_blast_tmhmm_hash {
    my ($handle,
        $href,
        $genome) = @_;
    my @params = qw(tcdbacc tcfull tcfam tcsubfam tcsubsubfam pidmult alnlen pident qcovs evalue bitscore tmdomains);
    my $params_str = join ", ", ('id', 'genome', @params);

    foreach my $id (keys %$href) {
        my @vals = ($id, $genome);
        foreach my $param (@params) {
            push @vals, $href->{$id}{$param};
        }
        my $sql = "INSERT OR REPLACE INTO main ($params_str) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)";
        $handle->do($sql,
                    undef,
                    @vals);
    }
}

# Read TMHMM results - count number of TM helices per protein
sub read_tmhmm_add_to_hash {
    my ($file, $href) = @_;
    open(IN, "<", $file) or die ("Cannot read file $file: $!");
    while (<IN>) {
        chomp;
        if ($_ =~ /^# (.+) Number of predicted TMHs:\s+(\d+)/) {
            $href->{$1}{'tmdomains'} = $2 if (defined $href->{$1}); # Only add if has a TCDB hit
        }
    }
    close(IN);
}

# Output uptake transporter report - all entries with hits
sub report_uptake_all {
    my ($handle,
        $fh) = @_;
    # print header to file:
    say $fh join "\t", qw(id genome tcdbacc pident qcovs tmdomains tcnum name);
    my $sql = << 'END_SQL';
SELECT id, genome, tcdbacc, pident, qcovs, tmdomains, tcnum, name FROM main JOIN uptakefam_names ON main.tcfam=uptakefam_names.tcnum
UNION ALL 
SELECT id, genome, tcdbacc, pident, qcovs, tmdomains, tcnum, name FROM main JOIN uptakefam_names ON main.tcsubfam=uptakefam_names.tcnum
UNION ALL 
SELECT id, genome, tcdbacc, pident, qcovs, tmdomains, tcnum, name FROM main JOIN uptakefam_names ON main.tcsubsubfam=uptakefam_names.tcnum;
END_SQL

    my $sth = $handle->prepare($sql);
    $sth->execute();
    while (my @row = $sth->fetchrow_array) {
        say $fh join "\t", @row;
    }
}

# Output uptake transporter report summary - only counts over thresholds per genome
sub report_uptake_summary {
    my ($handle,
        $fh) = @_;
    my %hash;
    # Prepare SQL queries
    my $sql_total = "SELECT genome, COUNT(*) FROM main WHERE pident >= ? AND qcovs >= ? GROUP BY genome";
    my $sql_totaltm = "SELECT genome, COUNT(*) FROM main WHERE pident >= ? AND qcovs >= ? AND tmdomains > 0 GROUP BY genome";
    my $sql_shortlist = << 'END_SQL';
SELECT genome, COUNT(*) FROM
                        (SELECT * FROM 
                                  (SELECT id, genome, tcdbacc, pident, qcovs, tmdomains, tcnum, name FROM main JOIN uptakefam_names ON main.tcfam=uptakefam_names.tcnum
                                  UNION ALL 
                                  SELECT id, genome, tcdbacc, pident, qcovs, tmdomains, tcnum, name FROM main JOIN uptakefam_names ON main.tcsubfam=uptakefam_names.tcnum
                                  UNION ALL 
                                  SELECT id, genome, tcdbacc, pident, qcovs, tmdomains, tcnum, name FROM main JOIN uptakefam_names ON main.tcsubsubfam=uptakefam_names.tcnum)
                        WHERE pident >= ? AND qcovs >= ?)
GROUP BY genome;
END_SQL

    my $sql_shortlisttm = << 'END_SQL';
SELECT genome, COUNT(*) FROM
                        (SELECT * FROM 
                                  (SELECT id, genome, tcdbacc, pident, qcovs, tmdomains, tcnum, name FROM main JOIN uptakefam_names ON main.tcfam=uptakefam_names.tcnum
                                  UNION ALL 
                                  SELECT id, genome, tcdbacc, pident, qcovs, tmdomains, tcnum, name FROM main JOIN uptakefam_names ON main.tcsubfam=uptakefam_names.tcnum
                                  UNION ALL 
                                  SELECT id, genome, tcdbacc, pident, qcovs, tmdomains, tcnum, name FROM main JOIN uptakefam_names ON main.tcsubsubfam=uptakefam_names.tcnum)
                        WHERE pident >= ? AND qcovs >= ? AND tmdomains > 0)
GROUP BY genome;
END_SQL

    # Execute queries and get values into hash
    my $sth1 = $handle->prepare($sql_total);
    $sth1->execute($blast_minpident,$blast_minqcovs);
    while (my @row = $sth1->fetchrow_array) {
        $hash{$row[0]}{'total'} = $row[1];
    }
    my $sth2 = $handle->prepare($sql_totaltm);
    $sth2->execute($blast_minpident,$blast_minqcovs);
    while (my @row = $sth2->fetchrow_array) {
        $hash{$row[0]}{'total_tm'} = $row[1];
    }
    my $sth3 = $handle->prepare($sql_shortlist);
    $sth3->execute($blast_minpident,$blast_minqcovs);
    while (my @row = $sth3->fetchrow_array) {
        $hash{$row[0]}{'shortlist'} = $row[1];
    }
    my $sth4 = $handle->prepare($sql_shortlisttm);
    $sth4->execute($blast_minpident,$blast_minqcovs);
    while (my @row = $sth4->fetchrow_array) {
        $hash{$row[0]}{'shortlist_tm'} = $row[1];
    }
    
    # Print header to output file
    say $fh join "\t", qw (genome total total_tm shortlist shortlist_tm);
    # Write report rows to output file
    foreach my $genome (sort {$a cmp $b} keys %hash) {
        my @out = ($genome);
        foreach my $value (qw(total total_tm shortlist shortlist_tm)) {
            if (defined $hash{$genome}{$value}) {
                push @out, $hash{$genome}{$value};
            } else {
                push @out, 0;
            }
        }
        say $fh join "\t", @out;
    }
}

# Output uptake transporter report summary - counts over thresholds per genome and transporter type
sub report_uptake_pergenome {
    my ($handle,
        $fh) = @_;
    my %hash;
    
    # Prepare SQL queries
    my $sql_countall = << 'END_SQL';
SELECT tcnum, genome, COUNT(*) FROM
                        (SELECT * FROM 
                                  (SELECT id, genome, tcdbacc, pident, qcovs, tmdomains, tcnum, name FROM main JOIN uptakefam_names ON main.tcfam=uptakefam_names.tcnum
                                  UNION ALL 
                                  SELECT id, genome, tcdbacc, pident, qcovs, tmdomains, tcnum, name FROM main JOIN uptakefam_names ON main.tcsubfam=uptakefam_names.tcnum
                                  UNION ALL 
                                  SELECT id, genome, tcdbacc, pident, qcovs, tmdomains, tcnum, name FROM main JOIN uptakefam_names ON main.tcsubsubfam=uptakefam_names.tcnum)
                        WHERE pident >= ? AND qcovs >= ? AND tmdomains >= 0)
GROUP BY tcnum,genome;
END_SQL

    my $sql_counttm = << 'END_SQL';
SELECT tcnum, genome, COUNT(*) FROM
                        (SELECT * FROM 
                                  (SELECT id, genome, tcdbacc, pident, qcovs, tmdomains, tcnum, name FROM main JOIN uptakefam_names ON main.tcfam=uptakefam_names.tcnum
                                  UNION ALL 
                                  SELECT id, genome, tcdbacc, pident, qcovs, tmdomains, tcnum, name FROM main JOIN uptakefam_names ON main.tcsubfam=uptakefam_names.tcnum
                                  UNION ALL 
                                  SELECT id, genome, tcdbacc, pident, qcovs, tmdomains, tcnum, name FROM main JOIN uptakefam_names ON main.tcsubsubfam=uptakefam_names.tcnum)
                        WHERE pident >= ? AND qcovs >= ? AND tmdomains > 0)
GROUP BY tcnum,genome;
END_SQL

    # Execute queries
    my $sth1 = $handle->prepare($sql_countall);
    $sth1->execute($blast_minpident, $blast_minqcovs);
    while (my @row = $sth1->fetchrow_array) {
        $hash{$row[0]}{$row[1]}{'all'} = $row[2];
    }
    my $sth2 = $handle->prepare($sql_counttm);
    $sth2->execute($blast_minpident, $blast_minqcovs);
    while (my @row = $sth2->fetchrow_array) {
        $hash{$row[0]}{$row[1]}{'tm'} = $row[2];
    }
    
    # Print header to report file
    say $fh join "\t", qw (tcnum genome counts counts_tm);
    # Write report rows to output file
    foreach my $tcnum (sort {$a cmp $b} keys %hash) {
        foreach my $genome (sort {$a cmp $b} keys %{$hash{$tcnum}}) {
            my @out = ($tcnum, $genome);
            if (defined $hash{$tcnum}{$genome}{'all'}) {
                push @out, $hash{$tcnum}{$genome}{'all'};
            } else {
                push @out, 0;
            }
            if (defined $hash{$tcnum}{$genome}{'tm'}) {
                push @out, $hash{$tcnum}{$genome}{'tm'};
            } else {
                push @out, 0;
            }
            say $fh join "\t", @out;
        }
    }
}

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2018 by Brandon Seah (kbseah@mpi-bremen.de)

LICENSE

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

=cut
