#!/usr/bin/env perl

=head1 NAME

do_blastp_tmhmm.pl - Run Blastp vs. TCDB and TMHMM for protein sequences

=head1 SYNOPSIS

perl do_blastp_tmhmm.pl --faa <fasta file list> --blast <blast filename list> --tmhmm <tmmhmm filename list> --outdir <output dir>

perl do_blastp_tmhmm.pl --help

perl do_blastp_tmhmm.pl --man

=head1 DESCRIPTION

Run parallelized Blastp vs TCDB database, and TMHMM on Fasta files of protein
sequences from genomes, to compare their transporter content with the script
tcdbparse_sqlite.pl. 

Requries blastp (NCBI Blast+), tmhmm v2.0, and GNU parallel

=cut

use strict;
use warnings;

use Getopt::Long;
use File::Basename;
use File::Spec;
use Pod::Usage;

my $faa_list;
my $blast_list;
my $tmhmm_list;
my $blast_eval = "1e-5";
my $blast_num_hits = 1;
my $num_cores = 8;
my $blast_suffix = "_v_tcdb.blastp.out6";
my $outdir = "./faa";

# Paths to local DBs and binaries - replace with your own
#my $tcdb_blastdb_path = "/data/db/tcdb/tcdb";
my $tcdb_blastdb_path = "./tcdb";
my ($blastp_bin,$tmhmm_bin,$parallel_bin);

# Hashes for input files
my %lib_faa_hash;

pod2usage(verbose=>0,exit_status=>1) if !@ARGV;

GetOptions("faa=s"=>\$faa_list,
           "blast=s"=>\$blast_list,
           "outdir=s"=>\$outdir,
           "tmhmm=s"=>\$tmhmm_list,
           "blastdb=s"=>\$tcdb_blastdb_path,
           "eval=s"=>\$blast_eval,
           "hits=i"=>\$blast_num_hits,
           "proc|cores=i"=>\$num_cores,
           "blastoutsuffix=s"=>\$blast_suffix,
           "blastp_path=s"=>\$blastp_bin,
           "tmhmm_path=s"=>\$tmhmm_bin,
           "parallel_path=s"=>\$parallel_bin,
           "help|h"=>sub{pod2usage(verbose=>1,exit_status=>1);},
           "man|m"=>sub{pod2usage(verbose=>2,exit_status=>1);},
            );

=head1 ARGUMENTS

=over 8

=item --faa <file>

Tab-separated table of genome names (column 1) and paths to Fasta files with 
corresponding protein sequences (column 2).

=item --blast <file>

File to write the table of paths to Blastp output files. The output is a tab-
separated table of genome names (column 1) and paths to output filenames for
Blast results (column 2)

=item --tmhmm <file>

File to write the table of paths to TMHMM output files. The output is a tab-
separated table of genome names (column 1) and paths to output filenames for 
TMHMM results (column 2)

=item --outdir <path>

Folder for output files (Default: ./faa)

=item --blastdb <path>

Path to TCDB blast database

=item --eval <numeric>

Max E-value for Blastp hit to be counted (Default: 1e-5)

=item --hits <integer>

Number of hits per Blast query to report (Default: 1)

=item --blastoutsuffix <string>

Filename suffix for Blastp output files (Default: _v_tcdb.blastp.out6)

=item --proc <integer>

Number of processors to use for parallelized Blastp (Default: 8)

=item --blastdb <path>

Path to TCDB Blast database (including DB filename prefix) (Default: ./tcdb)

=item --blastp_path <path>

Path to Blastp binary (including filename)

=item --tmhmm_path <path>

Path to TMHMM binary (including filename)

=item --parallel_path <path>

Path to GNU parallel binary (including filename)

=back

=cut

my $ref1 = hashTSV_KV($faa_list);
%lib_faa_hash = %$ref1;

my %lib_blastout_hash;
my %lib_tmhmmout_hash;

foreach my $lib (sort {$a cmp $b} keys %lib_faa_hash) {
    my $faa = $lib_faa_hash{$lib};
    my ($faafile,$faadirs,$faasuffix) = fileparse($faa,".faa");
    #my $blastout = $faadirs.$faafile.$blast_suffix;
    my $blastoutfile = $faafile.$blast_suffix;
    my $blastout = File::Spec->catfile ($outdir, $blastoutfile);
    #my $tmhmmout = $faadirs.$faafile."_tmhmm.out";
    my $tmhmmoutfile = $faafile."_tmhmm.out";
    my $tmhmmout = File::Spec->catfile ($outdir, $tmhmmoutfile);
    if (! -f $blastout) {
        # do blast
        # blastp -db /data/db/tcdb/tcdb_2015 -query $i -outfmt 6 -evalue 1e-5 -max_target_seqs 10 -out $l.vs_tcdb.blastp.out6
        print STDERR "Performing Blastp, results to file $blastout\n";
        my @blastparams = ("-db $tcdb_blastdb_path",
                           "-query -",
                           "-outfmt \\\'6 std qcovs\\\'", # Escape quote marks to work with GNU parallel
                           "-evalue $blast_eval",
                           "-max_target_seqs $blast_num_hits");
        my @parallelparams = ("--gnu",
                              "--no-notice",
                              "-j $num_cores",
                              "--recstart \'>\'",
                              "-N 100",
                              "--pipe");
        my $blastcommand = join " ", ("cat $faa \|",
                                      $parallel_bin,
                                      @parallelparams,
                                      $blastp_bin,
                                      @blastparams,
                                      "> $blastout");
        #my @blastparams = ("-db $tcdb_blastdb_path",
        #                   "-query $faa",
        #                   "-outfmt \'6 std qcovs\'",
        #                   "-evalue $blast_eval",
        #                   "-max_target_seqs $blast_num_hits",
        #                   "-out $blastout");
        #my $blastcommand = join " ", ($blastp_bin, @blastparams);
        my $exitblast = system ($blastcommand);
        $lib_blastout_hash{$lib} = $blastout if ($exitblast == 0);

    } else {
        print STDERR "Blast result file $blastout already exists, skipping...\n";
    }
    if (! -f $tmhmmout) {
        # cat $i | sed 's/\*//g' | tmhmm > $l.tmhmm_results
        print STDERR "Performing TMHMM, results to file $tmhmmout\n";
        my $tmhmmcommand = ("cat $faa | $tmhmm_bin > $tmhmmout");
        my $exittmhmm = system ($tmhmmcommand);
        $lib_tmhmmout_hash{$lib} = $tmhmmout if ($exittmhmm == 0);
        system("rm -r TMHMM_*"); # Remove bulky TMHMM output folder and keep only summary

    } else {
        print STDERR "TMHMM result file $tmhmmout already exists, skipping... \n";
    }
}

# Write to log files

open(BLOUT, ">>", $blast_list) or die ("Cannot open $blast_list: $!");
foreach my $thelib (sort {$a cmp $b} keys %lib_blastout_hash) {
    print BLOUT join "\t", ($thelib,$lib_blastout_hash{$thelib});
    print BLOUT "\n";
}
close(BLOUT);

open(TMOUT, ">>", $tmhmm_list) or die ("Cannot open $tmhmm_list:$!");
foreach my $thelib (sort {$a cmp $b} keys %lib_tmhmmout_hash) {
    print TMOUT join "\t", ($thelib, $lib_tmhmmout_hash{$thelib});
    print TMOUT "\n";
}
close(TMOUT);


# Open TSV file, split cols by tab, and hash with col1 as key and col2 as val
sub hashTSV_KV {
    my ($file) = @_;
    my %hash;
    open(IN, "<", $file) or die ("File $file not found: $!");
    while (<IN>) {
        chomp;
        my @splitline = split "\t";
        $hash{$splitline[0]} = $splitline[1];
    }
    close(IN);
    return (\%hash);
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
