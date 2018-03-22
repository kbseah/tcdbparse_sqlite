# tcdbparse_sqlite

Perl utilities to manage data for transporter analysis of prokaryotic genomes, using Blastp searches vs. the [Transporter Classification Database](http://tcdb.org/) (TCDB).

## Usage

Requires: sqlite3, NCBI Blast+, TMHMM v2.0c, GNU Parallel.

Download Fasta file of TCDB sequence data from [their website](http://tcdb.org/download.php), and make Blast database with `makeblastdb`.

Prepare Fasta files of protein sequences, one file per genome. Make a tab-separated table of genome names and paths to each Fasta file.

Run Blastp and TMHMM with `do_blastp_tmhmm.pl`

Prepare database with `tcdbparse_sqlite.pl -task setup` using the template `template_database.sqlite`, and load data from Blastp and TMHMM results with argument `-task read`. You can view the SQLite database e.g. with [DB Browser for SQLite](http://sqlitebrowser.org/)

Write summary reports from the database with `tcdbparse_sqlite.pl -task report`. 

Details are given in help messages for each script, run with `--help` or `--man` arguments.

