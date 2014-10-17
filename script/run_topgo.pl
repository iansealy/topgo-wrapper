#!/usr/bin/env perl

# PODNAME: run_topgo.pl
# ABSTRACT: Run topGO

## Author         : is1
## Maintainer     : is1
## Created        : 2014-10-17
## Last commit by : $Author$
## Last modified  : $Date$
## Revision       : $Revision$
## Repository URL : $HeadURL$

use warnings;
use strict;
use autodie;
use Carp;
use Try::Tiny;

use Getopt::Long;
use Pod::Usage;
use Readonly;
use TopGO;
use File::Spec;
use File::Path qw( make_path );
use File::Basename;

=head1 DESCRIPTION


=head1 EXAMPLES


=cut

# Constants
Readonly our %DOMAIN => (
    BP => 'biological_process',
    CC => 'cellular_component',
    MF => 'molecular_function',
);

# Default options
my $dir = q{.};
my $input_file;
my $gene_field    = 1;
my $p_value_field = 2;
my $sig_level     = 0.05;
my $go_terms_file;
my $r_binary   = 'R';
my $has_header = 0;
my ( $debug, $help, $man );

# Get and check command line options
get_and_check_options();

# Ensure working directory exists
make_path($dir);

# Get p values for every Ensembl gene in the gene universe
my $p_value_for =
  TopGO::read_p_values( $input_file, $has_header, $gene_field, $p_value_field );

# Get GO terms for every Ensembl gene
my $go_terms_for = TopGO::read_go_terms($go_terms_file);

# Write gene list
my $gene_list_file = File::Spec->catfile( $dir, 'gene_list.txt' );
TopGO::write_gene_list( $p_value_for, $gene_list_file );

# Write mapping file and run topGO for each domain
my $topgo_script = File::Spec->catfile( dirname(__FILE__), 'run_topgo.R' );
foreach my $domain ( keys %DOMAIN ) {
    my $gene_to_go_mapping_file =
      File::Spec->catfile( $dir,
        ( sprintf '%s_all.gene2go.txt', $DOMAIN{$domain} ) );
    TopGO::write_mapping_file( [ keys %{$p_value_for} ],
        $go_terms_for, $DOMAIN{$domain}, $gene_to_go_mapping_file );
    my $output_prefix =
      File::Spec->catfile( $dir, ( sprintf '%s_all', $domain ) );
    TopGO::run_topgo(
        {
            r_binary       => $r_binary,
            topgo_script   => $topgo_script,
            gene_list_file => $gene_list_file,
            mapping_file   => $gene_to_go_mapping_file,
            domain         => $domain,
            output_prefix  => $output_prefix,
            sig_level      => $sig_level,
        }
    );
}

# Get and check command line options
sub get_and_check_options {

    # Get options
    GetOptions(
        'dir=s'           => \$dir,
        'input_file=s'    => \$input_file,
        'gene_field=i'    => \$gene_field,
        'p_value_field=i' => \$p_value_field,
        'sig_level=f'     => \$sig_level,
        'go_terms_file=s' => \$go_terms_file,
        'r_binary=s'      => \$r_binary,
        'header'          => \$has_header,
        'debug'           => \$debug,
        'help'            => \$help,
        'man'             => \$man,
    ) or pod2usage(2);

    # Documentation
    if ($help) {
        pod2usage(1);
    }
    elsif ($man) {
        pod2usage( -verbose => 2 );
    }

    if ( !$input_file ) {
        pod2usage("--input_file must be specified\n");
    }
    if ( !$go_terms_file ) {
        pod2usage("--go_terms_file must be specified\n");
    }

    return;
}

=head1 USAGE

    run_topgo.pl
        [--dir directory]
        [--input_file file]
        [--gene_field int]
        [--p_value_field int]
        [--sig_level float]
        [--go_terms_file file]
        [--r_binary file]
        [--header]
        [--debug]
        [--help]
        [--man]

=head1 OPTIONS

=over 8

=item B<--dir DIRECTORY>

Working directory.

=item B<--input_file FILE>

The tab-separated file containing Ensembl gene IDs and associated p values.

=item B<--gene_field INT>

The field that specifies the Ensembl gene ID in the input file.

=item B<--p_value_field INT>

The field that specifies the p value in the input file.

=item B<--sig_level FLOAT>

The level at which p values are considered significant.

=item B<--go_terms_file FILE>

The file containing GO terms for all Ensembl genes (produced by
get_ensembl_go_terms.pl).

=item B<--r_binary FILE>

The R binary to run topGO with.

=item B<--header>

The input file has a header line.

=item B<--debug>

Print debugging information.

=item B<--help>

Print a brief help message and exit.

=item B<--man>

Print this script's manual page and exit.

=back

=cut
