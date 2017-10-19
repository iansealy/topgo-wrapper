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
Readonly our @ALGORITHM => qw( classic elim weight weight01 lea parentchild );
Readonly our %ALGORITHM => map { $_ => 1 } @ALGORITHM;

# Default options
my $dir = q{.};
my $input_file;
my $detct_file;
my $algorithm = 'elim';
my $genes_of_interest_file;
my $gene_field = 1;
my $p_value_field;
my $fold_change_field;
my $name_field;
my $description_field;
my $sig_level = 0.05;
my $input_sig_level;
my $output_sig_level;
my $go_terms_file;
my $r_binary   = 'R';
my $has_header = 0;
my ( $debug, $help, $man );

# Get and check command line options
get_and_check_options();

# Use default significance level for both, unless specifically specified
$input_sig_level  = $input_sig_level  || $sig_level;
$output_sig_level = $output_sig_level || $sig_level;

# Automatically configure for DETCT output
if ($detct_file) {

    # Set defaults for DETCT all.tsv file
    $input_file = $detct_file;
    ## no critic (ProhibitMagicNumbers)
    $gene_field        = 10;
    $p_value_field     = 8;
    $name_field        = 14;
    $description_field = 15;
    ## use critic
    $has_header = 1;

    # Get fold change field, Ensembl species and Ensembl version from input
    my $species;
    my %prefix2species = (
        ''    => 'homo_sapiens',
        'MUS' => 'mus_musculus',
        'DAR' => 'danio_rerio',
    );
    open my $fh, '<', $detct_file;
    my $header = <$fh>;
    while ( my $line = <$fh> ) {
        my @fields = split /\t/xms, $line;
        if ( $fields[9] =~ m/\A ENS(|MUS|DAR)G\d/xms ) {
            $species = $prefix2species{$1};
            last;
        }
    }
    close $fh;
    my @fields = split /\t/xms, $header;
    my ($ens) = $fields[9] =~ m/\A (e\d+) \s/xms;    # e.g. e76 Ensembl Gene ID
    foreach my $i ( 1 .. scalar @fields ) {
        if ( $fields[ $i - 1 ] =~ m/\A Log2 \s fold \s change/xms ) {
            $fold_change_field = $i;
            last;
        }
    }

    if ( !$go_terms_file ) {
        my $filename = sprintf '%s_%s_go.txt', $species, $ens;
        $go_terms_file =
          File::Spec->catfile( dirname(__FILE__), File::Spec->updir(), 'data',
            $filename );
    }
}

# Ensure working directory exists
make_path($dir);

# Get p values and other data for every Ensembl gene in the gene universe
my ( $p_value_for, $fold_change_for, $name_for, $description_for ) =
  TopGO::read_gene_info( $input_file, $genes_of_interest_file, $has_header,
    $gene_field, $p_value_field, $fold_change_field, $name_field,
    $description_field );

# Get GO terms for every Ensembl gene
my $go_terms_for = TopGO::read_go_terms($go_terms_file);

# Write gene list
my @sets = ('all');    # All genes
my ( %gene_list_file, %dir );
$gene_list_file{'all'} = File::Spec->catfile( $dir, 'gene_list.txt' );
$dir{'all'} = $dir;
TopGO::write_gene_list( $p_value_for, $gene_list_file{'all'} );

# If fold changes are available and, if genes of interest, not all are up or
# downregulated, also run up and downregulated subsets
my %sig_genes_for;
my $run_up_and_down = $fold_change_field ? 1 : 0;
if ($genes_of_interest_file) {
    @{ $sig_genes_for{'all'} } =
      grep { $p_value_for->{$_} == 1 } keys %{$p_value_for};
    @{ $sig_genes_for{'up'} } =
      grep { $fold_change_for->{$_} ne q{-} && $fold_change_for->{$_} > 0 }
      @{ $sig_genes_for{'all'} };
    @{ $sig_genes_for{'down'} } =
      grep { $fold_change_for->{$_} ne q{-} && $fold_change_for->{$_} < 0 }
      @{ $sig_genes_for{'all'} };
    if ( !@{ $sig_genes_for{'up'} } || !@{ $sig_genes_for{'down'} } ) {
        $run_up_and_down = 0;
    }
}
if ($run_up_and_down) {
    push @sets, 'up', 'down';

    $gene_list_file{'up'} = File::Spec->catfile( $dir, 'up', 'gene_list.txt' );
    $gene_list_file{'down'} =
      File::Spec->catfile( $dir, 'down', 'gene_list.txt' );

    $dir{'up'}   = File::Spec->catdir( $dir, 'up' );
    $dir{'down'} = File::Spec->catdir( $dir, 'down' );

    make_path( $dir{'up'} );
    make_path( $dir{'down'} );

    my @up_genes =
      grep { $fold_change_for->{$_} ne q{-} && $fold_change_for->{$_} > 0 }
      keys %{$fold_change_for};
    my @down_genes =
      grep { $fold_change_for->{$_} ne q{-} && $fold_change_for->{$_} < 0 }
      keys %{$fold_change_for};

    TopGO::write_gene_list( $p_value_for, $gene_list_file{'up'}, \@up_genes );
    TopGO::write_gene_list( $p_value_for, $gene_list_file{'down'},
        \@down_genes );
}

# Remove sets where genes of interest have no GO terms
if ($genes_of_interest_file) {
    my @extant_sets = @sets;
    foreach my $set (@extant_sets) {
        my $got_go_terms = 0;
        foreach my $sig_gene ( @{ $sig_genes_for{$set} } ) {
            if ( exists $go_terms_for->{$sig_gene} ) {
                $got_go_terms = 1;
                last;
            }
        }
        if ( !$got_go_terms ) {
            @sets = grep { $_ ne $set } @sets;
        }
    }
}

# Write mapping file and run topGO for each domain
my $topgo_script = File::Spec->catfile( dirname(__FILE__), 'run_topgo_ks.R' );
if ($genes_of_interest_file) {
    $topgo_script =
      File::Spec->catfile( dirname(__FILE__), 'run_topgo_fisher.R' );
}
foreach my $domain ( sort keys %DOMAIN ) {
    my $gene_to_go_mapping_file =
      File::Spec->catfile( $dir,
        ( sprintf '%s_all.gene2go.txt', $DOMAIN{$domain} ) );
    TopGO::write_mapping_file( [ sort keys %{$p_value_for} ],
        $go_terms_for, $DOMAIN{$domain}, $gene_to_go_mapping_file );
    foreach my $set (@sets) {
        my $output_prefix = File::Spec->catfile( $dir{$set}, $domain );
        TopGO::run_topgo(
            {
                r_binary         => $r_binary,
                topgo_script     => $topgo_script,
                algorithm        => $algorithm,
                gene_list_file   => $gene_list_file{$set},
                mapping_file     => $gene_to_go_mapping_file,
                domain           => $domain,
                output_prefix    => $output_prefix,
                input_sig_level  => $input_sig_level,
                output_sig_level => $output_sig_level,
            }
        );
        TopGO::annotate_with_genes(
            {
                input_file   => $output_prefix . '.all.tsv',
                output_file  => $output_prefix . '.all.genes.tsv',
                p_values     => $p_value_for,
                fold_changes => $fold_change_for,
                names        => $name_for,
                descriptions => $description_for,
            }
        );
        TopGO::filter_by_significance(
            $output_prefix . '.all.tsv',
            $output_prefix . '.sig.tsv',
            $output_sig_level
        );
        TopGO::filter_by_significance(
            $output_prefix . '.all.genes.tsv',
            $output_prefix . '.sig.genes.tsv',
            $output_sig_level
        );
    }
}

# Get and check command line options
sub get_and_check_options {

    # Get options
    GetOptions(
        'dir=s'                    => \$dir,
        'input_file=s'             => \$input_file,
        'detct_file=s'             => \$detct_file,
        'algorithm=s'              => \$algorithm,
        'genes_of_interest_file=s' => \$genes_of_interest_file,
        'gene_field=i'             => \$gene_field,
        'p_value_field=i'          => \$p_value_field,
        'fold_change_field=i'      => \$fold_change_field,
        'name_field=i'             => \$name_field,
        'description_field=i'      => \$description_field,
        'sig_level=f'              => \$sig_level,
        'input_sig_level=f'        => \$input_sig_level,
        'output_sig_level=f'       => \$output_sig_level,
        'go_terms_file=s'          => \$go_terms_file,
        'r_binary=s'               => \$r_binary,
        'header'                   => \$has_header,
        'debug'                    => \$debug,
        'help'                     => \$help,
        'man'                      => \$man,
    ) or pod2usage(2);

    # Documentation
    if ($help) {
        pod2usage(1);
    }
    elsif ($man) {
        pod2usage( -verbose => 2 );
    }

    if ( !$input_file && !$detct_file ) {
        pod2usage("--input_file or --detct_file must be specified\n");
    }
    if ( $input_file && $detct_file ) {
        pod2usage("--input_file and --detct_file must not both be specified\n");
    }
    if ( !$go_terms_file && !$detct_file ) {
        pod2usage("--go_terms_file must be specified\n");
    }

    if ( !$p_value_field && !$genes_of_interest_file ) {
        $p_value_field = 3;    # Default
    }

    if ( !exists $ALGORITHM{$algorithm} ) {
        pod2usage( sprintf "--algorithm must be one of %s\n",
            join ', ', @ALGORITHM );
    }

    return;
}

=head1 USAGE

    run_topgo.pl
        [--dir directory]
        [--input_file file]
        [--detct_file file]
        [--algorithm string]
        [--genes_of_interest_file file]
        [--gene_field int]
        [--p_value_field int]
        [--fold_change_field int]
        [--name_field int]
        [--description_field int]
        [--sig_level float]
        [--input_sig_level float]
        [--output_sig_level float]
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

The tab-separated file containing Ensembl gene IDs and (optional) associated p
values.

=item B<--detct_file FILE>

The all.tsv file output by the DETCT pipeline. If this option is specified then
any field options are ignored.

=item B<--algorithm ALGORITHM>

The algorithm used by topGO. Default to elim. Can also be classic, weight,
weight01, lea or parentchild.

=item B<--genes_of_interest_file FILE>

The file containing Ensembl gene IDs of interest (in which case p values aren't
needed).

=item B<--gene_field INT>

The field that specifies the Ensembl gene ID in the input file.

=item B<--p_value_field INT>

The field that specifies the p value in the input file.

=item B<--fold_change_field INT>

The field that specifies the gene's log fold change in the input file.

=item B<--name_field INT>

The field that specifies the gene's name in the input file.

=item B<--description_field INT>

The field that specifies the gene's description in the input file.

=item B<--sig_level FLOAT>

The level at which p values are considered significant.

=item B<--input_sig_level FLOAT>

The level at which input p values are considered significant.

=item B<--output_sig_level FLOAT>

The level at which output (i.e. topGO) p values are considered significant.

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
