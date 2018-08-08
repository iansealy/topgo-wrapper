#!/usr/bin/env perl

# PODNAME: get_go_terms.pl
# ABSTRACT: Get GO terms from an Ensembl database

## Author         : is1
## Maintainer     : is1
## Created        : 2014-11-13
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
use Bio::EnsEMBL::Registry;

=head1 DESCRIPTION


=head1 EXAMPLES

    perl \
        -Ibranch-ensembl-74/ensembl/modules \
        get_go_terms.pl \
        > data/e74_go.txt

    perl \
        -Ibranch-ensembl-74/ensembl/modules \
        get_go_terms.pl \
        --parent_accession GO:0006412 \
        > data/e74_go.txt

=cut

# Constants
Readonly our @ROOT_ACCESSIONS => ( 'GO:0008150', 'GO:0003674', 'GO:0005575' );

# Default options
my @parent_accessions;
my $ensembl_dbhost = 'ensembldb.ensembl.org';
my $ensembl_dbport;
my $ensembl_dbuser = 'anonymous';
my $ensembl_dbpass;
my ( $debug, $help, $man );

# Get and check command line options
get_and_check_options();

# Connnect to Ensembl database
Bio::EnsEMBL::Registry->load_registry_from_db(
    -host => $ensembl_dbhost,
    -port => $ensembl_dbport,
    -user => $ensembl_dbuser,
    -pass => $ensembl_dbpass,
);

# Get genebuild version
my $genebuild_version = 'e' . Bio::EnsEMBL::ApiVersion::software_version();
warn 'Genebuild version: ', $genebuild_version, "\n" if $debug;

# Get Ensembl adaptors
my $goa = Bio::EnsEMBL::Registry->get_adaptor( 'Multi', 'Ontology', 'GOTerm' );
if ( !$goa ) {
    $goa = Bio::EnsEMBL::Registry->get_adaptor( 'Multi', 'Ontology',
        'OntologyTerm' );
}

# Ensure database connection isn't lost; Ensembl 64+ can do this more elegantly
## no critic (ProhibitMagicNumbers)
if ( Bio::EnsEMBL::ApiVersion::software_version() < 64 ) {
## use critic
    Bio::EnsEMBL::Registry->set_disconnect_when_inactive();
}
else {
    Bio::EnsEMBL::Registry->set_reconnect_when_lost();
}

# Get all terms below parent term(s)
my %name_for;
my %namespace_for;
foreach my $parent_accession (@parent_accessions) {
    my $term = $goa->fetch_by_accession($parent_accession);
    get_child_terms($term);
}

# Output terms
foreach my $accession ( sort keys %name_for ) {
    printf "%s\t%s\t%s\n", $accession, $namespace_for{$accession},
      $name_for{$accession};
}

# Recursive function for getting child terms
sub get_child_terms {
    my ($term) = @_;

    return if exists $name_for{ $term->accession };

    warn $term->accession, "\n" if $debug;

    $name_for{ $term->accession }      = $term->name;
    $namespace_for{ $term->accession } = $term->namespace;

    my $child_terms = $term->children;
    foreach my $child_term ( @{$child_terms} ) {
        get_child_terms($child_term);
    }

    return;
}

# Get and check command line options
sub get_and_check_options {

    # Get options
    GetOptions(
        'parent_accession=s' => \@parent_accessions,
        'ensembl_dbhost=s'   => \$ensembl_dbhost,
        'ensembl_dbport=i'   => \$ensembl_dbport,
        'ensembl_dbuser=s'   => \$ensembl_dbuser,
        'ensembl_dbpass=s'   => \$ensembl_dbpass,
        'debug'              => \$debug,
        'help'               => \$help,
        'man'                => \$man,
    ) or pod2usage(2);

    # Documentation
    if ($help) {
        pod2usage(1);
    }
    elsif ($man) {
        pod2usage( -verbose => 2 );
    }

    if ( !@parent_accessions ) {
        @parent_accessions = @ROOT_ACCESSIONS;
    }

    return;
}

=head1 USAGE

    get_go_terms.pl
        [--parent_accession accession...]
        [--ensembl_dbhost host]
        [--ensembl_dbport port]
        [--ensembl_dbuser username]
        [--ensembl_dbpass password]
        [--debug]
        [--help]
        [--man]

=head1 OPTIONS

=over 8

=item B<--parent_accession ACCESSION...>

GO term(s) to begin query with.

=item B<--ensembl_dbhost HOST>

Ensembl MySQL database host.

=item B<--ensembl_dbport PORT>

Ensembl MySQL database port.

=item B<--ensembl_dbuser USERNAME>

Ensembl MySQL database username.

=item B<--ensembl_dbpass PASSWORD>

Ensembl MySQL database password.

=item B<--debug>

Print debugging information.

=item B<--help>

Print a brief help message and exit.

=item B<--man>

Print this script's manual page and exit.

=back

=cut
