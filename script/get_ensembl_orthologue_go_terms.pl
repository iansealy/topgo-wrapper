#!/usr/bin/env perl

# PODNAME: get_ensembl_orthologue_go_terms.pl
# ABSTRACT: Get GO terms for orthologues of all Ensembl genes

## Author         : is1
## Maintainer     : is1
## Created        : 2015-05-19
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
        -Ibranch-ensembl-74/ensembl-compara/modules \
        get_ensembl_orthologue_go_terms.pl \
        > data/orthologues/danio_rerio_homo_sapiens_e74_go.txt

=cut

# Default options
my $ensembl_species            = 'danio_rerio';
my $ensembl_orthologue_species = 'homo_sapiens';
my $ensembl_dbhost             = 'ensembldb.ensembl.org';
my $ensembl_dbport;
my $ensembl_dbuser = 'anonymous';
my $ensembl_dbpass;
my $slice_regexp;
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
my $sa =
  Bio::EnsEMBL::Registry->get_adaptor( $ensembl_species, 'core', 'Slice' );
my $goa = Bio::EnsEMBL::Registry->get_adaptor( 'Multi', 'Ontology', 'GOTerm' );
if ( !$goa ) {
    $goa = Bio::EnsEMBL::Registry->get_adaptor( 'Multi', 'Ontology',
        'OntologyTerm' );
}
my $gma =
  Bio::EnsEMBL::Registry->get_adaptor( 'Multi', 'Compara', 'GeneMember' );
my $ha = Bio::EnsEMBL::Registry->get_adaptor( 'Multi', 'Compara', 'Homology' );
my $gdba =
  Bio::EnsEMBL::Registry->get_adaptor( 'Multi', 'Compara', 'GenomeDB' );

# Get GenomeDB
my $gdb = $gdba->fetch_by_name_assembly($ensembl_species);

# Ensure database connection isn't lost; Ensembl 64+ can do this more elegantly
## no critic (ProhibitMagicNumbers)
if ( Bio::EnsEMBL::ApiVersion::software_version() < 64 ) {
## use critic
    Bio::EnsEMBL::Registry->set_disconnect_when_inactive();
}
else {
    Bio::EnsEMBL::Registry->set_reconnect_when_lost();
}

# Get all slices
my $slices = $sa->fetch_all('toplevel');
warn scalar @{$slices}, " slices\n" if $debug;

foreach my $slice ( @{$slices} ) {
    next
      if defined $slice_regexp
      && $slice->seq_region_name !~ m/$slice_regexp/xms;
    warn 'Slice: ', $slice->name, "\n" if $debug;

    # Get all genes
    my $genes = $slice->get_all_Genes( undef, 'core' );
    warn scalar @{$genes}, " genes\n" if $debug;
    foreach my $gene ( @{$genes} ) {
        warn ' Gene: ', $gene->stable_id, ' / ', $gene->biotype, "\n" if $debug;

        my $member;
        if ( $gma->can('fetch_by_stable_id') ) {
            $member = $gma->fetch_by_stable_id( $gene->stable_id );
        }
        elsif ( $gma->can('fetch_by_source_stable_id') ) {
            $member =
              $gma->fetch_by_source_stable_id( 'ENSEMBLGENE',
                $gene->stable_id );
        }
        elsif ( $gma->can('fetch_by_stable_id_GenomeDB') ) {
            $member =
              $gma->fetch_by_stable_id_GenomeDB( $gene->stable_id, $gdb );
        }
        if ( !$member ) {
            warn 'No member for ', $gene->stable_id, "\n";
            next;
        }

        my $homologies = $ha->fetch_all_by_Member( $member,
            -TARGET_SPECIES => $ensembl_orthologue_species );
        foreach my $homology ( @{$homologies} ) {
            foreach my $homology_member ( @{ $homology->gene_list() } ) {
                next if $homology_member->stable_id eq $gene->stable_id;
                my $links = $homology_member->get_Gene()->get_all_DBLinks();
                foreach my $link ( @{$links} ) {
                    next if ( ref $link ) !~ m/OntologyXref/xms;
                    my $term = $goa->fetch_by_accession( $link->primary_id );
                    next if !$term || !$term->namespace;
                    print join "\t", $gene->stable_id, $term->accession,
                      $term->namespace;
                    print "\n";
                }
            }
        }
    }
}

# Get and check command line options
sub get_and_check_options {

    # Get options
    GetOptions(
        'ensembl_species=s'            => \$ensembl_species,
        'ensembl_orthologue_species=s' => \$ensembl_orthologue_species,
        'ensembl_dbhost=s'             => \$ensembl_dbhost,
        'ensembl_dbport=i'             => \$ensembl_dbport,
        'ensembl_dbuser=s'             => \$ensembl_dbuser,
        'ensembl_dbpass=s'             => \$ensembl_dbpass,
        'slice_regexp=s'               => \$slice_regexp,
        'debug'                        => \$debug,
        'help'                         => \$help,
        'man'                          => \$man,
    ) or pod2usage(2);

    # Documentation
    if ($help) {
        pod2usage(1);
    }
    elsif ($man) {
        pod2usage( -verbose => 2 );
    }

    return;
}

=head1 USAGE

    get_ensembl_go_terms.pl
        [--ensembl_species species]
        [--ensembl_orthologue_species species]
        [--ensembl_dbhost host]
        [--ensembl_dbport port]
        [--ensembl_dbuser username]
        [--ensembl_dbpass password]
        [--slice_regexp regexp]
        [--debug]
        [--help]
        [--man]

=head1 OPTIONS

=over 8

=item B<--ensembl_species SPECIES>

Ensembl species (e.g. danio_rerio).

=item B<--ensembl_orthologue_species SPECIES>

Ensembl species from which to get orthologues (e.g. homo_sapiens).

=item B<--ensembl_dbhost HOST>

Ensembl MySQL database host.

=item B<--ensembl_dbport PORT>

Ensembl MySQL database port.

=item B<--ensembl_dbuser USERNAME>

Ensembl MySQL database username.

=item B<--ensembl_dbpass PASSWORD>

Ensembl MySQL database password.

=item B<--slice_regexp REGEXP>

Regular expression for limiting slices.

=item B<--debug>

Print debugging information.

=item B<--help>

Print a brief help message and exit.

=item B<--man>

Print this script's manual page and exit.

=back

=cut
