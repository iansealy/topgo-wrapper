## no critic (RequireUseStrict, RequireUseWarnings, RequireTidyCode)
package TopGO;
## use critic

# ABSTRACT: TopGO wrapper

## Author         : is1
## Maintainer     : is1
## Created        : 2014-10-16
## Last commit by : $Author$
## Last modified  : $Date$
## Revision       : $Revision$
## Repository URL : $HeadURL$

use warnings;
use strict;
use autodie;
use Carp;
use Try::Tiny;

use English qw( -no_match_vars );
use POSIX qw( WIFEXITED);

use base qw( Exporter );
our @EXPORT_OK = qw(
  read_p_values
  read_go_terms
  write_gene_list
  write_mapping_file
);

=head1 SYNOPSIS

    # Brief code examples

=cut

=func read_p_values

  Usage       : $p_value_for = read_p_values( $input_file, $has_header,
                    $gene_field, $p_value_field );
  Purpose     : Read p values from a file and organise by gene ID
  Returns     : Hashref {
                    String (gene ID) => Float (p value)
                }
  Parameters  : String (the input filename)
                Boolean (whether input file has header)
                Int (the gene field)
                Int (the p value field)
  Throws      : No exceptions
  Comments    : None

=cut

sub read_p_values {
    my ( $input_file, $has_header, $gene_field, $p_value_field ) = @_;

    my %p_value_for;

    open my $fh, '<', $input_file;
    if ($has_header) {
        <$fh>;
    }
    while ( my $line = <$fh> ) {
        chomp $line;
        my @fields = split /\t/xms, $line;
        my $gene_ids = $fields[ $gene_field - 1 ];
        next if $gene_ids eq q{-};
        my $p_value = $fields[ $p_value_field - 1 ];
        next if $p_value eq q{-} || $p_value eq 'NA';    # Ignore filtered genes
        foreach my $gene_id ( split /,/xms, $gene_ids ) {
            confess "Multiple p values for $gene_id"
              if exists $p_value_for{$gene_id};
            $p_value_for{$gene_id} = $p_value;
        }
    }
    close $fh;

    return \%p_value_for;
}

=func read_go_terms

  Usage       : $go_terms_for = read_go_terms( $go_terms_file );
  Purpose     : Read GO terms from a file and organise by gene ID and domain
  Returns     : Hashref {
                    String (gene ID) => Hashref {
                        String (domain) => Arrayref of strings (the GO terms)
                    }
                }
  Parameters  : String (the input filename)
  Throws      : No exceptions
  Comments    : None

=cut

sub read_go_terms {
    my ($input_file) = @_;

    my %go_terms_for;

    open my $fh, '<', $input_file;
    while ( my $line = <$fh> ) {
        chomp $line;
        my ( $gene_id, $go_term, $domain ) = split /\t/xms, $line;
        push @{ $go_terms_for{$gene_id}{$domain} }, $go_term;
    }
    close $fh;

    return \%go_terms_for;
}

=func write_gene_list

  Usage       : write_gene_list($p_value_for, $gene_list_file);
  Purpose     : Write gene list file
  Returns     : undef
  Parameters  : Hashref (of p values keyed by gene ID)
                String (the output filename)
  Throws      : No exceptions
  Comments    : None

=cut

sub write_gene_list {
    my ( $p_value_for, $gene_list_file ) = @_;

    open my $fh, '>', $gene_list_file;
    foreach my $gene_id ( keys %{$p_value_for} ) {
        printf {$fh} "%s\t%f\n", $gene_id, $p_value_for->{$gene_id};
    }
    close $fh;

    return;
}

=func write_mapping_file

  Usage       : write_mapping_file($gene_ids, $go_terms_for, $domain,
                    $gene_to_go_mapping_file);
  Purpose     : Write gene-to-GOs mapping file
  Returns     : undef
  Parameters  : Arrayref (of gene IDs)
                Hashref (of GO terms keyed by gene ID and domain)
                String (the domain)
                String (the output filename)
  Throws      : No exceptions
  Comments    : None

=cut

sub write_mapping_file {
    my ( $gene_ids, $go_terms_for, $domain, $gene_to_go_mapping_file ) = @_;

    open my $fh, '>', $gene_to_go_mapping_file;
    foreach my $gene_id ( @{$gene_ids} ) {
        if ( exists $go_terms_for->{$gene_id}->{$domain} ) {
            printf {$fh} "%s\t%s\n", $gene_id,
              ( join ', ', @{ $go_terms_for->{$gene_id}->{$domain} } );
        }
    }
    close $fh;

    return;
}

=func run_topgo

  Usage       : TopGO::run_topgo( {
                    r_binary       => 'R',
                    topgo_script   => 'script/run_topgo.R',
                    gene_list_file => 'gene_list.txt',
                    mapping_file   => 'gene2go.txt',
                    domain         => 'BP',
                    output_prefix  => 'BP_all',
                    sig_level      => 0.05,
                } );
  Purpose     : Run topGO
  Returns     : undef
  Parameters  : Hashref {
                    r_binary       => String (the R binary),
                    topgo_script   => String (the topGO script),
                    gene_list_file => String (the gene list file),
                    mapping_file   => String (the mapping file),
                    domain         => String (the domain),
                    output_prefix  => String (the output prefix),
                    sig_level      => Float (significance level threshold),
                }
  Throws      : If R binary is missing
                If topGO script is missing
                If gene list file is missing
                If mapping file is missing
                If domain is missing
                If output prefix is missing
                If significance level is missing
                If command line can't be run
  Comments    : None

=cut

sub run_topgo {
    my ($arg_ref) = @_;

    confess 'No R binary specified'     if !defined $arg_ref->{r_binary};
    confess 'No topGO script specified' if !defined $arg_ref->{topgo_script};
    confess 'No gene list file specified'
      if !defined $arg_ref->{gene_list_file};
    confess 'No mapping file specified'  if !defined $arg_ref->{mapping_file};
    confess 'No domain specified'        if !defined $arg_ref->{domain};
    confess 'No output prefix specified' if !defined $arg_ref->{output_prefix};
    confess 'No significance level specified' if !defined $arg_ref->{sig_level};

    my $output_file = $arg_ref->{output_prefix} . '.tsv';
    my $stdout_file = $arg_ref->{output_prefix} . '.o';
    my $stderr_file = $arg_ref->{output_prefix} . '.e';

    my $cmd = join q{ }, $arg_ref->{r_binary}, '--slave', '--args',
      $arg_ref->{gene_list_file}, $arg_ref->{mapping_file}, $arg_ref->{domain},
      $output_file, $arg_ref->{output_prefix}, $arg_ref->{sig_level}, '<',
      $arg_ref->{topgo_script};
    $cmd .= ' 1>' . $stdout_file;
    $cmd .= ' 2>' . $stderr_file;
    WIFEXITED( system $cmd) or confess "Couldn't run $cmd ($OS_ERROR)";

    return;
}

1;
