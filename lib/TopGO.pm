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
  read_gene_info
  read_go_terms
  write_gene_list
  write_mapping_file
  annotate_with_genes
);

=head1 SYNOPSIS

    # Brief code examples

=cut

=func read_gene_info

  Usage       : ($p_value_for, $name_for, $description_for) = read_gene_info(
                    $input_file, $has_header, $gene_field, $p_value_field,
                    $name_field, $description_field
                );
  Purpose     : Read p values, name and description from a file and organise by
                gene ID
  Returns     : Hashref { String (gene ID) => Float (p value) }
                Hashref { String (gene ID) => String (name) }
                Hashref { String (gene ID) => String (description) }
  Parameters  : String (the input filename)
                Boolean (whether input file has header)
                Int (the gene field)
                Int (the p value field)
                Int (the name field) or undef
                Int (the description field) or undef
  Throws      : No exceptions
  Comments    : None

=cut

sub read_gene_info {
    my (
        $input_file,    $has_header, $gene_field,
        $p_value_field, $name_field, $description_field
    ) = @_;

    my ( %p_value_for, %name_for, %description_for );

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
        my $name = $name_field ? $fields[ $name_field - 1 ] : q{};
        my $description =
          $description_field ? $fields[ $description_field - 1 ] : q{};
        next if $p_value eq q{-} || $p_value eq 'NA';    # Ignore filtered genes

        foreach my $gene_id ( split /,/xms, $gene_ids ) {
            if ( exists $p_value_for{$gene_id}
                && $p_value_for{$gene_id} < $p_value )
            {
                # Keep lowest p value if gene already seen
                $p_value = $p_value_for{$gene_id};
            }
            $p_value_for{$gene_id}     = $p_value;
            $name_for{$gene_id}        = $name;
            $description_for{$gene_id} = $description;
        }
    }
    close $fh;

    return \%p_value_for, \%name_for, \%description_for;
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

    my $stdout_file = $arg_ref->{output_prefix} . '.o';
    my $stderr_file = $arg_ref->{output_prefix} . '.e';

    my $cmd = join q{ }, $arg_ref->{r_binary}, '--slave', '--args',
      $arg_ref->{gene_list_file}, $arg_ref->{mapping_file}, $arg_ref->{domain},
      $arg_ref->{output_prefix},  $arg_ref->{sig_level},    '<',
      $arg_ref->{topgo_script};
    $cmd .= ' 1>' . $stdout_file;
    $cmd .= ' 2>' . $stderr_file;
    WIFEXITED( system $cmd) or confess "Couldn't run $cmd ($OS_ERROR)";

    return;
}

=func annotate_with_genes

  Usage       : TopGO::annotate_with_genes( {
                    input_file   => 'all.tsv',
                    output_file  => 'all.genes.tsv',
                    p_values     => $p_value_for,
                    names        => $name_for,
                    descriptions => $descriptions_for,
                } );
  Purpose     : Add gene annotation to topGO results
  Returns     : undef
  Parameters  : Hashref {
                    input_file   => String (the input file),
                    output_file  => String (the output file),
                    p_values     => Hashref (of p values keyed by gene ID),
                    names        => Hashref (of names keyed by gene ID),
                    descriptions => Hashref (of descriptions keyed by gene ID),
                }
  Throws      : If input file is missing
                If output file is missing
                If p values are missing
                If names are missing
                If descriptions are missing
  Comments    : None

=cut

sub annotate_with_genes {
    my ($arg_ref) = @_;

    confess 'No input file specified'   if !defined $arg_ref->{input_file};
    confess 'No output file specified'  if !defined $arg_ref->{output_file};
    confess 'No p values specified'     if !defined $arg_ref->{p_values};
    confess 'No names specified'        if !defined $arg_ref->{names};
    confess 'No descriptions specified' if !defined $arg_ref->{descriptions};

    my $p_value_for     = $arg_ref->{p_values};
    my $name_for        = $arg_ref->{names};
    my $description_for = $arg_ref->{descriptions};

    open my $fh_in,  '<', $arg_ref->{input_file};
    open my $fh_out, '>', $arg_ref->{output_file};

    # Get input header
    my $header = <$fh_in>;
    my @header_fields = split /\t/xms, $header;
    pop @header_fields;

    # Write output header
    push @header_fields, 'Gene', 'p value', 'Name', 'Description';
    print {$fh_out} ( join "\t", @header_fields ), "\n";

    # Rewrite input so GO terms are repeated for each gene they are annotatd to
    while ( my $line = <$fh_in> ) {
        chomp $line;
        my @fields = split /\t/xms, $line;
        my @genes  = split /,/xms,  pop @fields;
        foreach
          my $gene ( sort { $p_value_for->{$a} <=> $p_value_for->{$b} } @genes )
        {
            print {$fh_out} (
                join "\t", @fields, $gene, $p_value_for->{$gene},
                $name_for->{$gene}, $description_for->{$gene}
              ),
              "\n";
        }
    }

    close $fh_in;
    close $fh_out;

    return;
}

=func filter_by_significance

  Usage       : filter_by_significance($input_file, $output_file, $sig_level);
  Purpose     : Filter file by significance level
  Returns     : undef
  Parameters  : String (the input filename)
                String (the output filename)
                Float (the significance level)
  Throws      : No exceptions
  Comments    : None

=cut

sub filter_by_significance {
    my ( $input_file, $output_file, $sig_level ) = @_;

    open my $fh_in,  '<', $input_file;
    open my $fh_out, '>', $output_file;

    # Write header
    my $header = <$fh_in>;
    print {$fh_out} $header;

    # Write significant lines
    while ( my $line = <$fh_in> ) {
        my @fields = split /\t/xms, $line;
        next if $fields[5] >= $sig_level;
        print {$fh_out} $line;
    }

    close $fh_in;
    close $fh_out;

    return;
}

1;
