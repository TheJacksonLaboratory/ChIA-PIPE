#!/usr/bin/perl
# consider the most abound phased SNP infor at peak regions

use strict;
use warnings;
use Data::Dumper;
use File::Basename;

my $phased_snp_in = shift;
my $peak_in       = shift;
my $cluster_in    = shift;

my $phased_snp = Parsing_Phased_SNP ( $phased_snp_in );
my $peaks_anno = Peaks_Annotated_with_Phased_SNP ( $peak_in, $phased_snp );

&Clusters_Phasing_According_to_Peaks ( $cluster_in, $peaks_anno );




###########################################################################################################################################
sub Clusters_Phasing_According_to_Peaks {
    my ( $cluster_in, $peaks_anno ) = @_;

    my $clusters = Parsing_Cluster( $cluster_in );

    my $clusters_phasedType_fh = Get_Clusters_Phased_Types_and_FileHandle ( $cluster_in );

    my $i;
    print STDERR 'Phasing clusters ...', "\n";

    for my $chr ( sort { $a cmp $b } keys %{ $clusters } ){
	for my $headAnchor ( 
	    map { $_->[-1] }
	    sort { $a->[0] <=> $b->[0] || $a->[1] <=> $b->[1] }
	    map { [split(/-/, $_), $_] }
	    keys %{ $clusters->{ $chr } } ){

	for my $tailAnchor ( 
	    map { $_->[-1] }
	    sort { $a->[0] <=> $b->[0] || $a->[1] <=> $b->[1] }
	    map { [split(/-/, $_), $_] }
	    keys %{ $clusters->{ $chr }{ $headAnchor } } ){

	    $i++;
	    print STDERR "\[$i\]\r" unless $i % 1000;

	    my $PetCount = $clusters->{ $chr }{ $headAnchor }{ $tailAnchor };

	    my ( $headAnchor_s, $headAnchor_e ) = split( /-/, $headAnchor );
	    my ( $tailAnchor_s, $tailAnchor_e ) = split( /-/, $tailAnchor );

	    my $headAnchor_phased_type = Get_Anchor_Phasing_Infor ( $chr, $headAnchor_s, $headAnchor_e, $peaks_anno ); # return 'M', 'P' or 'NA'
	    my $tailAnchor_phased_type = Get_Anchor_Phasing_Infor ( $chr, $tailAnchor_s, $tailAnchor_e, $peaks_anno );

	    my $cluster_phasedType    = $headAnchor_phased_type.'_'.$tailAnchor_phased_type;
	    my $cluster_phasedGroup   = $clusters_phasedType_fh->{ $cluster_phasedType }{ 'phasedType' };
	    my $cluster_phasedType_fh = $clusters_phasedType_fh->{ $cluster_phasedType }{ 'fh'         };

	    print $cluster_phasedType_fh 
		join ( "\t", ( $chr, $headAnchor_s, $headAnchor_e, $chr, $tailAnchor_s, $tailAnchor_e, $PetCount, $cluster_phasedGroup, $cluster_phasedType ) ), "\n";

	}
	}
    }

    print STDERR "\n";
    &Close_fh ( $clusters_phasedType_fh );
    return 0;
}

sub Get_Clusters_Phased_Types_and_FileHandle {
    my $cluster_in = shift;

    my $file_name_prefix = fileparse( $cluster_in, ('.txt', '.bed') );

    my $output_clusters_M = $file_name_prefix . '.Maternal.txt';
    my $output_clusters_P = $file_name_prefix . '.Paternal.txt';
    my $output_clusters_N = $file_name_prefix . '.nonPhased.txt';
    my $output_clusters_C = $file_name_prefix . '.Crossed.txt';

    open ( my $fh_out_M, ">", $output_clusters_M ) or die $!;
    open ( my $fh_out_P, ">", $output_clusters_P ) or die $!;
    open ( my $fh_out_N, ">", $output_clusters_N ) or die $!;
    open ( my $fh_out_C, ">", $output_clusters_C ) or die $!;

    my %phaseTypes_fh = (
	'P_P'   => { 'phasedType' => 'Paternal', 'fh' => $fh_out_P },
	'P_M'   => { 'phasedType' => 'Crossed' , 'fh' => $fh_out_C },
	'P_NA'  => { 'phasedType' => 'Paternal', 'fh' => $fh_out_P },
	'M_P'   => { 'phasedType' => 'Crossed' , 'fh' => $fh_out_C },
	'M_M'   => { 'phasedType' => 'Maternal', 'fh' => $fh_out_M },
	'M_NA'  => { 'phasedType' => 'Maternal', 'fh' => $fh_out_M },
	'NA_P'  => { 'phasedType' => 'Paternal', 'fh' => $fh_out_P },
	'NA_M'  => { 'phasedType' => 'Maternal', 'fh' => $fh_out_M },
	'NA_NA' => { 'phasedType' => 'nonPhased', 'fh' => $fh_out_N },
    );

    return \%phaseTypes_fh;
}


sub Close_fh {
    my $clusters_phasedType_fh = shift;

    for my $t ( keys %{ $clusters_phasedType_fh } ){
	my $fh = $clusters_phasedType_fh->{ $t }{ 'fh' };
	close $fh;
    }
    return 0;
}


sub Get_Anchor_Phasing_Infor {
    my ( $chr, $anchor_s, $anchor_e, $peak_anno ) = @_;

#    my ( $anchor_s, $anchor_e ) = split( /-/, $anchor );

    for my $peak_s ( sort { $a <=> $b } keys %{ $peak_anno->{ $chr } } ){
	for my $peak_e ( sort { $a <=> $b } keys %{ $peak_anno->{ $chr }{ $peak_s } } ){

	    my $peak_phased_type = $peak_anno->{ $chr }{ $peak_s }{ $peak_e };
#case I
	    if ( $anchor_e >= $peak_s && $anchor_e <= $peak_e  ){
		return $peak_phased_type;
	    }
#case II
	    elsif ( $anchor_s >= $peak_s && $anchor_s <= $peak_e ){
		return $peak_phased_type;
	    }
#case III
	    elsif ( $anchor_s <= $peak_s && $anchor_e >= $peak_e ){
		return $peak_phased_type;
	    } 
	    else {
		next;
	    }
	}
    }
    return 'NA';
}



sub Peaks_Annotated_with_Phased_SNP {
    my $peak_in = shift;
    my $snp     = shift;

    my $peaks = Parsing_Peak( $peak_in );

    my $output_prefix = fileparse( $peak_in, ( '.txt', '.bed' ) );

    my $output_phased_peak_anno     = $output_prefix . '.SNP_Phased.annotation.txt';
    my $output_phased_peak_browser  = $output_prefix . '.SNP_Phased.browser.bed';
    my $output_phased_peak_maternal = $output_prefix . '.SNP_Phased.Maternal.txt';
    my $output_phased_peak_paternal = $output_prefix . '.SNP_Phased.Paternal.txt';

    open ( my $fh_out_anno, ">", $output_phased_peak_anno     ) or die $!;
    open ( my $fh_out_bed,  ">", $output_phased_peak_browser  ) or die $!;
    open ( my $fh_out_M,    ">", $output_phased_peak_maternal ) or die $!;
    open ( my $fh_out_P,    ">", $output_phased_peak_paternal ) or die $!;
    
    my %peak_anno;
    my $i;

    print STDERR 'Peaks annotation with phased SNPs ...', "\n";

    for my $chr ( sort { $a cmp $b } keys %{ $peaks } ){
	for my $peak_s ( sort { $a <=> $b } keys %{ $peaks->{ $chr } } ){
	    my $peak_e = $peaks->{ $chr }{ $peak_s };

	    $i++;
	    print STDERR "[$i]\r" unless $i % 100;
	
	    my ($total_phased_snp, $phased_type, $peak_anno) = Overlapping_with_Phased_SNP ( $chr, $peak_s, $peak_e, $phased_snp );
	
	    if ( $total_phased_snp ){

#print phased peak annotation
		print $fh_out_anno $chr.':'.$peak_s.'-'.$peak_e, "\t", $total_phased_snp, "\t", $phased_type, "\t", join( "\|", @{ $peak_anno } ), "\n";

#print Maternal peak
		if ( $phased_type eq 'M' ){

		    $peak_anno{ $chr }{ $peak_s }{ $peak_e } = 'M';

		    print $fh_out_M join( "\t", ( $chr, $peak_s, $peak_e ) ), "\n";
		    print $fh_out_bed join( "\t", ( $chr, $peak_s, $peak_e, 'Maternal', '215,48,39' ) ), "\n";
#print Maternal peak
		} elsif ( $phased_type eq 'P' ){

		    $peak_anno{ $chr }{ $peak_s }{ $peak_e } = 'P';

		    print $fh_out_P join( "\t", ( $chr, $peak_s, $peak_e ) ), "\n";
		    print $fh_out_bed join( "\t", ( $chr, $peak_s, $peak_e, 'Paternal', '26,152,80' ) ), "\n";
#print to both
		} elsif ( $phased_type eq 'Mix' ){
		    print $fh_out_M join( "\t", ( $chr, $peak_s, $peak_e ) ), "\n";
		    print $fh_out_P join( "\t", ( $chr, $peak_s, $peak_e ) ), "\n";
		    print $fh_out_bed join( "\t", ( $chr, $peak_s, $peak_e, 'Mix', '209,229,240' ) ), "\n";
#unknown case
		} else {
		    print STDERR 'Unknow case for peak phasing annotation!', "\n";
		    exit 1;
		}

#non phasing infor
#print to both
	    } else {
		print $fh_out_M join( "\t", ( $chr, $peak_s, $peak_e ) ), "\n";
		print $fh_out_P join( "\t", ( $chr, $peak_s, $peak_e ) ), "\n";
		print $fh_out_bed join( "\t", ( $chr, $peak_s, $peak_e, 'Unphased', '153,153,153' ) ), "\n";
	    }
	}
    }

    print STDERR "\n";

    close $fh_out_anno;
    close $fh_out_bed;
    close $fh_out_M;
    close $fh_out_P;

    return \%peak_anno;
}


sub Overlapping_with_Phased_SNP {
    my ( $chr, $peak_s, $peak_e, $phased_snp ) = @_;

    my ( @peak_anno, %anno_type, $total_phased_snp );
#    my ( @peak_anno_unphased, $total_unphased_snp );

    my $max_snp_cov = 0;
    my $max_snp_cov_type;

    for my $snp_pos ( sort { $a <=> $b } keys %{ $phased_snp->{ $chr } } ){

	my $snp_infor = $phased_snp->{ $chr }{ $snp_pos };
	my $snp_type  = $snp_infor->[0];

	my ( $snp_freq_M, $snp_freq_P ) = split( /:/, $snp_infor->[1] );


	if ( $snp_pos >= $peak_s && $snp_pos <= $peak_e && $snp_type ne 'N' ){

	    my $snp_infor_line = join( ",", @{ $snp_infor } );
	    push @peak_anno, $snp_pos .','. $snp_infor_line;
	    $anno_type{ $snp_type }++;
	    $total_phased_snp++;

# record max cov snp infor
	    if ( ( $snp_freq_M + $snp_freq_P ) > $max_snp_cov ){
		$max_snp_cov = $snp_freq_M + $snp_freq_P;
		$max_snp_cov_type = $snp_infor->[0];
	    }

	}
	elsif ( $snp_pos >= $peak_s && $snp_pos <= $peak_e && $snp_type eq 'N' ){

# record max cov snp infor
	    if ( ( $snp_freq_M + $snp_freq_P ) > $max_snp_cov ){
		$max_snp_cov = $snp_freq_M + $snp_freq_P;
		$max_snp_cov_type = $snp_infor->[0];
	    }
	}
    }


# return phased type
    if ( ! $max_snp_cov_type ){
	return ( 0, 'Unphased', 'NA' );
    }
    elsif ( $max_snp_cov_type eq 'N' ){
	return ( 0, 'Unphased', 'NA' );
    } else {

	if ( $total_phased_snp ){
	    if ( $anno_type{ 'M' } && $anno_type{ 'P' } ){
		my $phased_type = 'Mix';
		return ( $total_phased_snp, $phased_type, \@peak_anno );
	    } elsif ( $anno_type{ 'M' } && ! $anno_type{ 'P' } ){
		my $phased_type = 'M';
		return ( $total_phased_snp, $phased_type, \@peak_anno );
	    } elsif ( ! $anno_type{ 'M' } && $anno_type{ 'P' } ){
		my $phased_type = 'P';
		return ( $total_phased_snp, $phased_type, \@peak_anno );
	    } else {
		print STDERR 'Unknow case !';
		exit 1;
	    }
	} else {
	    return ( 0, 'Unphased', 'NA' );
	}
    }
}


sub Parsing_Phased_SNP {
    my $file = shift;

    my %phased_snp;

    open ( my $fh_in, $file ) or die $!;

    while ( my $line = <$fh_in> ){
	next if ( $line =~ /^chrom/ );
	chomp $line;
	my @lines = split(/\t/, $line);

	my $snp_freq = $lines[3].':'.$lines[4];

	my $snp = $lines[2];
	$snp =~ s/\|/\:/;

	if ( $lines[-1] eq 'No' ){

	    $phased_snp{ $lines[0] }{ $lines[1] } = ['N', $snp_freq, $snp];

	} elsif ( $lines[-1] eq 'Yes' && ( $lines[3] > $lines[4] ) ){

	    $phased_snp{ $lines[0] }{ $lines[1] } = ['P', $snp_freq, $snp];

	} elsif ( $lines[-1] eq 'Yes' && ( $lines[3] < $lines[4] ) ){

	    $phased_snp{ $lines[0] }{ $lines[1] } = ['M', $snp_freq, $snp];

	} else {
	    print STDERR 'Unknow case: ', $line, "\n";
	}
    }

    close $fh_in;
    return \%phased_snp;
}


sub Parsing_Peak {
    my $file = shift;

    my %peaks;

    open ( my $fh_in, $file ) or die $!;

    while ( my $line = <$fh_in> ){
	chomp $line;
	my @lines = split(/\t/, $line);
	$peaks{ $lines[0] }{ $lines[1] } = $lines[2];
    }
    close $fh_in;
    return \%peaks;
}


sub Parsing_Cluster {
    my $file = shift;
    my %clusters;

    open ( my $fh_in, $file ) or die $!;

    while ( my $line = <$fh_in> ){
	chomp $line;
	my @lines = split( /\t/, $line );
	$clusters{ $lines[0] }{ $lines[1].'-'.$lines[2] }{ $lines[4].'-'.$lines[5] } = $lines[6];
    }
    close $fh_in;
    return \%clusters;
}
