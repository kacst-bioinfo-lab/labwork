#!/usr/bin/perl
use strict;
use warnings;
no warnings 'uninitialized';
#use Archive::Extract;
use LWP::Simple;
use File::Copy;


unlink('../../data/amer/coding.bed');
unlink('../../data/amer/coding_CEs.bed');
unlink('../../data/amer/coding_spacer.bed');
unlink('../../data/amer/intergenic.bed');
unlink('../../data/amer/intergenic_CNEs.bed');
unlink('../../data/amer/intergenic_spacer.bed');
unlink('../../data/amer/intronic.bed');
unlink('../../data/amer/intronic_CNEs.bed');
unlink('../../data/amer/intronic_spacer.bed');
unlink('../../data/amer/noncoding.bed');
#===============================================================================================
# download chromInfo file
#my $chromInfo_url = 'http://hgdownload.soe.ucsc.edu/goldenPath/dm3/database/chromInfo.txt.gz';
#my $chromInfo_file = 'chromInfo.txt.gz';
#getstore($chromInfo_url, $chromInfo_file);
#my $extr_chrom = Archive::Extract->new( archive => 'chromInfo.txt.gz' );
#my $ok = $extr_chrom->extract;
#system("cut -f1-2 chromInfo.txt > allgenomefile");
#===============================================================================================

#===============================================================================================
# download most conserved elements file
#my $conserved_elements_url = 'http://hgdownload.soe.ucsc.edu/goldenPath/dm3/database/phastConsElements15way.txt.gz';
#my $conserved_elements_file = 'phastConsElements15way.txt.gz';
#getstore($conserved_elements_url, $conserved_elements_file);
#my $extr_conserved_elements = Archive::Extract->new( archive => 'phastConsElements15way.txt.gz' );
#my $ok = $extr_conserved_elements->extract;
#system("cut -f1-2 phastConsElements15way.txt > cefile.txt");
#===============================================================================================






#===============================================================================================
# sort the allgenomefile.txt into sorted_genome.txt
system("sort -k 1,1 -k2,2n ../../data/shared/chromInfo.txt > ../../data/amer/sorted_genome.txt");


#===============================================================================================
# sort the dm3_exons.txt to sorted_dm3_exons.bed
system("sort -k 1,1 -k2,2n ../../data/shared/dm3_exons.txt > ../../data/amer/sorted_dm3_exons.bed");


#===============================================================================================
# creat new genome file with valid entries  by perl
#open my $valid_g , '>', '../../data/amer/valid_genome.txt';

my $valid_genome = '../../data/amer/valid_genome.txt';
open(my $vg, '>>', $valid_genome) or die "Could not open file '$valid_genome' $!";




my $sorted_genome_file = '../../data/amer/sorted_genome.txt';
open my $sorted_genome, $sorted_genome_file or die "Could not open $sorted_genome_file: $!";
while( my $sorted_genome_line = <$sorted_genome>){   
	my @sorted_genome_line_split= split (/\t/ ,  $sorted_genome_line);
	 
	print $vg "$sorted_genome_line_split[0]\t$sorted_genome_line_split[1]\n" 
}

close $sorted_genome;
close $vg;


#==============================================================================================
# creat the noncoding file 
system("bedtools complement -i ../../data/amer/sorted_dm3_exons.bed -g ../../data/amer/valid_genome.txt > ../../data/amer/noncoding.bed");



#==============================================================================================
#  added 0 as start last genom file

open my $f_genome , '>', '../../data/amer/final_genome.bed';

my $final_genome = '../../data/amer/final_genome.bed';
open(my $fg, '>>', $final_genome) or die "Could not open file '$final_genome' $!";

 

my $valid_genome_file = '../../data/amer/valid_genome.txt';
open my $v_genome, $valid_genome_file or die "Could not open $valid_genome_file: $!";
while( my $v_genome_line = <$v_genome>){   
	my @v_genome_line_split= split (/\t/ ,  $v_genome_line);
	 

	print $fg "$v_genome_line_split[0]\t0\t$v_genome_line_split[1]"
}

close $fg;
close $v_genome;



#==============================================================================================
# creat the coding file 
#copy("../../data/amer/sorted_dm3_exons.bed" , "../../data/amer/coding.bed") ;

# or by final genom after after added to valid genom file 0 

system("bedtools subtract -a ../../data/amer/final_genome.bed -b ../../data/amer/noncoding.bed > ../../data/amer/coding.bed");




#==============================================================================================
open(my $fh, '>', '../../data/amer/intronic.bed') or die "Could not open file ''../../data/amer/intronic.bed'' $!";


open my $info1 , '<', '../../data/amer/coding.bed';
my @arry_line=<$info1>;
#print @arry_line;

open my $info, '<', '../../data/amer/coding.bed';
while( my $line = <$info>){   
  	#print $line; 
	##my $nextline=<$info>; 
	##print $nextline; 
	my $line_number=$.;
	my @current_split= split (/\t/ ,  $line);
	my @next_split= split (/\t/ ,  $arry_line[$line_number]);

	if ($current_split[0] eq $next_split[0])
	{
	my $start = $current_split[2]+1-1;
	my $end = $next_split[1];
	my $length = $end -  $start ; 
	print $fh "$current_split[0]\t$start\t$end\n";
	}
 

	#if ($current_split[0] ne $next_split[0])
	#{

	#}

}

close $info;
close $info1;
close $fh;

system("bedtools subtract -a ../../data/amer/noncoding.bed -b ../../data/amer/intronic.bed > ../../data/amer/intergenic.bed"); 




#========== calculate coding  length ===========
my $coding_file = '../../data/amer/coding.bed';
my $sum_coding=0;
open my $coding_length, $coding_file or die "Could not open $coding_file: $!";
while( my $coding_length_line = <$coding_length>){   
	my @coding_length_split= split (/\t/ ,  $coding_length_line);
	 
	$sum_coding= $sum_coding+($coding_length_split[2]-$coding_length_split[1]);
}

close $coding_length;
print "Total lenght of coding                         : $sum_coding\n";

#==================================================





#========== calculate noncoding  length ===========

 
my $noncoding_file = '../../data/amer/noncoding.bed';

my $sum_noncoding=0;
open my $noncoding_length, $noncoding_file or die "Could not open $noncoding_file: $!";
while( my $noncoding_length_line = <$noncoding_length>){   
	my @noncoding_length_split= split (/\t/ ,  $noncoding_length_line);
	 
	$sum_noncoding= $sum_noncoding+($noncoding_length_split[2]-$noncoding_length_split[1]);
}

close $noncoding_length;
print "Total lenght of noncoding                      : $sum_noncoding\n";
my $sum_coding_noncoding=$sum_noncoding+$sum_coding;
print "Total of lenght of noncoding and coding        : $sum_coding_noncoding \n";
#==================================================




#========== calculate genome  length ===========

 
my $genome_file = '../../data/amer/valid_genome.txt';

my $sum_genome=0;
open my $genome_length, $genome_file or die "Could not open $genome_file: $!";
while( my $genome_length_line = <$genome_length>){   
	my @genome_length_split= split (/\t/ ,  $genome_length_line);
	 
	$sum_genome= $sum_genome+$genome_length_split[1];
}

close $genome_length;
print "Total lenght of genome                         : $sum_genome\n \n\n";
#==================================================



#========== calculate intronic length ===========
my $intronic_file = '../../data/amer/intronic.bed';

my $sum_intronic=0;
open my $intronic_length, $intronic_file or die "Could not open $intronic_file: $!";
while( my $intronic_length_line = <$intronic_length>){   
	my @intronic_length_split= split (/\t/ ,  $intronic_length_line);
	 
	$sum_intronic= $sum_intronic+($intronic_length_split[2]-$intronic_length_split[1]);
}

close $intronic_length;
print "Total lenght of intronic                       : $sum_intronic\n";
#==================================================



#========== calculate intergenic length ===========
my $intergenic_file = '../../data/amer/intergenic.bed';

my $sum_intergenic=0;
open my $intergenic_length, $intergenic_file or die "Could not open $intergenic_file: $!";
while( my $intergenic_length_line = <$intergenic_length>){   
	my @intergenic_length_split= split (/\t/ ,  $intergenic_length_line);
	 
	$sum_intergenic= $sum_intergenic+($intergenic_length_split[2]-$intergenic_length_split[1]);
}

close $intergenic_length;
print "Total lenght of intergenic                     : $sum_intergenic\n";

my $sum_intergenic_intronic = $sum_intergenic+ $sum_intronic;
print "Total lenght of intergenic and intronic        : $sum_intergenic_intronic\n";
print "Total lenght of noncoding                      : $sum_noncoding\n";
 

 




#==============================================================================================
#  remove first column from  most conserved elements file

open my $f_c_e , '>', '../../data/amer/cnserved_elements.bed';

my $new_cnserved_elements = '../../data/amer/cnserved_elements.bed';
open(my $fce, '>>', $new_cnserved_elements) or die "Could not open file '$new_cnserved_elements' $!";



my $phastConsElements15way = '../../data/shared/phastConsElements15way.txt';
open my $orginal_c_e, $phastConsElements15way or die "Could not open $phastConsElements15way: $!";
while( my $orginal_c_e_line = <$orginal_c_e>){   
	my @orginal_c_e_line_split= split (/\t/ ,  $orginal_c_e_line);
	 

	print $fce "$orginal_c_e_line_split[1]\t$orginal_c_e_line_split[2]\t$orginal_c_e_line_split[3]\t$orginal_c_e_line_split[4]\t$orginal_c_e_line_split[5]"
}

close $fce;
close $orginal_c_e;



#==============================================================================================
# creat CEs coding file 
system("bedtools intersect -a ../../data/amer/coding.bed -b ../../data/amer/cnserved_elements.bed > ../../data/amer/coding_CEs.bed");


# creat Spacer coding file 
system("bedtools subtract -a ../../data/amer/coding.bed -b ../../data/amer/coding_CEs.bed > ../../data/amer/coding_spacer.bed");


#==============================================================================================
# creat CNEs intergenic file 
system("bedtools intersect -a ../../data/amer/intergenic.bed -b ../../data/amer/cnserved_elements.bed > ../../data/amer/intergenic_CNEs.bed");


# creat Spacer intergenic file 
system("bedtools subtract -a ../../data/amer/intergenic.bed -b ../../data/amer/intergenic_CNEs.bed > ../../data/amer/intergenic_spacer.bed");



#==============================================================================================
# creat CNEs intronic file 
system("bedtools intersect -a ../../data/amer/intronic.bed -b ../../data/amer/cnserved_elements.bed > ../../data/amer/intronic_CNEs.bed");


# creat Spacer intronic file 
system("bedtools subtract -a ../../data/amer/intronic.bed -b ../../data/amer/intronic_CNEs.bed > ../../data/amer/intronic_spacer.bed");







#========== calculate coding CNEs length ======================================================
my $coding_CEs_file = '../../data/amer/coding_CEs.bed';

my $sum_coding_CEs=0;
open my $coding_CEs_length, $coding_CEs_file or die "Could not open $coding_CEs_file: $!";
while( my $coding_CEs_length_line = <$coding_CEs_length>){   
	my @coding_CEs_split= split (/\t/ ,  $coding_CEs_length_line);
	 
	$sum_coding_CEs= $sum_coding_CEs+($coding_CEs_split[2]-$coding_CEs_split[1]);
}

close $coding_CEs_length;

print"\n*************************************************************************************\n";
print "Total lenght of coding CEs                     : $sum_coding_CEs\n";
#=============================================================================================
#========== calculate coding spacer length ===================================================
my $coding_spacer_file = '../../data/amer/coding_spacer.bed';

my $sum_coding_spacer=0;
open my $coding_spacer_length, $coding_spacer_file or die "Could not open $coding_spacer_file: $!";
while( my $coding_spacer_length_line = <$coding_spacer_length>){   
	my @coding_spacer_split= split (/\t/ ,  $coding_spacer_length_line);
	 
	$sum_coding_spacer= $sum_coding_spacer+($coding_spacer_split[2]-$coding_spacer_split[1]);
}

close $coding_spacer_length;
print "Total lenght of coding spacer                  : $sum_coding_spacer\n";
#==============================================================================================

my $sum_coding_CNEs_and_spacer= $sum_coding_CEs+$sum_coding_spacer ;
print "Total lenght of CNEs and spacer coding         : $sum_coding_CNEs_and_spacer\n";
print "Total lenght of coding                         : $sum_coding\n";
print"\n*************************************************************************************\n";






#========== calculate intergenic CNEs length ===========
my $intergenic_CNEs_file = '../../data/amer/intergenic_CNEs.bed';

my $sum_intergenic_CNEs=0;
open my $intergenic_CNEs_length, $intergenic_CNEs_file or die "Could not open $intergenic_CNEs_file: $!";
while( my $intergenic_CNEs_length_line = <$intergenic_CNEs_length>){   
	my @intergenic_CNEs_split= split (/\t/ ,  $intergenic_CNEs_length_line);
	 
	$sum_intergenic_CNEs= $sum_intergenic_CNEs+($intergenic_CNEs_split[2]-$intergenic_CNEs_split[1]);
}

close $intergenic_CNEs_length;
print "Total lenght of intergenic CNEs                : $sum_intergenic_CNEs\n";
#=================================================================================================
#========== calculate intergenic spacer length ===================================================
my $intergenic_spacer_file = '../../data/amer/intergenic_spacer.bed';

my $sum_intergenic_spacer=0;
open my $intergenic_spacer_length, $intergenic_spacer_file or die "Could not open $intergenic_spacer_file: $!";
while( my $intergenic_spacer_length_line = <$intergenic_spacer_length>){   
	my @intergenic_spacer_split= split (/\t/ ,  $intergenic_spacer_length_line);
	 
	$sum_intergenic_spacer= $sum_intergenic_spacer+($intergenic_spacer_split[2]-$intergenic_spacer_split[1]);
}

close $intergenic_spacer_length;
print "Total lenght of intergenic spacer              : $sum_intergenic_spacer\n";
#=================================================================================================



my $sum_intergenic_CNEs_and_spacer= $sum_intergenic_CNEs+$sum_intergenic_spacer ;
print "Total lenght of CNEs and spacer coding         : $sum_intergenic_CNEs_and_spacer\n";
print "Total lenght of intergenic                     : $sum_intergenic\n";
print"\n*************************************************************************************\n";






#========== calculate intronic CNEs length =========================================================
my $intronic_CNEs_file = '../../data/amer/intronic_CNEs.bed';

my $sum_intronic_CNEs=0;
open my $intronic_CNEs_length, $intronic_CNEs_file or die "Could not open $intronic_CNEs_file: $!";
while( my $intronic_CNEs_length_line = <$intronic_CNEs_length>){   
	my @intronic_CNEs_split= split (/\t/ ,  $intronic_CNEs_length_line);
	 
	$sum_intronic_CNEs= $sum_intronic_CNEs+($intronic_CNEs_split[2]-$intronic_CNEs_split[1]);
}

close $intronic_CNEs_length;
print "Total lenght of intronic CNEs                  : $sum_intronic_CNEs\n";
#=====================================================================================================
#========== calculate intronic spacer length =========================================================
my $intronic_spacer_file = '../../data/amer/intronic_spacer.bed';

my $sum_intronic_spacer=0;
open my $intronic_spacer_length, $intronic_spacer_file or die "Could not open $intronic_spacer_file: $!";
while( my $intronic_spacer_length_line = <$intronic_spacer_length>){   
	my @intronic_spacer_split= split (/\t/ ,  $intronic_spacer_length_line);
	 
	$sum_intronic_spacer= $sum_intronic_spacer+($intronic_spacer_split[2]-$intronic_spacer_split[1]);
}

close $intronic_spacer_length;
print "Total lenght of intronic spacer                : $sum_intronic_spacer\n";
#=====================================================================================================



my $sum_intronic_CNEs_and_spacer= $sum_intronic_CNEs+$sum_intronic_spacer ;
print "Total lenght of CNEs and spacer coding         : $sum_intronic_CNEs_and_spacer\n";
print "Total lenght of intronic                       : $sum_intronic\n";
print"\n*************************************************************************************\n";



unlink('../../data/amer/cnserved_elements.bed') or die "Could not delete the file!\n";
unlink('../../data/amer/final_genome.bed') or die "Could not delete the file!\n";
unlink('../../data/amer/sorted_dm3_exons.bed') or die "Could not delete the file!\n";
unlink('../../data/amer/sorted_genome.txt') or die "Could not delete the file!\n";
unlink('../../data/amer/valid_genome.txt') or die "Could not delete the file!\n";


