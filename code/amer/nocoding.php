<?php



  

$coding=file("C:/xampp/htdocs/labwork/data/amer/coding.txt");
$nocoding= fopen("C:/xampp/htdocs/labwork/data/amer/nocoding.txt", 'a') or exit("Unable to open file!");
  $chrom_new=file("C:/xampp/htdocs/labwork/data/amer/chromInfo_new.txt");

  
  

 
 $j=0;
 //non coding start
  $ncs=0;
foreach ($chrom_new as $chrom_line)
	{
		//echo"j :".$j."<br>";
		$split_string_chrom = explode('	', $chrom_line);

		echo $split_string_chrom[1] ."---".$split_string_chrom[2]."<br>";
		$i=0;
		foreach ($coding as $coding_line)
		{
			//$split_coding = explode('	', $coding_line);
			
			
			$split_coding_i = explode('	', $coding[$i]);
			$split_coding_next_i = explode('	', $coding[$i+1]);
		if (($ncs==0) &&($split_string_chrom[0] == $split_coding_i[0]) )
		{
			
			$coding_i_end =(int)$split_coding_i[1];
			$non_coding_end=$coding_i_end-1;
			$non_coding_start= 0;
			$non_coding_lengh=$non_coding_end-$non_coding_start;
			//echo $non_coding_start."##".$non_coding_end."##". $non_coding_lengh."<br>";					
			//write to txt non coding file ============================
			$new_nocoding_ine= $non_coding_start."\t".$non_coding_end."\t". $non_coding_lengh."\n";
			fwrite($nocoding, $new_nocoding_ine);
			$ncs=1;
			
		}
		
			if ( ($split_string_chrom[0] == $split_coding_i[0]) &&($split_coding_i[0]==$split_coding_next_i [0] ))
				{
					

					 $coding_i_start = (int)$split_coding_i[1];
					 $coding_i_end =(int)$split_coding_i[2];
					 
					 $coding_next_i_start = (int)$split_coding_next_i[1];
					 $coding_next_i_end =(int)$split_coding_next_i[2];
					

					$non_coding_start= $coding_i_end+1;
					$non_coding_end= $coding_next_i_start-1;
					$non_coding_lengh=$non_coding_end-$non_coding_start;
					//echo $non_coding_start."==".$non_coding_end."==". $non_coding_lengh."<br>";					
					$new_nocoding_ine= $non_coding_start."\t".$non_coding_end."\t". $non_coding_lengh."\n";
					fwrite($nocoding, $new_nocoding_ine);
				}
				
			// check if the next row is for new exon or not 	
			if ( ($split_string_chrom[0] != $split_coding_next_i[0]) && ($split_string_chrom[0]== $split_coding_i[0]))
			{
					 $coding_i_start = (int)$split_coding_i[1];
					 $coding_i_end =(int)$split_coding_i[2];
					 
					 $coding_next_i_start = (int)$split_coding_next_i[1];
					 $coding_next_i_end =(int)$split_coding_next_i[2];
					

					$non_coding_start= $coding_i_end+1;
					$non_coding_end= $split_string_chrom[2];
					$non_coding_lengh=$non_coding_end-$non_coding_start;
					//echo $non_coding_start."==".$non_coding_end."==". $non_coding_lengh."<br>";	
					$new_nocoding_ine= $non_coding_start."\t".$non_coding_end."\t". $non_coding_lengh."\n";
					fwrite($nocoding, $new_nocoding_ine);
	
			}
$i++; 
				 
		}
	  	  
		
$ncs=0;	

$j++;
  }
  
 
?>



  
