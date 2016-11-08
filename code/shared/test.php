


<?php
  
  

$dm3_exons=file("C:/xampp/htdocs/labwork/data/shared/dm3_exons_1.txt");
$new_exons=file("C:/xampp/htdocs/labwork/data/shared/new_exons.txt");
$chrom_new=file("C:/xampp/htdocs/labwork/data/shared/chromInfo_new.txt");
  
  
  
  
 
  foreach ($dm3_exons as $exons_line)
  {
    
	//echo htmlspecialchars($exons_line) . '<br>';
	$split_string_exons = explode('	', $exons_line);
	//echo $split_string_exons[1]."<br>";
	
 
	
	foreach ($chrom_new as $chrom_line)
		{
			$split_string_chrom = explode('	', $chrom_line);


		//echo htmlspecialchars($chrom_line) . '<br>';
			if ( $split_string_chrom[0] == $split_string_exons[0])
				{
					//echo $split_string_chrom[0]."==". $split_string_exons[0]."<br>";

					 $exons_start = (int)$split_string_exons[1];
					 $exons_end =(int)$split_string_exons[2];
					 $exons_length = $exons_end - $exons_start;
					 
					 $chrom_start=(int) $split_string_chrom[1];
					 $chrom_end= (int)$split_string_chrom[2];
					
					 $new_exons_start = $chrom_start+$exons_start;
					 $new_exons_end = $new_exons_start+$exons_length;
					
					
					
					 $new_exons_line= $split_string_chrom[0]."\t".$new_exons_start."\t".$new_exons_end;
					
				
				}
		 
		}
 
  }
 
?>

