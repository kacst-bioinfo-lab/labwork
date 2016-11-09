<?php

$file=fopen("C:/xampp/htdocs/labwork/data/amer/chromInfo.txt","r") or exit("Unable to open file!");
$handle = fopen("C:/xampp/htdocs/labwork/data/amer/chromInfo_new.txt", 'a') or exit("Unable to open file!");



while (!feof($file))
  {
 //echo fgetc($file). "<br>";  read a single character from a file.
  
  $line =fgets($file);
  $split_string = explode('	', $line);
  //echo $line. "<br>";
  $new_line= $split_string[0]. "\t" .  "0" . "\t" . $split_string[1] . "\t" .$split_string[2];
 //echo $new_line. "<br><br>";
 fwrite($handle, $new_line);
 

  
  }
fclose($file);
fclose($handle);
echo "New chromInfo file was successfully created <br>";


  

$dm3_exons=file("C:/xampp/htdocs/labwork/data/amer/dm3_exons.txt");
//$new_exons=file("C:/xampp/htdocs/labwork/data/amer/dm3_exons_new.txt");
$chrom_new=file("C:/xampp/htdocs/labwork/data/amer/chromInfo_new.txt");
$coding= fopen("C:/xampp/htdocs/labwork/data/amer/coding.txt", 'a') or exit("Unable to open file!");
  
  
  
 
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
					
					
					
					 $new_exons_line= $split_string_chrom[0]."\t".$new_exons_start."\t".$new_exons_end."\t".$exons_length."\n";
					
				    fwrite($coding, $new_exons_line);
				}
		 
		}
  

  }
 
 echo "New exons file was successfully created <br>";
 fclose($coding);
?>



  
