#!/usr/bin/php
<?php
/* 1. Consider joins within each scaffold only.
 * 2. Let P (denominator of sensitivity) = number of joins in answer.
 * 3. Let Y (denominator of precision) = TP + FP = number of joins in result.
 * 4. Coverage = 0.5 * (sequence length of correct-joined contigs) / total length of all contigs in answer
 *    (half lengths of the first and last contigs in a scaffold are counted in order to get 100% coverage in the result).
 */

$mesg = "Usage: php misJoin_intersect_v7.php <answer file> <.car file> [-co <contigMerged file>, -fw <output file>, -a (more info)]\n";
if($argc < 2){
	die($mesg);
}

$detail = 0;
$cov = 0;
$fw = 0;
$nga = 0;

foreach($argv as $k => $each){
	if($each == "-h"){
		die($mesg);
	}
	else if($each == "-a"){
		$detail = 1;
	}
	else if($each == "-co"){
		if(!is_file($argv[$k+1])){
			die($mesg);
		}
		$contigMer = $argv[$k+1];
		$cov = 1;
		
	}
	else if($each == "-fw"){
		$fw = 1;
		$outputFile = $argv[$k+1];
	}
	else if($each == "-nga"){
		$nga = 1;
		$nga_Dir = $argv[$k+1];
	}
}

$answer = file($argv[1], FILE_IGNORE_NEW_LINES);
$merged = file($argv[2], FILE_IGNORE_NEW_LINES);

$exist = 0;
if(count($merged) > 2){
	$exist = 1;
}



$mergePoint_ans = array();
$mergePoint_mer = array();
$merge_error = 0;

$common_ans = array();
$common_mer = array();
$ans_flag = array();
$mer_flag = array();

$num_ans = 0;
$num_mer = 0;
$num_common = 0;
$num_scaffold = 0;

// find common contigs
foreach($answer as $k => $ans){
	if($ans == "") continue;
	if(strstr($ans, "-----") || strstr($ans, ">")){
	//	array_push($common_ans, "-----");
		continue;
	}
	$tag = explode(" ", $ans);
	$ans_tag = $tag[0];
	foreach($merged as $i => $contig){
		if($contig == "" || strstr($contig, "-----") || strstr($contig, ">")) continue;
		$tag = explode(" ", $contig);
		$merged_tag = $tag[0];
		if($ans_tag == $merged_tag){
		//	$ans_flag[$k] = 1;
		//	$mer_flag[$i] = 1;
		//	array_push($common_ans, $ans);
			$num_common++;
			break;
		}
	}
//	$num_ans++;
}
/*
foreach($merged as $i => $contig){
	if($contig == "") continue;
	if(strstr($contig, "-----") || strstr($contig, ">")){
		array_push($common_mer, "-----");
	//	$num_scaffold++;
		continue;
	}
	if(isset($mer_flag[$i])){
		array_push($common_mer, $contig);
	}
	$num_mer++;
}
 */

foreach($answer as $k => $ans){
	if($ans == "") continue;
	if(strstr($ans, "-----") || strstr($ans, ">")){
		array_push($common_ans, "-----");
		continue;
	}
	array_push($common_ans, $ans);
	$num_ans++;
}
foreach($merged as $i => $contig){
	if($contig == "") continue;
	if(strstr($contig, "-----") || strstr($contig, ">")){
		array_push($common_mer, "-----");
		continue;
	}
	array_push($common_mer, $contig);
	$num_mer++;
}
//if($num_scaffold == 0) $num_scaffold = 1;

if($common_ans[0] != "-----"){
	array_unshift($common_ans, "-----");
}

if(@$common_mer[0] != "-----"){		
	array_unshift($common_mer, "-----");
}

if($common_ans[count($common_ans)-1] != "-----"){
	array_push($common_ans, "-----");
}
if($common_mer[count($common_mer)-1] != "-----"){		
	array_push($common_mer, "-----");
}

// Count the #scaffold in answer and result
$common_ans_scaffold = 0;
foreach($common_ans as $k => $each){
	if(!isset($common_ans[$k+1])) continue;
	if($each == "-----" && $common_ans[$k+1] != "-----"){
		$common_ans_scaffold++;
	}
}

$common_mer_scaffold = 0;
$singleton = 0;
foreach($common_mer as $k => $each){
	if(!isset($common_mer[$k+1])) continue;
	if($each == "-----" && $common_mer[$k+1] != "-----"){
		$common_mer_scaffold++;
	}
	if($each == "-----" && $common_mer[$k+1] != "-----" && $common_mer[$k+2] == "-----"){
		$singleton++;
	}
}

if($detail == 1){
	print "common contigs in answer:\n";
	print_r($common_ans);
	print "common contigs in result:\n";
	print_r($common_mer);
}

//================================================
//calculates coverage
//================================================
if($cov == 1){
	if(file_exists($contigMer)){
		$contigsFasta = $contigMer;
	}
	else{
		$contigsFasta = "ref0/$contigMer";
	}
	$seqLen = array();
	$sumLen = 0;
	$content = file_get_contents($contigsFasta);
	
	$lines = explode(">", $content);
	foreach($lines as $each){
		if($each == "") continue;
		$tag = substr($each, 0, strpos($each, " "));
		if($tag == ""){
			$tag = substr($each, 0, strpos($each, "\n"));
		}
		$seq = substr($each, strpos($each, "\n"));
		$seq = str_replace("\n", "", $seq);
		$contigLen = strlen($seq);
		$seqLen[$tag] = $contigLen;
	}
	if(0){
		foreach($common_mer as $each){
			if($each == "" || $each == "-----") continue;
			if($detail == 1){
				print "contig: $each\n";
				print "contig length: ".$seqLen[trim(substr($each, 0, -2))]."\n";
			}
			$sumLen += $seqLen[trim(substr($each, 0, -2))];
		}
	}
}

foreach($common_ans as $i => $each){
	//	if($each == "-----") continue;
	if(isset($common_ans[$i+1])){
		$next = $common_ans[$i+1];
		//	if($common_ans[$i+1] == "-----") continue;
		if($cov == 1 && $each != "-----"){
			$tag = substr($each, 0, strpos($each, " "));
			if($tag == ""){
				$tag = substr($each, 0, strpos($each, "\n"));
			}

			if($common_ans[$i-1] == "-----"){
				$sumLen += $seqLen[$tag]/2;
			}
			else if($next == "-----"){
				$sumLen += $seqLen[$tag]/2;
			}
			else{
				$sumLen += $seqLen[$tag];
			}

		}
		array_push($mergePoint_ans, array($each, $next));
	}
}
$num_mergePoint = 0;
foreach($common_mer as $i => $each){
	if($each != "-----" && $common_mer[$i+1] != "-----"){
		$num_mergePoint++;
	}
	if(isset($common_mer[$i+1])){
	//	if($common_mer[$i+1] == "-----") continue;
		array_push($mergePoint_mer, array($each, $common_mer[$i+1]));
	}
}


if($detail == 1){
	print_r($mergePoint_ans);
	print_r($mergePoint_mer);
}



/*
if($detail == 1){
	print "joins:\n";
	print_r($mergePoint_mer);
}
 */

//================================================
$coverLen = 0;
$s_start = 0;
$s_end = 0;
$s_end_single = 0;


// Filter redundant joins in result
$n = count($mergePoint_mer);
for($i = 0; $i < $n; $i++){
	if($mergePoint_mer[$i][0] == "-----" || $mergePoint_mer[$i][1] == "-----"){
		continue;
	}

	$tmp = explode(" ", $mergePoint_mer[$i][0]);
	$reverEach = $tmp[0]." ".(1-trim($tmp[1]));
	$tmp2 = explode(" ", $mergePoint_mer[$i][1]);
	$reverEach2 = $tmp2[0]." ".(1-trim($tmp2[1]));

	for($j = $i+1; $j < $n; $j++){
		if($mergePoint_mer[$i] == $mergePoint_mer[$j]
		|| array($reverEach2, $reverEach) == $mergePoint_mer[$j]){
			print "identical join: $i $j\n";
			print_r($mergePoint_mer[$i]);
			print_r($mergePoint_mer[$j]);
			array_splice($mergePoint_mer, $j, 1);
			$num_mergePoint--;
			$n--;
		}
	}
}
if ($detail == 1){
	print "After filter:\n";
	print_r($mergePoint_mer);
}

if(0){
	$correct_join = 0;
	foreach($mergePoint_mer as $k => $each){
		if($each[0] == "-----" || $each[1] == "-----"){
			continue;
		}

		$tmp = explode(" ", $each[0]);
		$reverEach = $tmp[0]." ".(1-trim($tmp[1]));
		$tmp2 = explode(" ", $each[1]);
		$reverEach2 = $tmp2[0]." ".(1-trim($tmp2[1]));

		if(in_array($each, $mergePoint_ans)
			|| in_array(array($reverEach2, $reverEach), $mergePoint_ans)){
				if($detail == 1){
					print "Correct Join: $each[0] -- $each[1]\n";
				}
				$correct_join++;	
			}
	}
}

foreach($mergePoint_mer as $k => $each){
	if($each[0] == "-----"){
		$s_start = 1;
		continue;
	}
	if($each[1] == "-----"){
		if($s_start == 1){
			$s_end_single = 1;
			continue;
		}
		else{
			continue;
		}
	}

	if(isset($mergePoint_mer[$k+1]) && $mergePoint_mer[$k+1][1] == "-----"){
		$s_end = 1;
	}
	
	if($each[0] != "-----"){
		$tmp = explode(" ", $each[0]);
		$reverEach = $tmp[0]." ".(1-trim($tmp[1]));
	}
	if($each[1] != "-----"){
		$tmp2 = explode(" ", $each[1]);
		$reverEach2 = $tmp2[0]." ".(1-trim($tmp2[1]));
	}

	if($each[0] != "-----" && $each[1] != "-----"
		&& !in_array($each, $mergePoint_ans)
		&& !in_array(array($reverEach2, $reverEach), $mergePoint_ans)){
		if($detail == 1){
			print "Mis-Join: $each[0] -- $each[1]\n";
		}
		$merge_error++;	
		$s_start = 0;
		$s_end = 0;
		$s_end_single = 0;
	}
	else{
		if($cov == 1){
			if(0){
				if($s_start == 1 && $s_end == 1){
					//	$len_tmp = $seqLen[$tmp[0]] + $seqLen[$tmp2[0]];
					//	$coverLen += $len_tmp;
					$s_start = 0;
					$s_end = 0;
					if($detail == 1){
						print "Start and End point\n";
					}
				}
				else if($s_start == 1 && $s_end_single == 1){
					//	$len_tmp = $seqLen[$tmp[0]];
					//	$coverLen += $len_tmp;
					$s_start = 0;
					$s_end_single = 0;
					if($detail == 1){
						print "Single contig\n";
					}
				}
				else if($s_start == 1){
					$len_tmp = $seqLen[$tmp[0]] + ($seqLen[$tmp2[0]]/2);
					$coverLen += $len_tmp;
					$s_start = 0;
					if($detail == 1){
						print "Start point\n";
					}
				}
				else if($s_end == 1){
					$len_tmp = ($seqLen[$tmp[0]]/2) + $seqLen[$tmp2[0]];
					$coverLen += $len_tmp;
					$s_end = 0;
					if($detail == 1){
						print "End point\n";
					}
				}
				else{
					$len_tmp = ($seqLen[$tmp[0]]/2) + ($seqLen[$tmp2[0]]/2);
					$coverLen += $len_tmp;
					if($detail == 1){
						print "Normal point\n";
					}
				}
			}
			$len_tmp = ($seqLen[$tmp[0]]/2) + ($seqLen[$tmp2[0]]/2);
			$coverLen += $len_tmp;
			
			if($detail == 1){
				print_r($each);
				print "len_tmp: $len_tmp\n";
				print "coverLen: $coverLen\n";
			}
		}
	//	print "cover: $coverLen\n";
	}
}

if($cov == 1){
	if($exist == 0){
		$coverage = 0;
	}
	else{
		$coverage = $coverLen/$sumLen;
	}
}

//$join_in_answer = $num_common-$common_ans_scaffold;
$join_in_answer = ($num_ans-1)-($common_ans_scaffold-1);
//$unused_contigs = $num_ans-$num_mer;
$unused_contigs = $num_ans-$num_common;
$multiContig_scaffold = $common_mer_scaffold-$singleton;
$all_scaffold = $common_mer_scaffold+$unused_contigs;


if($detail == 1){
	print "\n";
	print "Number of contigs in answer: $num_ans\n";
	print "Number of contigs in result: $num_mer\n";
	print "\n";
	print "After intersecting the contigs ...\n";
	print "Number of scaffolds in answer: $common_ans_scaffold\n";
	print "Number of scaffolds in result: $common_mer_scaffold\n";
	print "Number of contigs in common: $num_common\n";
	print "Number of unused contigs: $unused_contigs\n";
	print "Number of joins in answer: $join_in_answer\n";
	print "Number of joins in result: $num_mergePoint\n";
	print "Number of mis-joins: $merge_error\n";
//	print "Number of correct-joins: $correct_join\n";
}

//$join_num = $num_common-1;
//if($join_num <=0) $join_num = "N";

//print $merge_error."	".$num_mergePoint."\n";
//print $merge_error."	".$join_num."	".$num_mergePoint."\n";
//print $merge_error."	".$join_num."\n";

//@file_put_contents($argv[3], $merge_error."	".$join_num."	".$num_mergePoint."\n", FILE_APPEND);
//@file_put_contents($argv[3], $merge_error."	".$join_num."	".$coverage."\n", FILE_APPEND);
$sens = ($num_mergePoint - $merge_error) / $join_in_answer;

if($num_mergePoint > 0)
	$prec = ($num_mergePoint - $merge_error) / $num_mergePoint;
else
	$prec = "N/A";

if($sens + $prec > 0)
	$Fscore = 2 * $sens * $prec / ($sens + $prec);
else
	$Fscore = "N/A";

if ($nga == 1){
	include getenv("HOME")."/bin/parse_quast_report.php";
}

if($fw == 1){
	if($cov == 1){
		file_put_contents($outputFile, $merge_error."	".$num_mergePoint."	".$coverage."\n", FILE_APPEND);
	}
	else{
		file_put_contents($outputFile, $merge_error."	".$num_mergePoint."	"." "."\n", FILE_APPEND);
	}
}
if($cov == 1){
	if($detail == 1){
		print "Sum of contig length in result: $sumLen\n";
		print "Sum of contig length of correct-joins in result: $coverLen\n";
		print "Coverage: $coverage\n";
	}
//	print "<#MisJoins> <#Joins> <#Coverage> <#Joins in Answer>\n";
//	print "$merge_error\t$num_mergePoint\t$coverage\t$join_in_answer\n";
//	print "<#Scaffolds> <#MisJoins> <#Joins> <#Coverage> <#Joins in Answer>\n";
//	print "$common_mer_scaffold\t$merge_error\t$num_mergePoint\t$coverage\t$join_in_answer\n";
	if (isset($NGA50)){
		print "$sens\t$prec\t$Fscore\t$coverage\t$NGA50\t$all_scaffold\t$multiContig_scaffold\t$merge_error\t$num_mergePoint\n";
	}
	else{
	//	print "$all_scaffold\t$multiContig_scaffold\t$merge_error\t$num_mergePoint\t$sens\t$prec\t$Fscore\t$coverage\n";
		print "$sens\t$prec\t$Fscore\t$coverage\tN/A\t$all_scaffold\t$multiContig_scaffold\t$merge_error\t$num_mergePoint\n";
	}
//	print "$sens\t$prec\t$Fscore\t$coverage\t$all_scaffold\t$multiContig_scaffold\n";
//	print "$all_scaffold\t$multiContig_scaffold\n";
//	@file_put_contents($argv[4], $coverage."\n", FILE_APPEND);
}
else{
//	print "<#MisJoins> <#Joins> <#Joins in Answer>\n";
	//	print "$merge_error\t$num_mergePoint\t$join_in_answer\n";
	//print "<#Scaffolds> <#MisJoins> <#Joins> <#Joins in Answer>\n";
	//print "$common_mer_scaffold\t$merge_error\t$num_mergePoint\t$join_in_answer\n";
	
//	print "<#Scaffolds> <#Multi-contig scaffolds> <#MisJoins> <#Joins> <Sen.> <Pre.> <Fscore>\n";
//	print "$all_scaffold\t$multiContig_scaffold\t$merge_error\t$num_mergePoint\t$sens\t$prec\t$Fscore\n";
	printf("Scaffolds:		%5d\n", $all_scaffold);
	printf("Multi-contig:	%5d\n", $multiContig_scaffold);
	printf("MisJoins:		%5d\n", $merge_error);
	printf("Joins:			%5d\n", $num_mergePoint);
	printf("Sensitivity:	%2.3f\n", $sens);
	printf("Precision:		%2.3f\n", $prec);
	printf("Fscore:			%2.3f\n", $Fscore);

}
//print "Mis-joining rate: ".$merge_error/$num_mergePoint."\n";
?>
