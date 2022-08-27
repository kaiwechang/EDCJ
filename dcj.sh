reference_fna=$1
query_fna=$2
ans=$3

ts1=$(date +"%s")
./speedup_1 $reference_fna $query_fna
./markerReorder "refmod.all" "quemod.all"
./speedup_2 "refmod.all" "quemod.all"
./markerReorder "refmod.all" "quemod.all"
./speedup_3 "refmod.all" "quemod.all"
./ilp "refmod.all" "quemod.all"
ts2=$(date +"%s")
./findContigTelomere "quemod.all" joins.txt
ts3=$(date +"%s")
echo "time: "$((ts3-ts1))" s "
php misJoin_intersect_v7.php $ans myScaffold.txt
