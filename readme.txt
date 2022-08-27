(1)首先，請先確認相關程式是否正確安裝，如: Gurobi

(2)dcj.sh 包含了一次完整的 scaffolding 處理 (以下需先編譯)
Command: sh dcj.sh [reference genome 的 .fasta 之路徑加檔名] [target genome 的 .fasta 之路徑加檔名] [scaffold 的答案之路徑加檔名]

(3)speedup_1.cpp
論文中的定理一

(4)speedup_2.cpp
論文中的定理二

(5)speedup_3.cpp
論文中的定理三

(6)markerReorder.cpp
為了避免 family 重分配後有 family 沒有成員，將 family 編號重編
(ex: 1 2 4 => 1 2 3)

(7)findContigTelomere.cpp
還原 joins 對應 scaffolds

(2)~(7) 
透過 g++ -std=c++11 xxx.cpp 編譯

(8)ilp.cpp
ilp 主程式
編譯: g++ -std=c++11 -m64 -g ilpnewcons.cpp -o ilpv3 -I/home/帳號/gurobi版本/linux64/include -L/home/帳號/gurobi版本/linux64/lib -lgurobi_c++ -lgurobi91 -lm -O3

(9)misJoin_intersect_v7
評估 scaffolding 的結果


