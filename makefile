SHELL = /usr/bin/bash	# for array

CXX = g++
CXXFLAGS = -O3 -m64 -std=c++20 -lm -lfmt -fmax-errors=1
MAKEFLAGS = -j $$(nproc --all)
GUROBI_FLAGS = -I$(GUROBI_HOME)/include -L$(GUROBI_HOME)/lib -lgurobi_c++ -lgurobi100

SRC_DIR = src
BIN_DIR = bin
OUT_DIR = output
TOOL_DIR = tools
#TEST_DIR = ../testcase/EI_test
#TEST_DIR = ../testcase/sim_smaller/sim1
TEST_DIR = ../testcase/sim_1000/simdata1/sim_1000_30_5_100_1_30/5

TOOLS = $(patsubst $(TOOL_DIR)/%.cpp, $(BIN_DIR)/%, $(shell ls $(TOOL_DIR)/*.cpp))
OBJECTS = $(patsubst $(SRC_DIR)/%.cpp, $(BIN_DIR)/%.o, $(shell ls $(SRC_DIR)/*.cpp))
EBD_Scaffolder = ../related/EBD_Scaffolder/EBD_Scaffolder

.SILENT:

DCJ_Scaffolder: $(OBJECTS)
	$(CXX) -o $@ $^ $(CXXFLAGS) $(GUROBI_FLAGS)
$(BIN_DIR)/%.o: $(SRC_DIR)/%.cpp $(SRC_DIR)/utils.h | $(BIN_DIR)
	$(CXX) -c $< -o $@ $(CXXFLAGS) $(GUROBI_FLAGS)
$(BIN_DIR)/%: $(TOOL_DIR)/%.cpp | $(BIN_DIR)
	$(CXX) $< -o $@ $(CXXFLAGS)
$(BIN_DIR):
	mkdir -p $@
$(OUT_DIR):
	mkdir -p $@
$(EBD_Scaffolder):
	cd ../related/EBD_Scaffolder; make

.PHONY: all clean mkdir experiment analyze run_* gen_*

run_EDCJ: DCJ_Scaffolder
	./DCJ_Scaffolder -r $(TEST_DIR)/reference.all -t $(TEST_DIR)/target.all -o $(OUT_DIR)
	$(TOOL_DIR)/misJoin_eval $(TEST_DIR)/answerToAll $(OUT_DIR)/scaffolds.txt	\
	| tee $(OUT_DIR)/evaulate.txt

run_ilp: DCJ_Scaffolder
	./DCJ_Scaffolder -r $(TEST_DIR)/reference.all -t $(TEST_DIR)/query.all -o $(OUT_DIR) -n
	$(TOOL_DIR)/misJoin_eval $(TEST_DIR)/answerToAll $(OUT_DIR)/scaffolds.txt	\
	| tee $(OUT_DIR)/evaulate.txt

run_spd2E: DCJ_Scaffolder
	./DCJ_Scaffolder -r $(TEST_DIR)/reference.all -t $(TEST_DIR)/target.all -o $(OUT_DIR) -x
	$(TOOL_DIR)/misJoin_eval $(TEST_DIR)/answerToAll $(OUT_DIR)/scaffolds.txt	\
	| tee $(OUT_DIR)/evaulate.txt

run_EBD: $(EBD_Scaffolder) | $(OUT_DIR)
	cd ../related/EBD_Scaffolder;		\
	./EBD_Scaffolder -s _Sibelia_ -m 70 -i 1800 -e	\
	-cr ../$(TEST_DIR)/reference.all	\
	-ct ../$(TEST_DIR)/target.all		\
	-o ../../code/$(OUT_DIR)/result		\
	 > ../../code/$(OUT_DIR)/EBD_Scaffolder.log
	$(TOOL_DIR)/misJoin_eval $(TEST_DIR)/answerToAll $(OUT_DIR)/result/ScaffoldResult	\
	| tee $(OUT_DIR)/evaulate.txt

run_time: DCJ_Scaffolder | $(OUT_DIR)
	method=main;								\
	/usr/bin/time -f %e -o $(OUT_DIR)/time.txt	\
	make run_$$method OUT_DIR=$(OUT_DIR)		\
	TEST_DIR=$(TEST_DIR)

experiment: DCJ_Scaffolder $(EBD_Scaffolder)
	$(eval test_base="../testcase/ALIGN_50")
	for dir in $$(ls $(test_base)); do				\
		for sub in $$(ls $(test_base)/$$dir); do	\
			out_base=$(OUT_DIR)/$$dir/$$sub;		\
			for method in EBD EDCJ spd2E; do		\
				tar_dir=$$out_base/$$method;		\
				mkdir -p $$tar_dir;					\
				/usr/bin/time -f %e -o $$tar_dir/time.txt	\
				make run_$$method OUT_DIR=$$tar_dir	\
				TEST_DIR=$(test_base)/$$dir/$$sub;	\
			done;									\
			cp $$out_base/EDCJ/DCJ.txt $$out_base;	\
			make analyze OUT_DIR=$$out_base			\
			TEST_DIR=$(test_base)/$$dir/$$sub;		\
		done										\
	done
	# print_table read from analyze_data and gredu and then draw plot

analyze: | $(OUT_DIR)
	#../related/GREDU/bin/dcj $(TEST_DIR)/reference.all $(TEST_DIR)/target.all 100
	$(TOOL_DIR)/analyze_data $(TEST_DIR)/reference.all $(TEST_DIR)/target.all $(OUT_DIR)

GC: $(BIN_DIR)/count_fna
	#$(BIN_DIR)/count_fna ../testcase/LONG_data/NC_007651/ref00/NC_007651.fna
	tar_dir="../testcase/LONG_data/NC_007651";		\
	echo "tar: NC_007651 (draft)";					\
	$(BIN_DIR)/count_fna ../testcase/LONG_data/NC_007651/contigMerged.v3.randOrd;	\
	for dir in ref00 ref08 ref10 ref16 ref18; do	\
		echo $$dir": "$$(find $$tar_dir/$$dir/*.fna -printf "%f\n");				\
		$(BIN_DIR)/count_fna $$tar_dir/$$dir/*.fna;	\
	done;											\
	tar_dir="../testcase/LONG_data/NC_008149";		\
	echo "tar: NC_008149 (draft)";					\
	$(BIN_DIR)/count_fna ../testcase/LONG_data/NC_008149/contigMerged.v3.randOrd;	\
	for dir in ref00 ref09 ref11 ref12 ref13; do	\
		echo $$dir": "$$(find $$tar_dir/$$dir/*.fna -printf "%f\n");				\
		$(BIN_DIR)/count_fna $$tar_dir/$$dir/*.fna;	\
	done

human: $(BIN_DIR)/cut_human
	mkdir -p testcase
	$(BIN_DIR)/cut_human ../testcase/HS_100/HS/ch14 testcase

gen_semi: $(BIN_DIR)/cut_semi
	cmp=../testcase/semi_data/BT/GCA_001718635.1_ASM171863v1_genomic.fna;	\
	for cut in 500; do					\
		mkdir -p testcase/$$cut;			\
		$(BIN_DIR)/cut_semi $$cmp $$cut testcase/$$cut;	\
		mkdir -p testcase/$$cut/ref00;		\
		cp $$cmp testcase/$$cut/ref00;		\
		cmp_base=$$(dirname $$cmp);			\
		for dir in $$(ls $$cmp_base); do	\
			[ -d $$cmp_base/$$dir ] && cp -r $$cmp_base/$$dir testcase/$$cut/$$dir;	\
		done;	\
	done

gen_answer: $(BIN_DIR)/align
	$(eval test_base="../testcase/ALIGN_data")
	for organ in $$(ls $(test_base)); do				\
		comp=$$(ls $(test_base)/$$organ/ref00/*.fna);	\
		draft=$$(ls $(test_base)/$$organ/*.fna);		\
		mkdir -p testcase/$$organ;						\
		$(BIN_DIR)/align $$comp $$draft 90 80 testcase/$$organ;	\
		cp testcase/$$organ/answerToAll $(test_base)/$$organ/;	\
		cp testcase/$$organ/*.randOrd $(test_base)/$$organ/;	\
	done

gen_real: $(BIN_DIR)/fna2all
	$(eval test_base="../testcase/ALIGN_data")
	for organ in $$(ls $(test_base)); do				\
		tar=$$(ls $(test_base)/$$organ/*.randOrd);		\
		ans=$$(ls $(test_base)/$$organ/answerToAll);	\
		for dir in $$(ls $(test_base)/$$organ); do		\
			if [ -d $(test_base)/$$organ/$$dir ] && [ $$dir != "ext" ]; then	\
				ref=$$(ls $(test_base)/$$organ/$$dir/*.fna);					\
				mkdir -p testcase/$$organ/$$dir;								\
				$(BIN_DIR)/fna2all $$ref $$tar 50 testcase/$$organ/$$dir;		\
				cp $$ans testcase/$$organ/$$dir/answerToAll;					\
			fi;											\
		done;											\
	done

gen_sim: $(BIN_DIR)/simulator
	$(eval mkr=2000)
	$(eval dup=5)
	$(eval evo=200)
	$(eval ref=1)
	$(eval tar=100)
	# <# initial markers> <inverse rate> <duplicate length> <# evolutions> <# ref contigs> <# tar contigs> <output_dir>
	for inv in 0 10 20 30 40 50 60 70 80 90 100; do	\
		for sub in 1 2 3 4 5; do					\
			test_dir=testcase/sim_$(mkr)_$${inv}_$(dup)_$(evo)_$(ref)_$(tar)/$$sub;	\
			mkdir -p $$test_dir;					\
			$(BIN_DIR)/simulator $(mkr) $$inv $(dup) $(evo) $(ref) $(tar) $$test_dir;\
			rm $$test_dir/*_process;				\
		done										\
	done;


gen_sim_fix: $(BIN_DIR)/simulator
	$(eval mkr=1000)
	$(eval dup=5)
	$(eval ref=1)
	# <# initial markers> <inverse rate> <duplicate length> <# evolutions> <# ref contigs> <# tar contigs> <output_dir>
	# 1. dup: 100, inv: 0~200, contig: 100;				\
	#M=(100 120 140 160 180 200 220 240 260 280 300);	\
	#R=(  0  17  29  38  44  50  55  58  62  64  67);	\
	#T=(100 100 100 100 100 100 100 100 100 100 100);	\
	# 2. dup: 0~200, inv: 100, contig: 100;				\
	#M=(100 120 140 160 180 200 220 240 260 280 300);	\
	#R=(100  83  71  63  56  50  45  42  38  36  33);	\
	#T=(100 100 100 100 100 100 100 100 100 100 100);	\
	# 3. dup: 100, inv: 100, contig: 50~250;			\
	#M=(200 200 200 200 200 200 200 200 200 200 200);	\
	#R=( 50  50  50  50  50  50  50  50  50  50  50);	\
	#T=( 50  70  90 110 130 150 170 190 210 230 250);	\
	# 4. dup: 50, inv: 0~100, contig: 50;				\
	#M=( 50  60  70  80  90 100 110 120 130 140 150);	\
	#R=(  0  17  29  38  44  50  55  58  62  64  67);	\
	#T=( 50  50  50  50  50  50  50  50  50  50  50);	\
	# 5. dup: 0~100, inv: 50, contig: 50;				\
	M=( 50  60  70  80  90 100 110 120 130 140 150);	\
	R=(100  83  71  63  56  50  45  42  38  36  33);	\
	T=( 50  50  50  50  50  50  50  50  50  50  50);	\
	# 6. dup: 50, inv: 50, contig: 50~150;				\
	#M=(100 100 100 100 100 100 100 100 100 100 100);	\
	#R=( 50  50  50  50  50  50  50  50  50  50  50);	\
	#T=( 50  60  70  80  90 100 110 120 130 140 150);	\
	for idx in {0..10}; do			\
		evo=$${M[$$idx]};			\
		inv=$${R[$$idx]};			\
		tar=$${T[$$idx]};			\
		for sub in 1 2 3 4 5; do	\
			test_dir=testcase/sim_$(mkr)_$${inv}_$(dup)_$${evo}_$(ref)_$$tar/$$sub;\
			mkdir -p $$test_dir;	\
			$(BIN_DIR)/simulator $(mkr) $$inv $(dup) $$evo $(ref) $$tar $$test_dir;\
			rm $$test_dir/*_process;\
		done						\
	done;

clean:
	@-rm -rf $(BIN_DIR) $(OUT_DIR) testcase DCJ_Scaffolder
