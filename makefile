CC = gcc
CXX = g++
CFLAGS = -O3 -lm
CXXFLAGS = $(CFLAGS) -std=c++20 -lfmt
MAKEFLAGS = -j $$(nproc --all)
GUROBI_FLAGS = -m64 -g -I$(GUROBI_HOME)/include -L$(GUROBI_HOME)/lib -lgurobi_c++ -lgurobi100

BIN_DIR = bin
OUT_DIR = output
#TEST_DIR = ../testcase/EI_test
#TEST_DIR = ../testcase/sim_smaller/sim1
#TEST_DIR = ../testcase/sim_1000/simdata1/sim_1000_30_5_100_1_30/1
TEST_DIR = ../testcase/sim_1000/simdata1/sim_1000_30_5_100_1_30/5
#TEST_DIR = ../testcase/sim_2000_50/sim_2000_30_5_200_10_50/2

TARGETS = $(patsubst %.cpp, %, $(shell ls *.cpp))

.PHONY: all
all: mkdir $(addprefix $(BIN_DIR)/, $(TARGETS))

mkdir:
	@mkdir -p $(BIN_DIR) $(OUT_DIR)
$(BIN_DIR)/%: %.cpp
	$(CXX) $< -o $@ $(CXXFLAGS) $(GUROBI_FLAGS)

.SILENT:

run_main: all
	$(BIN_DIR)/main -r $(TEST_DIR)/reference.all -t $(TEST_DIR)/query.all	\
					-m -i -e -x -p 16 -l 1800 -g 0.001 #-h

run_time: all
	method=ilp;							\
	time -f %e -o $(OUT_DIR)/time.txt		\
	make run_$$method OUT_DIR=$(OUT_DIR)	\
	TEST_DIR=$(TEST_DIR)

run_ilp: all
	$(BIN_DIR)/speedup_1E	$(TEST_DIR)/reference.all	$(TEST_DIR)/query.all		$(OUT_DIR)
	$(BIN_DIR)/speedup_3E	$(OUT_DIR)/ref_spd1.all		$(OUT_DIR)/tar_spd1.all		$(OUT_DIR)
	$(BIN_DIR)/ilp_R		$(OUT_DIR)/ref_spd3.all		$(OUT_DIR)/tar_spd3.all		$(OUT_DIR)
	$(BIN_DIR)/postprocess	$(OUT_DIR)/ref_spd3.all		$(OUT_DIR)/tar_spd3.all		$(OUT_DIR)
	./misJoin_eval.php		$(TEST_DIR)/answerToAll		$(OUT_DIR)/scaffolds.txt	\
	| tee $(OUT_DIR)/evaulate.txt

run_cycle: all
	$(BIN_DIR)/speedup_1E	$(TEST_DIR)/reference.all	$(TEST_DIR)/query.all		$(OUT_DIR)
	$(BIN_DIR)/speedup_3E	$(OUT_DIR)/ref_spd1.all		$(OUT_DIR)/tar_spd1.all		$(OUT_DIR)
	$(BIN_DIR)/ilp_nocap	$(OUT_DIR)/ref_spd3.all		$(OUT_DIR)/tar_spd3.all		$(OUT_DIR)	> $(OUT_DIR)/ilp.log
	$(BIN_DIR)/postprocess	$(OUT_DIR)/ref_spd3.all		$(OUT_DIR)/tar_spd3.all		$(OUT_DIR)
	./misJoin_eval.php		$(TEST_DIR)/answerToAll		$(OUT_DIR)/scaffolds.txt				> $(OUT_DIR)/evaulate.txt

run_hack: all
	$(BIN_DIR)/speedup_1E	$(TEST_DIR)/reference.all	$(TEST_DIR)/query.all		$(OUT_DIR)
	$(BIN_DIR)/speedup_3ER	$(OUT_DIR)/ref_spd1.all		$(OUT_DIR)/tar_spd1.all		$(OUT_DIR)
	$(BIN_DIR)/ilp_nocap	$(OUT_DIR)/ref_spd3.all		$(OUT_DIR)/tar_spd3.all		$(OUT_DIR)	> $(OUT_DIR)/ilp.log
	$(BIN_DIR)/postprocess	$(OUT_DIR)/ref_spd3.all		$(OUT_DIR)/tar_spd3.all		$(OUT_DIR)
	./misJoin_eval.php		$(TEST_DIR)/answerToAll		$(OUT_DIR)/scaffolds.txt				> $(OUT_DIR)/evaulate.txt

run_spd3: all
	$(BIN_DIR)/speedup_1E	$(TEST_DIR)/reference.all	$(TEST_DIR)/query.all		$(OUT_DIR)
	$(BIN_DIR)/speedup_3E	$(OUT_DIR)/ref_spd1.all		$(OUT_DIR)/tar_spd1.all		$(OUT_DIR)
	$(BIN_DIR)/ilp_nocap	$(OUT_DIR)/ref_spd3.all		$(OUT_DIR)/tar_spd3.all		$(OUT_DIR)	> $(OUT_DIR)/ilp.log
	$(BIN_DIR)/postprocess	$(OUT_DIR)/ref_spd3.all		$(OUT_DIR)/tar_spd3.all		$(OUT_DIR)
	./misJoin_eval.php		$(TEST_DIR)/answerToAll		$(OUT_DIR)/scaffolds.txt				> $(OUT_DIR)/evaulate.txt

run_spd3E: all
	$(BIN_DIR)/speedup_1E	$(TEST_DIR)/reference.all	$(TEST_DIR)/query.all		$(OUT_DIR)
	$(BIN_DIR)/speedup_3E	$(OUT_DIR)/ref_spd1.all		$(OUT_DIR)/tar_spd1.all		$(OUT_DIR)	extended
	$(BIN_DIR)/ilp_nocap	$(OUT_DIR)/ref_spd3.all		$(OUT_DIR)/tar_spd3.all		$(OUT_DIR)	> $(OUT_DIR)/ilp.log
	$(BIN_DIR)/postprocess	$(OUT_DIR)/ref_spd3.all		$(OUT_DIR)/tar_spd3.all		$(OUT_DIR)
	./misJoin_eval.php		$(TEST_DIR)/answerToAll		$(OUT_DIR)/scaffolds.txt				> $(OUT_DIR)/evaulate.txt

run_EBD: mkdir
	cd ../related_software/EBD_Scaffolder;			\
	./EBD_Scaffolder -s _Sibelia_ -m 70 -i 1800 -e	\
	-cr ../$(TEST_DIR)/reference.all				\
	-ct ../$(TEST_DIR)/query.all					\
	-o ../../src/$(OUT_DIR)/result					\
	 > ../../src/$(OUT_DIR)/EBD_Scaffolder.log
	./misJoin_eval.php		$(TEST_DIR)/answerToAll		$(OUT_DIR)/result/ScaffoldResult		> $(OUT_DIR)/evaulate.txt

experiment: all
	$(eval test_base="../testcase/sim_3000")
	$(eval out_base="output")
	for dir in $$(ls $(test_base)); do						\
		for sub in $$(ls $(test_base)/$$dir); do			\
			for method in cycle ilp; do						\
				tar_dir=$(out_base)/$$dir/$$sub/$$method;	\
				mkdir -p $$tar_dir;							\
				time -f %e -o $$tar_dir/time.txt			\
				make run_$$method OUT_DIR=$$tar_dir			\
				TEST_DIR=$(test_base)/$$dir/$$sub;			\
			done											\
		done												\
	done
	./print_table $(out_base)

gen_test: all
	$(eval ori_mkr=3000)
	$(eval dup_len=5)
	$(eval evo_num=300)
	$(eval ref_num=1)
	$(eval tar_num=150)
	$(eval test_base="testcase/")
	# <# initial markers> <inverse rate> <duplicate length> <# evolutions> <# ref contigs> <# tar contigs> <output_dir>
	for inv in 10 20 30 40 50 60 70 80 90 100; do	\
		for sub in 1 2 3 4 5; do					\
			test_dir=$(test_base)/sim_$(ori_mkr)_$${inv}_$(dup_len)_$(evo_num)_$(ref_num)_$(tar_num)/$$sub;	\
			mkdir -p $$test_dir;					\
			$(BIN_DIR)/simulator $(ori_mkr) $${inv} $(dup_len) $(evo_num) $(ref_num) $(tar_num)	$$test_dir;	\
			rm $$test_dir/*_process;				\
		done										\
	done

.PHONY: clean
clean:
	-rm -rf $(BIN_DIR) $(OUT_DIR)
