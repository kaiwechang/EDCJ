CC = gcc
CXX = g++
CFLAGS = -O3 -lm
CXXFLAGS = $(CFLAGS) -std=c++20 -lfmt

LOCAL_LIB = -lgurobi100
SERVER_LIB = -I$(HOME)/gurobi952/linux64/include -L$(HOME)/gurobi952/linux64/lib -lgurobi95
GUROBI_FLAGS = -m64 -g $(LOCAL_LIB) -lgurobi_c++

BIN_DIR = bin
OUT_DIR = output
#TEST_DIR = ../testcase/EI_test
#TEST_DIR = ../testcase/sim_1000/simData_smaller/sim1
TEST_DIR = ../testcase/sim_1000/simdata1/sim_1000_80_5_100_1_30/1

TARGETS = $(patsubst %.cpp, %, $(shell ls *.cpp))

.PHONY: all
all: mkdir $(addprefix $(BIN_DIR)/, $(TARGETS))

mkdir:
	@mkdir -p $(BIN_DIR) $(OUT_DIR)
$(BIN_DIR)/%: %.cpp
	$(CXX) $(CXXFLAGS) $< -o $@ $(GUROBI_FLAGS)

.SILENT:

run_cycle: all
	$(BIN_DIR)/speedup_1E	$(TEST_DIR)/reference.all	$(TEST_DIR)/query.all		$(OUT_DIR)
	$(BIN_DIR)/speedup_3	$(OUT_DIR)/ref_spd1.all		$(OUT_DIR)/tar_spd1.all		$(OUT_DIR)	> $(OUT_DIR)/speedup_3.log
	$(BIN_DIR)/ilp_nocap	$(OUT_DIR)/ref_spd3.all		$(OUT_DIR)/tar_spd3.all		$(OUT_DIR)	> $(OUT_DIR)/ilp.log
	$(BIN_DIR)/postprocess	$(OUT_DIR)/tar_spd3.all		$(OUT_DIR)/joins.txt		$(OUT_DIR)
	./misJoin_eval.php		$(TEST_DIR)/answerToAll		$(OUT_DIR)/scaffolds.txt				> $(OUT_DIR)/evaulate.txt

run_pure: all
	$(BIN_DIR)/ilp_nocap	$(TEST_DIR)/reference.all	$(TEST_DIR)/query.all		$(OUT_DIR)	> $(OUT_DIR)/ilp.log
	$(BIN_DIR)/postprocess	$(OUT_DIR)/query.all		$(OUT_DIR)/joins.txt		$(OUT_DIR)	> $(OUT_DIR)/postprocess.log
	./misJoin_eval.php		$(TEST_DIR)/answerToAll		$(OUT_DIR)/scaffolds.txt				> $(OUT_DIR)/evaulate.txt

run_spd1: all
	$(BIN_DIR)/speedup_1	$(TEST_DIR)/reference.all	$(TEST_DIR)/query.all		$(OUT_DIR)	> $(OUT_DIR)/speedup_1.log
	$(BIN_DIR)/ilp_nocap	$(OUT_DIR)/ref_spd1.all		$(OUT_DIR)/tar_spd1.all		$(OUT_DIR)	> $(OUT_DIR)/ilp.log
	$(BIN_DIR)/postprocess	$(OUT_DIR)/tar_spd1.all		$(OUT_DIR)/joins.txt		$(OUT_DIR)	> $(OUT_DIR)/postprocess.log
	./misJoin_eval.php		$(TEST_DIR)/answerToAll		$(OUT_DIR)/scaffolds.txt				> $(OUT_DIR)/evaulate.txt

run_spd1E: all
	$(BIN_DIR)/speedup_1E	$(TEST_DIR)/reference.all	$(TEST_DIR)/query.all		$(OUT_DIR)	> $(OUT_DIR)/speedup_1.log
	$(BIN_DIR)/ilp_nocap	$(OUT_DIR)/ref_spd1.all		$(OUT_DIR)/tar_spd1.all		$(OUT_DIR)	> $(OUT_DIR)/ilp.log
	$(BIN_DIR)/postprocess	$(OUT_DIR)/tar_spd1.all		$(OUT_DIR)/joins.txt		$(OUT_DIR)	> $(OUT_DIR)/postprocess.log
	./misJoin_eval.php		$(TEST_DIR)/answerToAll		$(OUT_DIR)/scaffolds.txt				> $(OUT_DIR)/evaulate.txt

experiment: all
	$(eval test_base="../testcase/sim_5000")
	$(eval out_base="output")
	for dir in $$(ls $(test_base)); do						\
		for sub in $$(ls $(test_base)/$$dir); do			\
			for method in spd1E; do							\
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
	$(eval ori_mkr=2000)
	$(eval dup_len=5)
	$(eval evo_num=200)
	$(eval ref_num=10)
	$(eval tar_num=100)
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
