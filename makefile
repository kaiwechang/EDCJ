CC = gcc
CXX = g++
CFLAGS = -O3 -lm
CXXFLAGS = $(CFLAGS) -std=c++11 -lfmt

LOCAL_LIB = -lgurobi100
SERVER_LIB = -I$(HOME)/gurobi952/linux64/include -L$(HOME)/gurobi952/linux64/lib -lgurobi95
GUROBI_FLAGS = -m64 -g $(LOCAL_LIB) -lgurobi_c++

BIN_DIR = bin
OUT_DIR = output
#TEST_DIR = ../testcase/EI_test
#TEST_DIR = ../testcase/simData_smaller/sim1
TEST_DIR = ../testcase/simdata1/sim_1000_80_5_100_1_30/1

TARGETS = $(patsubst %.cpp, %, $(shell ls *.cpp))

.PHONY: all
all: mkdir $(addprefix $(BIN_DIR)/, $(TARGETS))

mkdir:
	@mkdir -p $(BIN_DIR) $(OUT_DIR)
$(BIN_DIR)/%: %.cpp
	$(CXX) $(CXXFLAGS) $< -o $@ $(GUROBI_FLAGS)

.SILENT:

run_cycle: mkdir
	$(BIN_DIR)/speedup_1E	$(TEST_DIR)/reference.all	$(TEST_DIR)/query.all		$(OUT_DIR)	> $(OUT_DIR)/speedup_1.log
	$(BIN_DIR)/speedup_3	$(OUT_DIR)/ref_spd1.all		$(OUT_DIR)/tar_spd1.all		$(OUT_DIR)	> $(OUT_DIR)/speedup_3.log
	$(BIN_DIR)/ilp_nocap	$(OUT_DIR)/ref_spd3.all		$(OUT_DIR)/tar_spd3.all		$(OUT_DIR)	> $(OUT_DIR)/ilp.log
	$(BIN_DIR)/postprocess	$(OUT_DIR)/tar_spd3.all		$(OUT_DIR)/joins.txt		$(OUT_DIR)	> $(OUT_DIR)/postprocess.log
	./misJoin_eval.php		$(TEST_DIR)/answerToAll		$(OUT_DIR)/scaffolds.txt				> $(OUT_DIR)/evaulate.txt

run_pure: mkdir
	$(BIN_DIR)/ilp_nocap	$(TEST_DIR)/reference.all	$(TEST_DIR)/query.all		$(OUT_DIR)	> $(OUT_DIR)/ilp.log
	$(BIN_DIR)/postprocess	$(OUT_DIR)/query.all		$(OUT_DIR)/joins.txt		$(OUT_DIR)	> $(OUT_DIR)/postprocess.log
	./misJoin_eval.php		$(TEST_DIR)/answerToAll		$(OUT_DIR)/scaffolds.txt				> $(OUT_DIR)/evaulate.txt

run_spd1: mkdir
	$(BIN_DIR)/speedup_1	$(TEST_DIR)/reference.all	$(TEST_DIR)/query.all		$(OUT_DIR)	> $(OUT_DIR)/speedup_1.log
	$(BIN_DIR)/ilp_nocap	$(OUT_DIR)/ref_spd1.all		$(OUT_DIR)/tar_spd1.all		$(OUT_DIR)	> $(OUT_DIR)/ilp.log
	$(BIN_DIR)/postprocess	$(OUT_DIR)/tar_spd1.all		$(OUT_DIR)/joins.txt		$(OUT_DIR)	> $(OUT_DIR)/postprocess.log
	./misJoin_eval.php		$(TEST_DIR)/answerToAll		$(OUT_DIR)/scaffolds.txt				> $(OUT_DIR)/evaulate.txt

run_spd1E: mkdir
	$(BIN_DIR)/speedup_1E	$(TEST_DIR)/reference.all	$(TEST_DIR)/query.all		$(OUT_DIR)	> $(OUT_DIR)/speedup_1.log
	$(BIN_DIR)/ilp_nocap	$(OUT_DIR)/ref_spd1.all		$(OUT_DIR)/tar_spd1.all		$(OUT_DIR)	> $(OUT_DIR)/ilp.log
	$(BIN_DIR)/postprocess	$(OUT_DIR)/tar_spd1.all		$(OUT_DIR)/joins.txt		$(OUT_DIR)	> $(OUT_DIR)/postprocess.log
	./misJoin_eval.php		$(TEST_DIR)/answerToAll		$(OUT_DIR)/scaffolds.txt				> $(OUT_DIR)/evaulate.txt

experiment:
	$(eval test_base="../testcase/simdata1")
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

.PHONY: clean
clean:
	-rm -rf $(BIN_DIR) $(OUT_DIR)
