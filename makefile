CC = gcc
CXX = g++
CFLAGS = -O3 -lm
CXXFLAGS = $(CFLAGS) -std=c++11

SERVER_ONLY = -I$(HOME)/gurobi952/linux64/include -L$(HOME)/gurobi952/linux64/lib -lgurobi95
LOCAL_ONLY = -lgurobi100

#TARGETS = speedup_1E speedup_2 speedup_3 findContigTelomere ilp_nocap
TARGETS = speedup_1 speedup_1R speedup_1E ilp_nocap
BIN_DIR = bin
OUT_DIR = output
#TEST_DIR = ../testcase/EI_test
#TEST_DIR = ../testcase/simData_smaller/sim1
TEST_DIR = ../testcase/simdata1/sim_1000_80_5_100_1_30/1

.PHONY: all
all: mkdir $(TARGETS)

mkdir:
	mkdir -p $(BIN_DIR) $(OUT_DIR)
%: %.cpp
	$(CXX) $(CXXFLAGS) $< -o $(BIN_DIR)/$@
ilp: ilp.cpp
	$(CXX) $(CXXFLAGS) $< -o $(BIN_DIR)/$@ -m64 -g $(LOCAL_ONLY) -lgurobi_c++
ilp_nocap: ilp_nocap.cpp
	$(CXX) $(CXXFLAGS) $< -o $(BIN_DIR)/$@ -m64 -g $(LOCAL_ONLY) -lgurobi_c++

.SILENT:

run_cycle: mkdir
	bin/speedup_1E			$(TEST_DIR)/reference.all	$(TEST_DIR)/query.all		$(OUT_DIR)	> $(OUT_DIR)/speedup_1.log
	bin/speedup_2			$(OUT_DIR)/ref_spd1.all		$(OUT_DIR)/tar_spd1.all		$(OUT_DIR)	> $(OUT_DIR)/speedup_2.log
	bin/speedup_3			$(OUT_DIR)/ref_spd2.all		$(OUT_DIR)/tar_spd2.all		$(OUT_DIR)	> $(OUT_DIR)/speedup_3.log
	bin/ilp_nocap			$(OUT_DIR)/ref_spd3.all		$(OUT_DIR)/tar_spd3.all		$(OUT_DIR)	> $(OUT_DIR)/ilp.log
	bin/findContigTelomere	$(OUT_DIR)/tar_spd3.all		$(OUT_DIR)/joins.txt		$(OUT_DIR)	> $(OUT_DIR)/findContigTelomere.log
	./misJoin_eval.php		$(TEST_DIR)/answerToAll		$(OUT_DIR)/myScaffold.txt				> $(OUT_DIR)/evaulate.txt

run_pure: mkdir
	bin/ilp_nocap			$(TEST_DIR)/reference.all	$(TEST_DIR)/query.all		$(OUT_DIR)	> $(OUT_DIR)/ilp.log
	bin/findContigTelomere	$(TEST_DIR)/query.all		$(OUT_DIR)/joins.txt		$(OUT_DIR)	> $(OUT_DIR)/findContigTelomere.log
	./misJoin_eval.php		$(TEST_DIR)/answerToAll		$(OUT_DIR)/myScaffold.txt				> $(OUT_DIR)/evaulate.txt

run_spd1: mkdir
	bin/speedup_1			$(TEST_DIR)/reference.all	$(TEST_DIR)/query.all		$(OUT_DIR)	> $(OUT_DIR)/speedup_1.log
	bin/ilp_nocap			$(OUT_DIR)/ref_spd1.all		$(OUT_DIR)/tar_spd1.all		$(OUT_DIR)	> $(OUT_DIR)/ilp.log
	bin/findContigTelomere	$(OUT_DIR)/tar_spd1.all		$(OUT_DIR)/joins.txt		$(OUT_DIR)	> $(OUT_DIR)/findContigTelomere.log
	./misJoin_eval.php		$(TEST_DIR)/answerToAll		$(OUT_DIR)/myScaffold.txt				> $(OUT_DIR)/evaulate.txt

run_spd1R: mkdir
	bin/speedup_1R			$(TEST_DIR)/reference.all	$(TEST_DIR)/query.all		$(OUT_DIR)	> $(OUT_DIR)/speedup_1.log

run_spd1E: mkdir
	bin/speedup_1E			$(TEST_DIR)/reference.all	$(TEST_DIR)/query.all		$(OUT_DIR)	> $(OUT_DIR)/speedup_1.log
#	bin/ilp_nocap			$(OUT_DIR)/ref_spd1.all		$(OUT_DIR)/tar_spd1.all		$(OUT_DIR)	> $(OUT_DIR)/ilp.log
#	bin/findContigTelomere	$(OUT_DIR)/tar_spd1.all		$(OUT_DIR)/joins.txt		$(OUT_DIR)	> $(OUT_DIR)/findContigTelomere.log
#	./misJoin_eval.php		$(TEST_DIR)/answerToAll		$(OUT_DIR)/myScaffold.txt				> $(OUT_DIR)/evaulate.txt

experiment:
	$(eval test_base="../testcase/simdata1")
	$(eval out_base="output")
	for dir in $$(ls $(test_base)); do						\
		for sub in $$(ls $(test_base)/$$dir); do			\
			for method in pure spd1E cycle; do				\
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
