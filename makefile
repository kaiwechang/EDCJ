CC = gcc
CXX = g++
CFLAGS = -O3 -lm
CXXFLAGS = $(CFLAGS) -std=c++11

SERVER_ONLY = -I$(HOME)/gurobi952/linux64/include -L$(HOME)/gurobi952/linux64/lib -lgurobi95
LOCAL_ONLY = -lgurobi100

TARGETS = speedup_1 speedup_2 speedup_3 markerReorder findContigTelomere ilp ilp_nocap
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

.SILENT: run_pure run_nocap

run: run_nocap

run_pure: mkdir
	bin/speedup_1			$(TEST_DIR)/reference.all	$(TEST_DIR)/query.all		$(OUT_DIR)	> $(OUT_DIR)/speedup_1.log
	bin/markerReorder		$(OUT_DIR)/ref_spd1.all		$(OUT_DIR)/tar_spd1.all		$(OUT_DIR)	> $(OUT_DIR)/markerReorder.log
	bin/speedup_2			$(OUT_DIR)/ref_spd1.all		$(OUT_DIR)/tar_spd1.all		$(OUT_DIR)	> $(OUT_DIR)/speedup_2.log
	bin/markerReorder		$(OUT_DIR)/ref_spd2.all		$(OUT_DIR)/tar_spd2.all		$(OUT_DIR)	> $(OUT_DIR)/markerReorder.log
	bin/speedup_3			$(OUT_DIR)/ref_spd2.all		$(OUT_DIR)/tar_spd2.all		$(OUT_DIR)	> $(OUT_DIR)/speedup_3.log
	bin/ilp					$(OUT_DIR)/ref_spd3.all		$(OUT_DIR)/tar_spd3.all		$(OUT_DIR)	> $(OUT_DIR)/ilp.log
	bin/findContigTelomere	$(OUT_DIR)/tar_spd3.all		$(OUT_DIR)/joins.txt		$(OUT_DIR)	> $(OUT_DIR)/findContigTelomere.log
	./misJoin_eval.php		$(TEST_DIR)/answerToAll		$(OUT_DIR)/myScaffold.txt				> $(OUT_DIR)/evaulate.txt

run_nocap: mkdir
	bin/speedup_1			$(TEST_DIR)/reference.all	$(TEST_DIR)/query.all		$(OUT_DIR)	> $(OUT_DIR)/speedup_1.log
	bin/markerReorder		$(OUT_DIR)/ref_spd1.all		$(OUT_DIR)/tar_spd1.all		$(OUT_DIR)	> $(OUT_DIR)/markerReorder.log
	bin/speedup_2			$(OUT_DIR)/ref_spd1.all		$(OUT_DIR)/tar_spd1.all		$(OUT_DIR)	> $(OUT_DIR)/speedup_2.log
	bin/markerReorder		$(OUT_DIR)/ref_spd2.all		$(OUT_DIR)/tar_spd2.all		$(OUT_DIR)	> $(OUT_DIR)/markerReorder.log
	bin/speedup_3			$(OUT_DIR)/ref_spd2.all		$(OUT_DIR)/tar_spd2.all		$(OUT_DIR)	> $(OUT_DIR)/speedup_3.log
	bin/ilp_nocap			$(OUT_DIR)/ref_spd3.all		$(OUT_DIR)/tar_spd3.all		$(OUT_DIR)	> $(OUT_DIR)/ilp.log
	bin/findContigTelomere	$(OUT_DIR)/tar_spd3.all		$(OUT_DIR)/joins.txt		$(OUT_DIR)	> $(OUT_DIR)/findContigTelomere.log
	./misJoin_eval.php		$(TEST_DIR)/answerToAll		$(OUT_DIR)/myScaffold.txt				> $(OUT_DIR)/evaulate.txt

cycle_test:
	$(eval test_base="../testcase/rough_test")
	$(eval out_base="output")
	for dir in $$(ls $(test_base)); do				\
		for sub in $$(ls $(test_base)/$$dir); do	\
			nocap_dir=$(out_base)/$$dir/nocap/$$sub;\
			mkdir -p $$nocap_dir;					\
			time -f %e -o $$nocap_dir/time.txt		\
			make run_nocap OUT_DIR=$$nocap_dir		\
			TEST_DIR=$(test_base)/$$dir/$$sub;		\
			pure_dir=$(out_base)/$$dir/pure/$$sub;	\
			mkdir -p $$pure_dir;					\
			time -f %e -o $$pure_dir/time.txt		\
			make run_pure OUT_DIR=$$pure_dir		\
			TEST_DIR=$(test_base)/$$dir/$$sub;		\
		done										\
	done
	#time -f %e make run_nocap | tee $(OUT_DIR)/time.txt 
#	tick=$$(date +%s); \
#	tock=$$(date +%s); make run_nocap; \
#	echo "time: $$((tock-tick)) sec" | tee $(OUT_DIR)/time.txt

.PHONY: clean
clean:
	-rm -rf $(BIN_DIR) $(OUT_DIR)
