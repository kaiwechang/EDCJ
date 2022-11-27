CC = gcc
CXX = g++
CFLAGS = -O3 -lm
CXXFLAGS = $(CFLAGS) -std=c++11

SERVER_ONLY = -I$(HOME)/gurobi952/linux64/include -L$(HOME)/gurobi952/linux64/lib -lgurobi95
LOCAL_ONLY = -lgurobi100

#TARGETS = speedup_1 speedup_2 speedup_3 markerReorder findContigTelomere ilp ilp_nocap
#TARGETS = findContigTelomere ilp ilp_nocap
TARGETS = findContigTelomere ilp ilp_nocap
BIN_DIR = bin
OUT_DIR = output
#TEST_DIR = testcase/EI_test
TEST_DIR = testcase/simData_smaller/sim1
#TEST_DIR = testcase/simdata1/sim_1000_80_5_100_1_30/1

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

run:
	bin/speedup_1			$(TEST_DIR)/reference.all		$(TEST_DIR)/query.all			$(OUT_DIR)
	bin/markerReorder		$(OUT_DIR)/refmod_spd1.all		$(OUT_DIR)/quemod_spd1.all		$(OUT_DIR)
	bin/speedup_2			$(OUT_DIR)/refmod_spd1_mr.all	$(OUT_DIR)/quemod_spd1_mr.all	$(OUT_DIR)
	bin/markerReorder		$(OUT_DIR)/refmod_spd2.all		$(OUT_DIR)/quemod_spd2.all		$(OUT_DIR)
	bin/speedup_3			$(OUT_DIR)/refmod_spd2_mr.all	$(OUT_DIR)/quemod_spd2_mr.all	$(OUT_DIR)
	bin/ilp					$(OUT_DIR)/refmod_spd3.all		$(OUT_DIR)/quemod_spd3.all		$(OUT_DIR)
	bin/findContigTelomere	$(OUT_DIR)/quemod_spd3.all		$(OUT_DIR)/joins.txt			$(OUT_DIR)
	related_software/misJoin_intersect_v7/misJoin_intersect_v7.php $(TEST_DIR)/answerToAll $(OUT_DIR)/myScaffold.txt > $(OUT_DIR)/evaulate.txt

run_pure:
	bin/ilp					$(TEST_DIR)/reference.all		$(TEST_DIR)/query.all 			$(OUT_DIR)	| tee $(OUT_DIR)/stdout.log
	bin/findContigTelomere	$(TEST_DIR)/query.all			$(OUT_DIR)/joins.txt			$(OUT_DIR)
	related_software/misJoin_intersect_v7/misJoin_intersect_v7.php $(TEST_DIR)/answerToAll $(OUT_DIR)/myScaffold.txt > $(OUT_DIR)/evaulate.txt

run_nocap:
	bin/ilp_nocap			$(TEST_DIR)/reference.all 		$(TEST_DIR)/query.all			$(OUT_DIR) | tee $(OUT_DIR)/stdout.log
	bin/findContigTelomere	$(TEST_DIR)/query.all			$(OUT_DIR)/joins.txt			$(OUT_DIR)
	related_software/misJoin_intersect_v7/misJoin_intersect_v7.php $(TEST_DIR)/answerToAll $(OUT_DIR)/myScaffold.txt > $(OUT_DIR)/evaulate.txt

cycle_test:
	$(eval test_base="testcase/simData_smaller")
	$(eval out_base="output")
	for dir in $$(ls $(test_base)); do		\
		nocap_dir=$(out_base)/$$dir/nocap;	\
		mkdir -p $$nocap_dir;				\
		time -f %e -o $$nocap_dir/time.txt	\
		make run_nocap OUT_DIR=$$nocap_dir	\
		TEST_DIR=$(test_base)/$$dir;		\
		pure_dir=$(out_base)/$$dir/pure;	\
		mkdir -p $$pure_dir;				\
		time -f %e -o $$pure_dir/time.txt	\
		make run_pure OUT_DIR=$$pure_dir	\
		TEST_DIR=$(test_base)/$$dir;		\
	done
	#time -f %e make run_nocap | tee $(OUT_DIR)/time.txt 
#	tick=$$(date +%s); \
#	tock=$$(date +%s); make run_nocap; \
#	echo "time: $$((tock-tick)) sec" | tee $(OUT_DIR)/time.txt

.PHONY: clean
clean:
	-rm -rf $(BIN_DIR) $(OUT_DIR)
