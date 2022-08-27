CC = gcc
CXX = g++
CFLAGS = -O3 -lm
CXXFLAGS = $(CFLAGS) -std=c++11

TARGETS = speedup_1 speedup_2 speedup_3 markerReorder findContigTelomere ilp
BIN_DIR = bin
OUT_DIR = output
TEST_DIR = testcase/simData_平常實驗用/sim1
#TEST_DIR = testcase/simdata1/sim_1000_0_5_100_1_30/1

.PHONY: all
all: mkdir $(TARGETS)
mkdir:
	mkdir -p $(BIN_DIR) $(OUT_DIR)
%: %.cpp
	$(CXX) $(CXXFLAGS) $< -o $(BIN_DIR)/$@
ilp: ilp.cpp
	$(CXX) $(CXXFLAGS) $< -o $(BIN_DIR)/$@ -m64 -g -lgurobi_c++ -lgurobi95
asdf:
	echo asdf
run:
	bin/speedup_1 $(TEST_DIR)/reference.all $(TEST_DIR)/query.all $(OUT_DIR)
	bin/markerReorder $(OUT_DIR)/refmod_spd1.all $(OUT_DIR)/quemod_spd1.all $(OUT_DIR)
	bin/speedup_2 $(OUT_DIR)/refmod_spd1_mr.all $(OUT_DIR)/quemod_spd1_mr.all $(OUT_DIR)
	bin/markerReorder $(OUT_DIR)/refmod_spd2.all $(OUT_DIR)/quemod_spd2.all $(OUT_DIR)
	bin/speedup_3 $(OUT_DIR)/refmod_spd2_mr.all $(OUT_DIR)/quemod_spd2_mr.all $(OUT_DIR)
	bin/ilp $(OUT_DIR)/refmod_spd3.all $(OUT_DIR)/quemod_spd3.all
	bin/findContigTelomere $(OUT_DIR)/quemod_spd3.all $(OUT_DIR)/joins.txt $(OUT_DIR)
	related_software/misJoin_intersect_v7/misJoin_intersect_v7.php $(TEST_DIR)/answerToAll $(OUT_DIR)/myScaffold.txt

test:
	./dcj.sh $(TEST_DIR)/reference.all $(TEST_DIR)/query.all $(TEST_DIR)/answerToAll
.PHONY: clean
clean:
	-rm -rf $(TARGETS) $(BIN_DIR) $(OUT_DIR)
