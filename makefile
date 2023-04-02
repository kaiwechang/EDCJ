CXX = g++
CXXFLAGS = -O3 -m64 -std=c++20 -lm -lfmt
MAKEFLAGS = -j $$(nproc --all)
GUROBI_FLAGS = -I$(GUROBI_HOME)/include -L$(GUROBI_HOME)/lib -lgurobi_c++ -lgurobi100

BIN_DIR = bin
OUT_DIR = output
TOOL_DIR = tools
#TEST_DIR = ../testcase/EI_test
#TEST_DIR = ../testcase/sim_smaller/sim1
#TEST_DIR = ../testcase/sim_1000/simdata1/sim_1000_30_5_100_1_30/1
TEST_DIR = ../testcase/sim_1000/simdata1/sim_1000_30_5_100_1_30/5
#TEST_DIR = ../testcase/sim_2000_50/sim_2000_30_5_200_10_50/2

TARGETS = $(patsubst %.cpp, %, $(shell ls *.cpp))
OBJECTS = $(patsubst %, $(BIN_DIR)/%.o, $(TARGETS))

.PHONY: all
all: mkdir DCJ_Scaffolder

mkdir:
	@mkdir -p $(BIN_DIR) $(OUT_DIR)
$(BIN_DIR)/%.o: %.cpp utils.h
	$(CXX) -c $< -o $@ $(CXXFLAGS) $(GUROBI_FLAGS)
DCJ_Scaffolder: $(OBJECTS)
	$(CXX) -o $@ $^ $(CXXFLAGS) $(GUROBI_FLAGS)

.SILENT:

run_main: all
	./DCJ_Scaffolder -r $(TEST_DIR)/reference.all -t $(TEST_DIR)/query.all -o $(OUT_DIR)
	$(TOOL_DIR)/misJoin_eval.php $(TEST_DIR)/answerToAll $(OUT_DIR)/scaffolds.txt	\
	| tee $(OUT_DIR)/evaulate.txt

run_ilp: all
	./DCJ_Scaffolder -r $(TEST_DIR)/reference.all -t $(TEST_DIR)/query.all -o $(OUT_DIR) -n
	$(TOOL_DIR)/misJoin_eval.php $(TEST_DIR)/answerToAll $(OUT_DIR)/scaffolds.txt	\
	| tee $(OUT_DIR)/evaulate.txt

run_hack: all
	./DCJ_Scaffolder -r $(TEST_DIR)/reference.all -t $(TEST_DIR)/query.all -o $(OUT_DIR) -ea
	$(TOOL_DIR)/misJoin_eval.php $(TEST_DIR)/answerToAll $(OUT_DIR)/scaffolds.txt	\
	| tee $(OUT_DIR)/evaulate.txt

run_EBD: mkdir
	cd ../related_software/EBD_Scaffolder;			\
	./EBD_Scaffolder -s _Sibelia_ -m 70 -i 1800 -e	\
	-cr ../$(TEST_DIR)/reference.all				\
	-ct ../$(TEST_DIR)/query.all					\
	-o ../../src/$(OUT_DIR)/result					\
	 > ../../src/$(OUT_DIR)/EBD_Scaffolder.log
	$(TOOL_DIR)/misJoin_eval.php $(TEST_DIR)/answerToAll $(OUT_DIR)/result/ScaffoldResult	\
	| tee $(OUT_DIR)/evaulate.txt

run_time: all
	method=main;							\
	time -f %e -o $(OUT_DIR)/time.txt		\
	make run_$$method OUT_DIR=$(OUT_DIR)	\
	TEST_DIR=$(TEST_DIR)

experiment: all
	$(eval test_base="../testcase/sim_1000/simdata1")
	$(eval out_base="output")
	for dir in $$(ls $(test_base)); do						\
		for sub in $$(ls $(test_base)/$$dir); do			\
			for method in EBD main; do					\
				tar_dir=$(out_base)/$$dir/$$sub/$$method;	\
				mkdir -p $$tar_dir;							\
				time -f %e -o $$tar_dir/time.txt			\
				make run_$$method OUT_DIR=$$tar_dir			\
				TEST_DIR=$(test_base)/$$dir/$$sub;			\
			done											\
		done												\
	done
	$(TOOL_DIR)/print_table $(out_base)

gen_test:
	$(eval ori_mkr=3000)
	$(eval dup_len=5)
	$(eval evo_num=300)
	$(eval ref_num=1)
	$(eval tar_num=150)
	$(eval test_base="testcase/")
	# <# initial markers> <inverse rate> <duplicate length> <# evolutions> <# ref contigs> <# tar contigs> <output_dir>
	cd $(TOOL_DIR);									\
	g++ simulator.cpp -o simulator;					\
	for inv in 10 20 30 40 50 60 70 80 90 100; do	\
		for sub in 1 2 3 4 5; do					\
			test_dir=$(test_base)/sim_$(ori_mkr)_$${inv}_$(dup_len)_$(evo_num)_$(ref_num)_$(tar_num)/$$sub;	\
			mkdir -p $$test_dir;					\
			./simulator $(ori_mkr) $${inv} $(dup_len) $(evo_num) $(ref_num) $(tar_num)	$$test_dir;	\
			rm $$test_dir/*_process;				\
		done										\
	done;											\
	rm -f simulator

.PHONY: clean
clean:
	@-rm -rf $(BIN_DIR) $(OUT_DIR) DCJ_Scaffolder
