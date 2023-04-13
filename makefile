CXX = g++
CXXFLAGS = -O3 -m64 -std=c++20 -lm -lfmt
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
EBD_Scaffolder = ../senior/EBD_Scaffolder/EBD_Scaffolder

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
	cd ../senior/EBD_Scaffolder; make

.PHONY: all clean mkdir experiment run_* gen_*

run_main: DCJ_Scaffolder
	./DCJ_Scaffolder -r $(TEST_DIR)/reference.all -t $(TEST_DIR)/target.all -o $(OUT_DIR)
	$(TOOL_DIR)/misJoin_eval.php $(TEST_DIR)/answerToAll $(OUT_DIR)/scaffolds.txt	\
	| tee $(OUT_DIR)/evaulate.txt

run_ilp: DCJ_Scaffolder
	./DCJ_Scaffolder -r $(TEST_DIR)/reference.all -t $(TEST_DIR)/query.all -o $(OUT_DIR) -n
	$(TOOL_DIR)/misJoin_eval.php $(TEST_DIR)/answerToAll $(OUT_DIR)/scaffolds.txt	\
	| tee $(OUT_DIR)/evaulate.txt

run_spd3E: DCJ_Scaffolder
	./DCJ_Scaffolder -r $(TEST_DIR)/reference.all -t $(TEST_DIR)/target.all -o $(OUT_DIR) -x
	$(TOOL_DIR)/misJoin_eval.php $(TEST_DIR)/answerToAll $(OUT_DIR)/scaffolds.txt	\
	| tee $(OUT_DIR)/evaulate.txt

run_EBD: $(EBD_Scaffolder) | $(OUT_DIR)
	cd ../senior/EBD_Scaffolder;		\
	./EBD_Scaffolder -s _Sibelia_ -m 70 -i 1800 -e	\
	-cr ../$(TEST_DIR)/reference.all	\
	-ct ../$(TEST_DIR)/target.all		\
	-o ../../code/$(OUT_DIR)/result		\
	 > ../../code/$(OUT_DIR)/EBD_Scaffolder.log
	$(TOOL_DIR)/misJoin_eval.php $(TEST_DIR)/answerToAll $(OUT_DIR)/result/ScaffoldResult	\
	| tee $(OUT_DIR)/evaulate.txt

run_time: DCJ_Scaffolder
	method=main;							\
	time -f %e -o $(OUT_DIR)/time.txt		\
	make run_$$method OUT_DIR=$(OUT_DIR)	\
	TEST_DIR=$(TEST_DIR)

experiment: DCJ_Scaffolder $(EBD_Scaffolder)
	$(eval test_base="../testcase/real_test")
	$(eval out_base=$(OUT_DIR))
	for dir in $$(ls $(test_base)); do						\
		for sub in $$(ls $(test_base)/$$dir); do			\
			for method in EBD main spd3E; do				\
				tar_dir=$(out_base)/$$dir/$$sub/$$method;	\
				mkdir -p $$tar_dir;							\
				time -f %e -o $$tar_dir/time.txt			\
				make run_$$method OUT_DIR=$$tar_dir			\
				TEST_DIR=$(test_base)/$$dir/$$sub;			\
			done											\
		done												\
	done
	$(TOOL_DIR)/print_table $(out_base)

gen_real: $(BIN_DIR)/fna2all
	$(eval test_base="../testcase/real_data")
	for organ in $$(ls $(test_base)); do				\
		tar=$$(ls $(test_base)/$$organ/*.randOrd);		\
		ans=$$(ls $(test_base)/$$organ/answerToAll);	\
		for dir in $$(ls $(test_base)/$$organ); do		\
			if [ -d $(test_base)/$$organ/$$dir ] && [ $$dir != "ext" ]; then	\
				ref=$$(ls $(test_base)/$$organ/$$dir/*.fna);					\
				mkdir -p testcase/$$organ/$$dir;								\
				$(BIN_DIR)/fna2all $$ref $$tar testcase/$$organ/$$dir/sibelia;	\
				cp $$ans testcase/$$organ/$$dir/;								\
				mv testcase/$$organ/$$dir/sibelia/reference.all testcase/$$organ/$$dir/reference.all;	\
				mv testcase/$$organ/$$dir/sibelia/target.all	testcase/$$organ/$$dir/target.all;		\
				rm -r testcase/$$organ/$$dir/sibelia;	\
			fi;											\
		done;											\
	done

gen_sim: $(BIN_DIR)/simulator
	$(eval ori_mkr=3000)
	$(eval dup_len=5)
	$(eval evo_num=300)
	$(eval ref_num=1)
	$(eval tar_num=150)
	# <# initial markers> <inverse rate> <duplicate length> <# evolutions> <# ref contigs> <# tar contigs> <output_dir>
	for inv in 10 20 30 40 50 60 70 80 90 100; do	\
		for sub in 1 2 3 4 5; do					\
			test_dir=testcase/sim_$(ori_mkr)_$${inv}_$(dup_len)_$(evo_num)_$(ref_num)_$(tar_num)/$$sub;	\
			mkdir -p $$test_dir;					\
			$(BIN_DIR)/simulator $(ori_mkr) $$inv $(dup_len) $(evo_num) $(ref_num) $(tar_num) $$test_dir;\
			rm $$test_dir/*_process;				\
		done										\
	done;											\

clean:
	@-rm -rf $(BIN_DIR) $(OUT_DIR) testcase DCJ_Scaffolder
