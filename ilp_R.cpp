#include "gurobi_c++.h"
#include "utils.h"
#include <utility>
using std::pair;

void readGenome(string refFile, string tarFile, auto& ref, auto& tar, auto& refFamilySize, auto& tarFamilySize, int& maxFamily, auto& refTeloIdx, auto& tarTeloIdx) {
	string contig;
	int id, family, tmp;
	ifstream fin(refFile);
	while (fin >> id >> family >> contig >> tmp) {
		ref.push_back(Marker(id, family, contig));
		maxFamily = max(maxFamily, abs(family));
	}	fin.close();
	fin.open(tarFile);
	while (fin >> id >> family >> contig >> tmp) {
		tar.push_back(Marker(id, family, contig));
		maxFamily = max(maxFamily, abs(family));
	}	fin.close();

	// count family size
	refFamilySize.assign(maxFamily+1, 0);
	tarFamilySize.assign(maxFamily+1, 0);
	for (Marker& m: ref)
		refFamilySize[m.absFamily]++;
	for (Marker& m: tar)
		tarFamilySize[m.absFamily]++;

	// store telo indices
	for (int i = 0; i < ref.size(); i++) {
		if (i == 0 || ref[i-1].contig != ref[i].contig)
			refTeloIdx.push_back(2*i);
		if (i == ref.size()-1 || ref[i].contig != ref[i+1].contig)
			refTeloIdx.push_back(2*i+1);
	}
	for (int i = 0; i < tar.size(); i++) {
		if (i == 0 || tar[i-1].contig != tar[i].contig)
			tarTeloIdx.push_back(2*i);
		if (i == tar.size()-1 || tar[i].contig != tar[i+1].contig)
			tarTeloIdx.push_back(2*i+1);
	}
}
void GRBInit(GRBEnv*& env, GRBModel*& model, string out_dir) {
	env = new GRBEnv(true);
	env->set("LogToConsole", "0");
	env->set("LogFile", out_dir+"/gurobi.log");
	env->set("TimeLimit", "1800");
	env->start();

	model = new GRBModel(env);
}
void matchConstr(auto& ref, auto& tar, auto& refFamilySize, auto& tarFamilySize, int maxFamily, GRBModel* model, auto& refMatchVars, auto& tarMatchVars) {
	// variables (Note: 0 for matched, 1 for unmatched)
	for (Marker& m: ref)
		refMatchVars.push_back(
			refFamilySize[m.absFamily] == 1 ?
			model->addVar(0, 0, 0, GRB_BINARY) :
			model->addVar(0, 1, 0, GRB_BINARY));
	for (Marker& m: tar)
		tarMatchVars.push_back(
			tarFamilySize[m.absFamily] == 1 ?
			model->addVar(0, 0, 0, GRB_BINARY) :
			model->addVar(0, 1, 0, GRB_BINARY));

	// constraints
	vector<GRBLinExpr> refMatchSums(maxFamily+1, 0), tarMatchSums(maxFamily+1, 0);
	for (int i = 0; i < ref.size(); i++)
		refMatchSums[ref[i].absFamily] += (1 - refMatchVars[i]);
	for (int i = 0; i < tar.size(); i++)
		tarMatchSums[tar[i].absFamily] += (1 - tarMatchVars[i]);
	for (int f = 1; f <= maxFamily; f++) {
		// same gene content
		model->addConstr(refMatchSums[f] == tarMatchSums[f]);
		//model->addConstr(refMatchSums[f] == min(refFamilySize[f], tarFamilySize[f]));	// maximum matching
		//model->addConstr(refMatchSums[f] >= 1);										// intermediate
		model->addConstr(refMatchSums[f] == 1);											// exemplar
	}
}
void pairConstr(auto& ref, auto& tar, auto& refFamilySize, auto& tarFamilySize, GRBModel* model, auto& pairVars, auto& refMatchVars, auto& tarMatchVars) {
	// variables
	for (int i = 0; i < ref.size(); i++)
		for (int j = 0; j < tar.size(); j++) {
			if (ref[i].absFamily != tar[j].absFamily)
				continue;
			pairVars[pair(i, j)] = 
				refFamilySize[ref[i].absFamily] == 1 &&
				tarFamilySize[tar[j].absFamily] == 1 ?
				model->addVar(1, 1, 0, GRB_BINARY) :
				model->addVar(0, 1, 0, GRB_BINARY) ;
		}

	// constraints
	vector<GRBLinExpr> refPairSums(ref.size(), 0), tarPairSums(tar.size(), 0);
	for (auto& mp: pairVars) {
		auto& [i, j] = mp.first;
		refPairSums[i] += mp.second;
		tarPairSums[j] += mp.second;
		//fmt::print("({}, {}): {}\n", i, j, 9);
	}
	for (int i = 0; i < ref.size(); i++)
		model->addConstr(refPairSums[i] == 1 - refMatchVars[i]);
	for (int i = 0; i < tar.size(); i++)
		model->addConstr(tarPairSums[i] == 1 - tarMatchVars[i]);
}
void adjConstr(auto& ref, auto& tar, auto& refTeloIdx, auto& tarTeloIdx, GRBModel* model, auto& refRAVars, auto& tarRAVars, auto& refPAVars, auto& tarPAVars) {
	// real adjacencies variables
	for (int i = 0; i < ref.size()-1; i++)
		if (ref[i].contig == ref[i+1].contig)
			refRAVars[pair(2*i+1, 2*i+2)] = model->addVar(1, 1, 0, GRB_BINARY);
	for (int i = 0; i < tar.size()-1; i++)
		if (tar[i].contig == tar[i+1].contig)
			tarRAVars[pair(2*i+1, 2*i+2)] = model->addVar(1, 1, 0, GRB_BINARY);

	// potential adjacencies variables
	for (int i = 0; i < refTeloIdx.size(); i++)
		for (int j = i+1; j < refTeloIdx.size(); j++)
			refPAVars[pair(refTeloIdx[i], refTeloIdx[j])] = model->addVar(0, 1, 0, GRB_BINARY);
	for (int i = 0; i < tarTeloIdx.size(); i++)
		for (int j = i+1; j < tarTeloIdx.size(); j++)
			tarPAVars[pair(tarTeloIdx[i], tarTeloIdx[j])] = model->addVar(0, 1, 0, GRB_BINARY);

	// potential adjacencies constraints
	map<int, GRBLinExpr> refPASums, tarPASums;
	for (auto& [p, v]: refPAVars) {
		auto& [i, j] = p;
		refPASums[i] += v;
		refPASums[j] += v;
	}
	for (auto& [p, v]: tarPAVars) {
		auto& [i, j] = p;
		tarPASums[i] += v;
		tarPASums[j] += v;
	}
	for (auto& [t, v]: refPASums)
		model->addConstr(refPASums[t] == 1);
	for (auto& [t, v]: tarPASums)
		model->addConstr(tarPASums[t] == 1);
}
void cycleNodeConstr(auto& ref, auto& tar, GRBModel* model, auto& refLabelVars, auto& tarLabelVars, auto& refCycleVars, auto& tarCycleVars) {
	// label/cycle variables & upper bound constraints
	int upper = 0;
	for (int i = 0; i < 2*ref.size(); i++) {
		GRBVar reach = model->addVar(0, 1, 0, GRB_BINARY);
		GRBVar label = model->addVar(1, ++upper, 0, GRB_INTEGER);
		//model->addConstr((1 - reach) <= upper - label);
		model->addConstr(reach * upper <= label);
		refLabelVars.push_back(label);
		refCycleVars.push_back(reach);
	}
	for (int i = 0; i < 2*tar.size(); i++) {
		GRBVar reach = model->addVar(0, 1, 0, GRB_BINARY);
		GRBVar label = model->addVar(1, ++upper, 0, GRB_INTEGER);
		//model->addConstr((1 - reach) <= upper - label);
		model->addConstr(reach * upper <= label);
		tarLabelVars.push_back(label);
		tarCycleVars.push_back(reach);
	}
}
void cycleEdgeConstr(auto& ref, auto& tar, GRBModel* model, auto& refMatchVars, auto& tarMatchVars, auto& pairVars, auto& refRAVars, auto& tarRAVars, auto& refPAVars, auto& tarPAVars, auto& refLabelVars, auto& tarLabelVars) {
	auto cycleAdjConstr = [&](auto p, auto& v, auto& labelVars, int offset) {
		auto& [i, j] = p;
		model->addConstr(labelVars[i] <= labelVars[j] + (1 - v) * (offset+i+1));
		model->addConstr(labelVars[j] <= labelVars[i] + (1 - v) * (offset+j+1));
	};
	// RA constraints
	for (auto& [p, v]: refRAVars)	// ref upper bound = idx + 1
		cycleAdjConstr(p, v, refLabelVars, 0);
	for (auto& [p, v]: tarRAVars)	// tar upper bound = ref.size() + idx + 1
		cycleAdjConstr(p, v, tarLabelVars, ref.size());

	// PA constraints
	for (auto& [p, v]: refPAVars)
		cycleAdjConstr(p, v, refLabelVars, 0);
	for (auto& [p, v]: tarPAVars)
		cycleAdjConstr(p, v, tarLabelVars, ref.size());

	// internal edge constraints
	for (int i = 0; i < ref.size(); i++)
		cycleAdjConstr(pair(2*i, 2*i+1), refMatchVars[i], refLabelVars, 0);
	for (int i = 0; i < tar.size(); i++)
		cycleAdjConstr(pair(2*i, 2*i+1), tarMatchVars[i], tarLabelVars, ref.size());

	// pairing constraints
	auto cyclePairConstr = [&](int i, int j, auto& v) {
		model->addConstr(refLabelVars[i] <= tarLabelVars[j] + (1 - v) * (i+1));
		model->addConstr(tarLabelVars[j] <= refLabelVars[i] + (1 - v) * (j+1+ref.size()));
	};
	for (auto& [p, v]: pairVars) {
		auto& [i, j] = p;	// i for ref, j for tar
		if (ref[i].family == tar[j].family) {
			cyclePairConstr(2*i  , 2*j  , v);
			cyclePairConstr(2*i+1, 2*j+1, v);
		} else {
			cyclePairConstr(2*i  , 2*j+1, v);
			cyclePairConstr(2*i+1, 2*j  , v);
		}
	}
}
void objective(GRBModel* model, auto& refCycleVars, auto& tarCycleVars, auto& refMatchVars, auto& tarMatchVars, auto& cycleSum, auto& indelSum) {
	for (auto& v: refCycleVars)	cycleSum += v;
	for (auto& v: tarCycleVars)	cycleSum += v;

	for (auto& v: refMatchVars)	indelSum += v;
	for (auto& v: tarMatchVars)	indelSum += v;

	//model->setObjective(-indelSum-cycleSum);	// intermediate only (?)
	model->setObjective(-cycleSum);
	model->update();
	model->optimize();
}
void outputJoins(string filename, auto& ref, auto& tar, auto& refPAVars, auto& tarPAVars) {
	auto i2t = [](auto& genome, int i) {
		int f = genome[int(i/2)].family;
		int id = genome[int(i/2)].id;
		return f > 0 ^ i%2 == 1 ? -id : id ;
	};

	ofstream fout(filename);
	for (auto& [p, v]: tarPAVars) {
		auto& [i, j] = p;
		if (v.get(GRB_DoubleAttr_X) == 1)
			fout << format("marker {} {}\n", i2t(tar, i), i2t(tar, j));
	}	fout.close();
}
void debug(auto& ref, auto& tar, auto& refTeloIdx, auto& tarTeloIdx, auto& refMatchVars, auto& tarMatchVars, auto& pairVars, auto& refPAVars, auto& tarPAVars, auto& refLabelVars, auto& tarLabelVars, auto& cycleSum, auto& indelSum) {
	logging("basic infos:\n");
	logging("-------+-------+-------+\n");
	logging("   num |   ref |   tar |\n");
	logging("-------+-------+-------+\n");
	logging("marker | {:5d} | {:5d} |\n", ref.size(), tar.size());
	logging("contig | {:5d} | {:5d} |\n", refTeloIdx.size()/2, tarTeloIdx.size()/2);
	logging("-------+-------+-------+\n");
	logging("    PA | {:5d} | {:5d} |\n", refPAVars.size(), tarPAVars.size());
	logging("-------+-------+-------+\n");
	logging("  pair | {:13d} |\n", pairVars.size());
	logging(" cycle | {:13d} |\n", (int)cycleSum.getValue());
	logging(" indel | {:13d} |\n", (int)indelSum.getValue());
	logging("-------+---------------+\n");
	logging("\nref/tar genome:\n");
	logging("-----+--------------------+-------------+\n");
	logging(" idx |  id family  contig | del lab  ub |\n");
	logging("-----+--------------------+-------------+\n");
	for (int i = 0; i < 2*ref.size(); i++)
		logging(" {:3d} | {:3d} {:3d} {:>10s} | {:3d} {:3d} {:3d} |\n", i, ref[i/2].id, ref[i/2].family, ref[i/2].contig, (int)refMatchVars[i/2].get(GRB_DoubleAttr_X), (int)refLabelVars[i].get(GRB_DoubleAttr_X), (int)refLabelVars[i].get(GRB_DoubleAttr_UB));
	logging("-----+--------------------+-------------+\n");
	for (int i = 0; i < 2*tar.size(); i++)
		logging(" {:3d} | {:3d} {:3d} {:>10s} | {:3d} {:3d} {:3d} |\n", i, tar[i/2].id, tar[i/2].family, tar[i/2].contig, (int)tarMatchVars[i/2].get(GRB_DoubleAttr_X), (int)tarLabelVars[i].get(GRB_DoubleAttr_X), (int)tarLabelVars[i].get(GRB_DoubleAttr_UB));
	logging("-----+--------------------+-------------+\n");
	logging("\npairing variables:\n");
	logging("+---------+-----+\n");
	logging("| ref tar | var |\n");
	logging("+---------+-----+\n");
	for (auto& [p, v]: pairVars)
		logging("| {:3d} {:3d} | {:3d} |\n", p.first, p.second, (int)v.get(GRB_DoubleAttr_X));
	logging("+---------+-----+\n");
	logging("\npotential adjacencies:\n");
	logging("+---------+-----+\n");
	logging("| idx jdx | var |\n");
	logging("+---------+-----+\n");
	for (auto& [p, v]: refPAVars)
		logging("| {:3d} {:3d} | {:3d} |\n", p.first, p.second, (int)v.get(GRB_DoubleAttr_X));
	logging("+---------+-----+\n");
	for (auto& [p, v]: tarPAVars)
		logging("| {:3d} {:3d} | {:3d} |\n", p.first, p.second, (int)v.get(GRB_DoubleAttr_X));
	logging("+---------+-----+\n");
}
int main(int argc, char *argv[]) {
	if (argc < 4) {
		fmt::print("[error] Usage:\n>>> ilp <ref genome> <tar genome> <output_dir>\n");
		return 0;
	}	string out_dir(argv[3]);
	logFile.open(out_dir+"/ilp.log");

	// 
	int maxFamily = 0;
	vector<Marker> ref, tar;
	vector<int> refTeloIdx, tarTeloIdx;
	vector<int> refFamilySize, tarFamilySize;

	// gurobi variables
	GRBEnv* env;
	GRBModel* model;
	vector<GRBVar> refMatchVars, tarMatchVars;
	vector<GRBVar> refLabelVars, tarLabelVars;
	vector<GRBVar> refCycleVars, tarCycleVars;
	map<pair<int, int>, GRBVar> refRAVars, tarRAVars;
	map<pair<int, int>, GRBVar> refPAVars, tarPAVars;
	map<pair<int, int>, GRBVar> pairVars;
	GRBLinExpr cycleSum, indelSum;

	// read files
	readGenome(argv[1], argv[2], ref, tar, refFamilySize, tarFamilySize, maxFamily, refTeloIdx, tarTeloIdx);

	// ILP
	GRBInit(env, model, out_dir);
	matchConstr(ref, tar, refFamilySize, tarFamilySize, maxFamily, model, refMatchVars, tarMatchVars);
	pairConstr(ref, tar, refFamilySize, tarFamilySize, model, pairVars, refMatchVars, tarMatchVars);
	adjConstr(ref, tar, refTeloIdx, tarTeloIdx, model, refRAVars, tarRAVars, refPAVars, tarPAVars);
	cycleNodeConstr(ref, tar, model, refLabelVars, tarLabelVars, refCycleVars, tarCycleVars);
	cycleEdgeConstr(ref, tar, model, refMatchVars, tarMatchVars, pairVars, refRAVars, tarRAVars, refPAVars, tarPAVars, refLabelVars, tarLabelVars);
	objective(model, refCycleVars, tarCycleVars, refMatchVars, tarMatchVars, cycleSum, indelSum);

	// output
	outputJoins(out_dir+"/joins.txt", ref, tar, refPAVars, tarPAVars);
	debug(ref, tar, refTeloIdx, tarTeloIdx, refMatchVars, tarMatchVars, pairVars, refPAVars, tarPAVars, refLabelVars, tarLabelVars, cycleSum, indelSum);

	logFile.close();
	return 0;
}
