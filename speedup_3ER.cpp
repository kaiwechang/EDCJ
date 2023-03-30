#include "utils.h"
#include <utility>
#include <deque>
using std::pair;
using std::deque;
using std::make_pair;

void readGenome(string refFile, string tarFile, auto& ref, auto& tar, auto& refFamilySize, auto& tarFamilySize, int& maxFamily, auto& refTelos, auto& tarTelos, auto& refOrder, auto& tarOrder) {
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

	//  store telos
	auto storeTelos = [](auto& genome, auto& contig2telos) {
		for (int i = 0; i < genome.size(); i++) {
			if (i == 0 || genome[i-1].contig != genome[i].contig)
				contig2telos[genome[i].contig].lhs = -1 * sign(genome[i].family) * genome[i].id;
			if (i == genome.size()-1 || genome[i].contig != genome[i+1].contig)
				contig2telos[genome[i].contig].rhs =      sign(genome[i].family) * genome[i].id;
		}
	};
	storeTelos(ref, refTelos);
	storeTelos(tar, tarTelos);

	logging("refTelos:\n");
	for (auto& p: refTelos)
		logging("  {}, l: ({}, {}), r: ({}, {})\n", p.first, p.second.lhs, p.second.ljs, p.second.rhs, p.second.rjs);
	logging("tarTelos:\n");
	for (auto& p: tarTelos)
		logging("  {}, l: ({}, {}), r: ({}, {})\n", p.first, p.second.lhs, p.second.ljs, p.second.rhs, p.second.rjs);

	// contig order
	for (int i = 0; i < ref.size()-1; i++) {
		if (i == 0)
			refOrder.push_back(ref[i].contig);
		if (ref[i].contig != ref[i+1].contig)
			refOrder.push_back(ref[i+1].contig);
	}
	for (int i = 0; i < tar.size()-1; i++) {
		if (i == 0)
			tarOrder.push_back(tar[i].contig);
		if (tar[i].contig != tar[i+1].contig)
			tarOrder.push_back(tar[i+1].contig);
	}
}
void speedup(auto& ref, auto& tar, auto& refFamilySize, auto& tarFamilySize, auto& refTelos, auto& tarTelos, bool extended) {
	auto concat = [&](auto& solid, int i, int f1, int f2, int hs1, int hs2, int& js1, int& js2) {
		if (extended &&
			// at least one singleton
			!(refFamilySize[solid[i  ].absFamily] == 1 && tarFamilySize[solid[i  ].absFamily] == 1 ||
			  refFamilySize[solid[i+1].absFamily] == 1 && tarFamilySize[solid[i+1].absFamily] == 1))
			return;
		if (!extended &&
			// all singleton
			!(refFamilySize[solid[i  ].absFamily] == 1 && tarFamilySize[solid[i  ].absFamily] == 1 &&
			  refFamilySize[solid[i+1].absFamily] == 1 && tarFamilySize[solid[i+1].absFamily] == 1))
			return;
		if (// not joined yet
			//js1 == 0 && js2 == 0 &&
			// solid ?
			solid[i].contig == solid[i+1].contig &&
			// adjacency ?
			(solid[i].family ==  f1 && solid[i+1].family ==  f2 ||
			 solid[i].family == -f2 && solid[i+1].family == -f1)) {
			js1 = hs2, js2 = hs1;
		}
	};
	auto iterTelos = [&](auto& g1, auto& g2, auto& contig2telos) {
		for (int i = 0; i < g1.size()-1; i++)
			for (auto it1 = tarTelos.begin(); std::next(it1, 1) != tarTelos.end(); it1++)
				for (auto it2 = std::next(it1, 1); it2 != tarTelos.end(); it2++) {
					auto& [c1, tp1] = *it1;
					auto& [c2, tp2] = *it2;
					if (c1 == c2)
						continue;
					//fmt::print("ref[{}]: ({}, {})\n", i, cp1.first, cp2.first);
					int f1l = g2[idx(tp1.lhs)].family;
					int f1r = g2[idx(tp1.rhs)].family;
					int f2l = g2[idx(tp2.lhs)].family;
					int f2r = g2[idx(tp2.rhs)].family;
					concat(g1, i,  f1l, -f2l, tp1.lhs, tp2.lhs, tp1.ljs, tp2.ljs);
					concat(g1, i, -f1l, -f2r, tp1.lhs, tp2.rhs, tp1.ljs, tp2.rjs);
					concat(g1, i,  f1r,  f2l, tp1.rhs, tp2.lhs, tp1.rjs, tp2.ljs);
					concat(g1, i,  f1r, -f2r, tp1.rhs, tp2.rhs, tp1.rjs, tp2.rjs);
				}
	};
	iterTelos(ref, tar, tarTelos);
	iterTelos(tar, ref, refTelos);

	logging("spd3:\n");
	logging("refTelos:\n");
	for (auto& p: refTelos)
		logging("  {}, l: ({}, {}), r: ({}, {})\n", p.first, p.second.lhs, p.second.ljs, p.second.rhs, p.second.rjs);
	logging("tarTelos:\n");
	for (auto& p: tarTelos)
		logging("  {}, l: ({}, {}), r: ({}, {})\n", p.first, p.second.lhs, p.second.ljs, p.second.rhs, p.second.rjs);
}
auto joinGenome(auto& genome, auto& contig2telos, auto& draft, auto& genomeOrder) {
	map<string, bool> visited;
	auto travse = [&](auto& curScaffold, Telos tp, bool right) {
		int hs = right ? tp.rhs : tp.lhs ;
		int js = right ? tp.rjs : tp.ljs ;
		while (js != 0) {
			string nc = genome[idx(js)].contig;
			if (visited[nc])
				break;

			Telos ntp = contig2telos[genome[idx(js)].contig];
			if (right)
				curScaffold.push_back(make_pair(nc, js == ntp.lhs ? 0 : 1));
			else
				curScaffold.push_front(make_pair(nc, js == ntp.lhs ? 1 : 0));
			visited[nc] = true;

			hs = js == ntp.lhs ? ntp.rhs : ntp.lhs ;
			js = js == ntp.lhs ? ntp.rjs : ntp.ljs ;
		}
	};
	for (Marker& m: genome)
		visited[m.contig] = false;
	for (auto& [c, v]: visited) {
		deque<pair<string, int>> curScaffold;
		if (v)
			continue;
		Telos tp = contig2telos[c];

		curScaffold.push_back(make_pair(c, 0));
		v = true;

		travse(curScaffold, tp, true);
		travse(curScaffold, tp, false);

		draft.push_back(curScaffold);
	}
	// make new genome
	int uidx = 1;
	vector<Marker> newGenome;
	for (auto& sv: draft) {
		string leader = sv[0].first;
		for (auto& p: sv) {
			Telos tp = contig2telos[p.first];
			int s = idx(p.second ? tp.rhs : tp.lhs);
			int e = idx(p.second ? tp.lhs : tp.rhs);
			//fmt::print("contig: {}, s: {}, e: {}\n", p.first, s, e);
			if (p.second == 0) {
				for (int i = s; i <= e; i++)
					newGenome.push_back(Marker(uidx++, genome[i].family, leader));
			} else {
				for (int i = s; i >= e; i--)
					newGenome.push_back(Marker(uidx++, -1 * genome[i].family, leader));
			}
		}
	}
	//fmt::print("g: {}, ng: {}\n", genome.size(), newGenome.size());
	logging("{} contigs in draft\n", draft.size());
	return newGenome;
}
void outputMergeContigs(string filename, auto& draft) {
	ofstream fout(filename);
	for (auto& sv: draft)
		if (sv.size() > 1) {
			fout << format("{}\n", sv[0].first);
			for (auto& p: sv)
				fout << format("{} {}\n", p.first, p.second);
			fout <<  format("end\n");
		}
}
void outputNewGenome(string refFile, string tarFile, auto& ref, auto& tar, auto& reorder) {
	ofstream fout(refFile);
	for (Marker& m: ref) {
		fout << format("{} {} {} 1\n", m.id, m.family, m.contig);
	}	fout.close();
	fout.open(tarFile);
	for (Marker& m: tar) {
		fout << format("{} {} {} 1\n", m.id, m.family, m.contig);
	}	fout.close();
}
int main(int argc, char *argv[]) {
	if (argc < 4) {
		fmt::print("[error] Usage:\n>>> speedup_1 <ref genome> <tar genome> <output_dir>\n");
		return 0;
	}	string out_dir(argv[3]);
	logFile.open(out_dir+"/speedup_3.log");

	int maxFamily = 0;
	vector<Marker> ref, tar;
	vector<string> refOrder, tarOrder;
	vector<int> refFamilySize, tarFamilySize, reorder;
	vector<deque<pair<string, int>>> refDraft, tarDraft;
	map<string, Telos> refTelos, tarTelos;

	// read ref/tar
	readGenome(argv[1], argv[2], ref, tar, refFamilySize, tarFamilySize, maxFamily,
								 refTelos, tarTelos, refOrder, tarOrder);

	// speedup 1
	bool extended = argc == 4 ? false : true ;
	speedup(ref, tar, refFamilySize, tarFamilySize, refTelos, tarTelos, extended);
	ref = joinGenome(ref, refTelos, refDraft, refOrder);
	tar = joinGenome(tar, tarTelos, tarDraft, tarOrder);

	// output new ref/tar & mergeContigs
	outputNewGenome(out_dir+"/ref_spd3.all", out_dir+"/tar_spd3.all", ref, tar, reorder);
	outputMergeContigs(out_dir+"/ref_merge.txt", refDraft);
	outputMergeContigs(out_dir+"/tar_merge.txt", tarDraft);

	logFile.close();
	return 0;
}
