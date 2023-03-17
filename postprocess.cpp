#include "utils.h"
#include <utility>
using std::pair;
using std::make_pair;

void readGenome(string refFile, string tarFile, auto& ref, auto& tar) {
	string contig;
	int id, family, tmp;
	ifstream fin(refFile);
	while (fin >> id >> family >> contig >> tmp) {
		ref.push_back(Marker(id, family, contig));
	}	fin.close();
	fin.open(tarFile);
	while (fin >> id >> family >> contig >> tmp) {
		tar.push_back(Marker(id, family, contig));
	}	fin.close();
}
void readJoins(string filename, auto& tar, auto& joins, auto& contig2telos) {
	int t1, t2;
	string contig;
	ifstream fin(filename);
    while (fin >> contig >> t1 >> t2) {
		joins[t1] = t2, joins[t2] = t1;
		//Marker m1 = tar[abs(t1)-1], m2 = tar[abs(t2)-1];
		//t1 * m1.family > 0 ? contig2telos[m1.contig].rhs = t1 : contig2telos[m1.contig].lhs = t1;
		//t2 * m2.family > 0 ? contig2telos[m2.contig].rhs = t2 : contig2telos[m2.contig].lhs = t2;
	}	fin.close();

	// make contig2telos
	for (int i = 0; i < tar.size(); i++) {
		int s1 = tar[i  ].family > 0 ? 1 : -1 ;
		int s2 = tar[i+1].family > 0 ? 1 : -1 ;
		if (i == 0)
			contig2telos[tar[i].contig].lhs = -s1 * tar[i].id;
		if (i == tar.size()-1) {
			contig2telos[tar[i].contig].rhs =  s1 * tar[i].id;
			break;
		}
		if (tar[i].contig != tar[i+1].contig) {
			contig2telos[tar[i  ].contig].rhs =  s1 * tar[i  ].id;
			contig2telos[tar[i+1].contig].lhs = -s2 * tar[i+1].id;
		}
	}
	for (auto& p: contig2telos)
		logging("{}: ({}, {})\n", p.first, p.second.lhs, p.second.rhs);
}
void readMergeContigs(string filename, auto& tar, auto& mergeContig) {
	int tmp;
	string key, contig;
	ifstream fin(filename);
    while (fin >> key) {
		vector<pair<string, int>> seg;
		while (fin >> contig) {
			if (contig == "end")
				break;
			fin >> tmp;
			seg.push_back(make_pair(contig, tmp));
		}	mergeContig[key] = seg;
	}	fin.close();

	for (auto& p: mergeContig) {
		logging("{}:\n", p.first);
		for (auto& cp: p.second)
			logging("  {} {}\n", cp.first, cp.second);
	}

	// 2023_0317: merge with joins
	/*for (auto& p: mergeContig) {
		fmt::print("joins {}:\n", p.first);
		for (int i = 0; i < p.second.size()-1; i++) {
			auto& tp1 = contig2telos[p.second[i  ].first];
			auto& tp2 = contig2telos[p.second[i+1].first];
			int t1 = p.second[i  ].second == 0 ? tp1.rhs: tp1.lhs;
			int t2 = p.second[i+1].second == 0 ? tp2.lhs: tp2.rhs;
			fmt::print("  {}, {}\n", t1, t2);
			joins[t1] = t2, joins[t2] = t1;
		}
	}*/
}
void scaffolding(auto& tar, auto& joins, auto& contig2telos, auto& scaffolds) {
	map<string, bool> visited;
	for (Marker& m: tar)
		visited[m.contig] = false;

	for (auto& p: visited) {
		Telos tp = contig2telos[p.first];
		vector<pair<string, int>> curScaffold;
		if (joins.count(tp.lhs) == 0 && joins.count(tp.rhs) == 0 || p.second)
			continue;
		// traverse scaffolds
		string curContig = p.first;
		int curTelo = joins.count(tp.rhs) ? tp.rhs : tp.lhs;
		while (!p.second) {
			int nextTelo = joins[curTelo];
			string nextContig = tar[abs(nextTelo)-1].contig;

			Telos ctp = contig2telos[curContig];
			Telos ntp = contig2telos[nextContig];
			curScaffold.push_back(make_pair(curContig, curTelo == ctp.lhs ? 1 : 0));
			//curScaffold.push_back(make_pair(nextContig, nextTelo == ntp.lhs ? 0 : 1));
			curTelo = nextTelo == ntp.lhs ? ntp.rhs : ntp.lhs;

			visited[nextContig] = true;
			curContig = nextContig;
		}	scaffolds.push_back(curScaffold);
	}
}
void addReducedContigs(string filename, auto& scaffolds) {
	string contig;
	int id, family, tmp;
	ifstream fin(filename);
    while (fin >> contig) {
		scaffolds.push_back(vector<pair<string, int>>(1, make_pair(contig, 0)));
	}	fin.close();
}
void cutCycle(auto& ref, auto& tar, auto& contig2telos, auto& scaffolds) {
	// store ref adjacencies (in "signed extremities" format)
	set<pair<int, int>> refAdj;
	for (int i = 0; i < ref.size()-1; i++)
		if (ref[i].contig == ref[i+1].contig)
			refAdj.insert(ref[i].absFamily < ref[i+1].absFamily ?
				make_pair( ref[i  ].family, -ref[i+1].family) :
				make_pair(-ref[i+1].family,  ref[i  ].family));
	for (auto& p: refAdj)
		fmt::print("adj: ({}, {})\n", p.first, p.second);

	// collect candidate cuts
	for (auto& scaff: scaffolds) {
		if (scaff.size() < 2)
			continue;
		vector<pair<int, int>> candicates;

		fmt::print("scaff\n");
		for (int i = 0; i < scaff.size(); i++) {
			int idx1 = i, idx2 = i == scaff.size()-1 ? 0 : i+1 ;
			Telos tp1 = contig2telos[scaff[idx1].first];
			Telos tp2 = contig2telos[scaff[idx2].first];
			int t1 = scaff[idx1].second == 0 ? tp1.rhs : tp1.lhs ;
			int t2 = scaff[idx2].second == 0 ? tp2.lhs : tp2.rhs ;
			int s1 = t1 > 0 ? 1 : -1, s2 = t2 > 0 ? 1 : -1 ;
			Marker& m1 = tar[abs(t1)-1], m2 = tar[abs(t2)-1];
			auto p = m1.absFamily < m2.absFamily ?
				make_pair(s1 * m1.absFamily, s2 * m2.absFamily): 
				make_pair(s2 * m2.absFamily, s1 * m1.absFamily);

			if (refAdj.count(p) == 0) {
				fmt::print("  ({}, {}): ({}, {})\n", p.first, p.second, idx1, idx2);
				candicates.push_back(make_pair(idx1, idx2));
			}
		}

		// determine which join to cut
		int cut = 2;
		// map<absFamily, map<contig, pair<indexSum, familyNum>>>
		map<int, map<string, pair<int, int>>> score;
		for (int i = 0; i < ref.size(); i++) {
			auto& p = score[ref[i].absFamily][ref[i].contig];
			p.first += i, p.second++;
		}
		for (Marker& m: tar) {

		}
		for (auto& fp: score) {
			fmt::print("family {}:\n", fp.first);
			for (auto& cp: fp.second)
				fmt::print("    {}: {}/{}\n", cp.first, cp.second.first, cp.second.second);
		}

		// cut scaffold
		vector<pair<string, int>> newScaffold;
		for (int i = candicates[cut].second; i < scaff.size(); i++)
			newScaffold.push_back(scaff[i]);
		for (int i = 0; i < candicates[cut].second; i++)
			newScaffold.push_back(scaff[i]);
		scaff = newScaffold;
	}
}
void outputScaffold(string filename, auto& scaffolds, auto& mergeContig) {
	auto mergeContigSign = [&](string target, auto merged) {
		// find target(leader) contig in the merged contigs
		// and then return its sign
		for (auto& mp: merged)
			if (mp.first == target)
				return mp.second;
		return -1;
	};

	ofstream fout(filename);
	for (int i = 0; i < scaffolds.size(); i++) {
		fout << format("> Scaffold_{}\n", i+1);
		for (auto& sp: scaffolds[i]) {
			auto it = mergeContig.find(sp.first);
			if (it != mergeContig.end()) {
				int sign = mergeContigSign(sp.first, it->second);
				logging("{}: {}, {}\n", sp.first, sp.second, sign);
				// leader contig sign XOR its sign in the merged contigs
				if (sp.second + sign == 1) {
					for (int i = it->second.size()-1; i > -1; i--)
						fout << format("{} {}\n", it->second[i].first, 1-it->second[i].second);
				} else {
					for (auto& mp: it->second)
						fout << format("{} {}\n", mp.first, mp.second);
				}
			} else
				fout << format("{} {}\n", sp.first, sp.second);
		}	fout << "\n";
	}	fout.close();
}
int main(int argc, char *argv[]) {
	if (argc < 4) {
		fmt::print("[error] Usage:\n>>> postprocess <tar genome> <joins.txt> <output_dir>\n");
		return 0;
	}	string out_dir(argv[3]);
	logFile.open(out_dir+"/postprocess.log");

	vector<Marker> ref, tar;
	vector<vector<pair<string, int>>> scaffolds;

	map<int, int> joins;
	map<string, Telos> contig2telos;
	map<string, vector<pair<string, int>>> mergeContig;

	// read files
	readGenome(argv[1], argv[2], ref, tar);
	readJoins(out_dir+"/joins.txt", tar, joins, contig2telos);
	readMergeContigs(out_dir+"/tar_merge.txt", tar, mergeContig);

	// scaffolding
	scaffolding(tar, joins, contig2telos, scaffolds);
	addReducedContigs(out_dir+"/removed_spd1.txt", scaffolds);

	// cut cycles
	//cutCycle(ref, tar, contig2telos, scaffolds);

	// output final scaffold
	outputScaffold(out_dir+"/scaffolds.txt", scaffolds, mergeContig);

	logFile.close();
	return 0;
}
