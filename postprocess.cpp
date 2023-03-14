#include "utils.h"
#include <utility>
using std::pair;

void readGenome(string filename, auto& tar, auto& visited) {
	string contig;
	int id, family, tmp;
	ifstream fin(filename);
	while (fin >> id >> family >> contig >> tmp) {
		tar.push_back(Marker(id, family, contig));
		visited[contig] = 0;
	}	fin.close();
}
void readJoins(string filename, auto& tar, auto& joins, auto& contig2telos) {
	int t1, t2;
	string contig;
	ifstream fin(filename);
    while (fin >> contig >> t1 >> t2) {
		joins[t1] = t2, joins[t2] = t1;
		Marker m1 = tar[abs(t1)-1], m2 = tar[abs(t2)-1];
		t1 * m1.family > 0 ? contig2telos[m1.contig].rhs = t1 : contig2telos[m1.contig].lhs = t1;
		t2 * m2.family > 0 ? contig2telos[m2.contig].rhs = t2 : contig2telos[m2.contig].lhs = t2;
	}	fin.close();

	for (auto& p: contig2telos)
		logging("{}: ({}, {})\n", p.first, p.second.lhs, p.second.rhs);
}
auto scaffolding(auto& tar, auto& joins, auto& visited, auto& contig2telos) {
	vector<vector<pair<string, int>>> scaffolds;

	for (auto& p: visited) {
		Telos tp = contig2telos[p.first];
		vector<pair<string, int>> curScaffold;
		if (joins.count(tp.lhs) == 0 && joins.count(tp.rhs) == 0 || p.second)
			continue;
		// traverse scaffolds
		string curContig = p.first;
		int curTelo = joins.count(tp.lhs) ? tp.rhs : tp.lhs;
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
	return scaffolds;
}
auto addReducedContigs(string filename, auto scaffolds) {
	string contig;
	int id, family, tmp;
	ifstream fin(filename);
    while (fin >> contig) {
		scaffolds.push_back(vector<pair<string, int>>(1, make_pair(contig, 0)));
	}	fin.close();
	return scaffolds;
}
void readMergeContigs(string filename, auto& mergeContig) {
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
}
int mergeContigSign(string target, auto merged) {
	// find target(leader) contig in the merged contigs
	// and then return its sign
	for (auto& mp: merged)
		if (mp.first == target)
			return mp.second;
	return -1;
}
void outputScaffold(string filename, auto& scaffolds, auto& mergeContig) {
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

	vector<Marker> tar;
	vector<vector<pair<string, int>>> scaffolds;

	map<int, int> joins;
	map<string, bool> visited;
	map<string, Telos> contig2telos;
	map<string, vector<pair<string, int>>> mergeContig;

	// read files
	readGenome(argv[1], tar, visited);
	readJoins(argv[2], tar, joins, contig2telos);
	readMergeContigs(out_dir+"/tar_merge.txt", mergeContig);

	// scaffolding
	scaffolds = scaffolding(tar, joins, visited, contig2telos);
	scaffolds = addReducedContigs(out_dir+"/removed_spd1.txt", scaffolds);

	// cut cycles

	// output final scaffold
	outputScaffold(out_dir+"/scaffolds.txt", scaffolds, mergeContig);

	logFile.close();
	return 0;
}
