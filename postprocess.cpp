#include "utils.cpp"
#include <utility>
using std::pair;

vector<Marker> tar;
vector<vector<pair<string, int>>> scaffolds;

map<int, int> joins;
map<string, bool> visited;
map<string, Telos> contig2telos;
map<string, vector<pair<string, int>>> mergeContig;

int main(int argc, char *argv[])
{
	if (argc < 4) {
		print("[error] Usage:\n>>> postprocess <tar genome> <joins.txt> <output_dir>\n");
		return 0;
	}	string out_dir(argv[3]);

	string contig;
	int id, family, tmp;
	ifstream fin(argv[1]);
	while (fin >> id >> family >> contig >> tmp) {
		tar.push_back(Marker(id, family, contig));
		visited[contig] = 0;
	}	fin.close();

	// read joins
	int t1, t2;
	fin.open(argv[2]);
    while (fin >> contig >> t1 >> t2) {
		joins[t1] = t2, joins[t2] = t1;
		Marker m1 = tar[abs(t1)-1], m2 = tar[abs(t2)-1];
		t1 * m1.family > 0 ? contig2telos[m1.contig].rhs = t1 : contig2telos[m1.contig].lhs = t1;
		t2 * m2.family > 0 ? contig2telos[m2.contig].rhs = t2 : contig2telos[m2.contig].lhs = t2;
	}	fin.close();

	for (auto& p: contig2telos)
		print("{}: ({}, {})\n", p.first, p.second.lhs, p.second.rhs);

	// scaffolding
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

	// add reduced contigs
	fin.open(out_dir+"/removed_spd1.txt");
    while (fin >> contig) {
		scaffolds.push_back(vector<pair<string, int>>(1, make_pair(contig, 0)));
	}	fin.close();

	// merge contigs
	fin.open(out_dir+"/tar_merge.txt");
	string key;
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
		print("{}:\n", p.first);
		for (auto& cp: p.second)
			print("  {} {}\n", cp.first, cp.second);
	}

	// cut cycles

	// output final scaffolds
	ofstream fout(out_dir+"/scaffolds.txt");
	for (int i = 0; i < scaffolds.size(); i++) {
		fout << format("> Scaffold_{}\n", i+1);
		for (auto& sp: scaffolds[i]) {
			auto it = mergeContig.find(sp.first);
			if (it != mergeContig.end()) {
				int sign;
				for (auto& mp: it->second)
					if (mp.first == sp.first)
						sign = mp.second;
				print("{}: {}, {}\n", sp.first, sp.second, sign);
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
	return 0;
}
