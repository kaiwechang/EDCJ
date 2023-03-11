#include "utils.cpp"

vector<Marker> tar;
vector<vector<string>> scaffolds;

map<int, int> joins;
map<string, int> visited;
map<string, Telos> contig2telos;

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
		if (joins.count(tp.lhs) == 0 && joins.count(tp.rhs) == 0 ||
			p.second != 0)
			continue;
		vector<string> curScaffold;

		// traverse scaffolds
		string curContig = p.first;
		int curTelo = joins.count(tp.lhs) ? tp.rhs : tp.lhs;
		while (p.second == 0) {
			int nextTelo = joins[curTelo];
			string nextContig = tar[abs(nextTelo)-1].contig;

			Telos ntp = contig2telos[nextContig];
			if (nextTelo == ntp.lhs) {
				visited[nextContig] = 1;
				curTelo = ntp.rhs;
			} else {
				visited[nextContig] = -1;
				curTelo = ntp.lhs;
			}	curScaffold.push_back(curContig);
			curContig = nextContig;
		}	scaffolds.push_back(curScaffold);
	}

	// add reduced contigs
	fin.open(out_dir+"/removed_spd1.txt");
    while (fin >> contig) {
		scaffolds.push_back(vector<string>(1, contig));
	}	fin.close();

	// merge contigs
	fin.open(out_dir+"/tar_merge.txt");
    while (fin >> contig) {
	
	}	fin.close();

	// cut cycles

	// output final scaffolds
	ofstream fout(out_dir+"/scaffolds.txt");
	for (int i = 0; i < scaffolds.size(); i++) {
		fout << format("> Scaffold_{}\n", i+1);
		for (string& s: scaffolds[i])
			fout << format("{} {}\n", s, visited[s] > -1 ? 0 : 1);
		fout << "\n";
	}	fout.close();
	return 0;
}
