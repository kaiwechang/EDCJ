#include "utils.cpp"
#include <fstream>
#include <map>

using std::map;

vector<Marker> tar;
vector<vector<string>> scaffolds;

map<int, int> joins;
map<string, int> contigSize;
map<string, Telos> contig2telos;

int main(int argc, char *argv[])
{
	if (argc < 4) {
		print("[error] Usage:\n>>> postprocess <tar genome> <joins.txt> <output_dir>\n");
		return 0;
	}
	string out_dir(argc < 4 ? "output" : argv[3]);

	string contig;
	int id, family, tmp;
	std::ifstream fin(argv[1]);
	while (fin >> id >> family >> contig >> tmp) {
		tar.push_back(Marker(id, family, contig));
		contigSize[contig]++;
	}	fin.close();

	// joins.txt
	int t1, t2;
	fin.open(argv[2]);
    while (fin >> contig >> t1 >> t2) {
		joins[t1] = t2, joins[t2] = t1;
		Marker m1 = tar[abs(t1)-1], m2 = tar[abs(t2)-1];
		t1 * m1.family > 0 ? contig2telos[m1.contig].rhs = t1 : contig2telos[m1.contig].lhs = t1;
		t2 * m2.family > 0 ? contig2telos[m2.contig].rhs = t2 : contig2telos[m2.contig].lhs = t2;
	}	fin.close();

	/*for (auto& p: contigSize)
		print("{} {}\n", p.first, p.second);*/
	for (auto& p: contig2telos)
		print("{}: ({}, {})\n", p.first, p.second.lhs, p.second.rhs);

	fin.open(out_dir+"/duplicated.txt");
    while (fin >> contig) {
		contigSize[contig]--;
	}	fin.close();

	/*print("remove duplicated\n");
	for (auto& p: contigSize)
		print("{} {}\n", p.first, p.second);*/

	// add removed contigs into scaffold first
	// if (contigSize[contig] == 0)
	//	scaffolds.push_back(vector(contig))

	int uid = 1;
	// scaffolding
	for (auto& p: contigSize) {
		if (p.second == 0)
			continue;
		print(">scaffold_{}\n", uid++);
		// traverse scaffolds
		string curContig = p.first, nextContig;
		int curTelo = contig2telos[curContig].rhs, nextTelo;
		while (p.second > 0) {
			nextTelo = joins[curTelo];
			nextContig = tar[abs(nextTelo)-1].contig;

			Telos ctp = contig2telos[curContig];
			Telos ntp = contig2telos[nextContig];
			if (curTelo == ctp.rhs) {
				print("{} 0\n", curContig);
			} else {
				print("{} 1\n", curContig);
			}	contigSize[nextContig] = 0;
			curTelo = nextTelo == ntp.lhs ? ntp.rhs : ntp.lhs ;
			curContig = nextContig;
		}	print("\n");
	}

	// tarContigMerge.txt
	// 
	std::ofstream fout(out_dir+"/myScaffold.txt");
	for (int i = 0; i < 10; i++) {
		/*fout << ">scaffold_" << scaffoldID++ << endl;
		for (auto contig : currScaffold)
		{
			fout << contig.first << " " << contig.second << endl;
		}
		fout << endl;*/
	}	fout.close();
	return 0;
}
