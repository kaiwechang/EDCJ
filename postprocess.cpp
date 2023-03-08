#include "utils.cpp"
#include <fstream>
#include <map>

using std::map;

map<int, string> telo2contig;
map<string, int> contigSize;

int main(int argc, char *argv[])
{
	string target_dir(argc < 4 ? "output" : argv[3]);

	ifstream fin;
	fin.open(argv[1]);

	string contig;
	int id, family, tmp;
	while (fin >> id >> family >> contig >> tmp) {

	}	fin.close();

	fin.open(target_dir+"/duplicated.txt");
    while (fin >> contig) {
		contigSize[contig]--;
	}	fin.close();


	// joins.txt
	// tarContigMerge.txt
	// 
	ofstream fout(target_dir+"/myScaffold.txt");
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
