#include <vector>
#include <fstream>
#include <sstream>
#include <fmt/core.h>

using fmt::print;
using fmt::format;

using std::vector;
using std::string;

using std::getline;
using std::ifstream;
using std::ofstream;
using std::stringstream;

void readFNA(string inFile, int& sizeMb, int& numGC, int& numContigs) {
	int uid = 0;
	string line, token;
	ifstream fin(inFile);

	while (getline(fin, line)) {
		if (line[0] == '>') {
			// print("{}\n", line);
			numContigs++;
		} else {
			for (auto& c: line)
				switch(c) {
					case 'G': case 'C':
						numGC++;
						sizeMb++;
						break;
					case 'A': case 'T': case 'U':
						sizeMb++;
						break;
					default:
						//print("Unknown character: {}\n", c);
						//exit(1);
						sizeMb++;
				}
				//print("{}", line[i]);
				continue;
		}
	}	fin.close();
}
int main(int argc, char* argv[]) {
	if (argc < 2) {
		print("[error] Usage:\n>>> count_FNA <FNA file>\n");
		return 0;
	}
	int sizeMb, numGC, numContigs;
	string inFile(argv[1]);
	readFNA(inFile, sizeMb, numGC, numContigs);
	print("  GC content  : {:.2f}%\n", (double)numGC/sizeMb*100);
	print("  genome size : {:.2f}Mb\n", (double)sizeMb/1000000);
	print("  # of contigs: {}\n", numContigs);
	print("\n");

	return 0;
}
