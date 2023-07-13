#include <vector>
#include <utility>

#include <random>
#include <algorithm>

#include <sstream>
#include <fstream>
#include <fmt/core.h>

using fmt::print;
using fmt::format;

using std::pair;
using std::vector;
using std::string;

using std::getline;
using std::ifstream;
using std::ofstream;
using std::stringstream;

// random pick k numbers from range [0, n)
void rand_pick(auto& picked, int k, int n) {
	vector<int> nums(n, 0);
	for (int i = 0; i < n; i++)
		nums[i] = i;

	std::random_device rd;
    std::mt19937 gen(rd());
	shuffle(nums.begin(), nums.end(), gen);

	for (int i = 0; i < k; i++)
		picked.push_back(nums[i]);
	std::sort(picked.begin(), picked.end());
}
void readFNA(string cmpFNA, auto& seqLines, int& scaff_num) {
	scaff_num = 0;
	string line, token;
	ifstream fin(cmpFNA);

	while (getline(fin, line)) {
		if (line[0] == '>') {
			scaff_num++;
		} else {
			seqLines.push_back(pair(scaff_num, line));
		}
	}	fin.close();
}
void outputFNA(string file, auto& seqLines, auto& cutIdx, auto& ansIdx) {
	int uid = 0, cidx = 0;

	ofstream fout(file);
	fout << format(">contig_{}\n", ++uid);
	for (int i = 0; i < seqLines.size(); i++) {
		if (cidx >= cutIdx.size() || i < cutIdx[cidx]) {
			if (i > 0 && seqLines[i-1].first != seqLines[i].first) {
				fout << format("\n>contig_{} original breakpoint here\n", ++uid);
				ansIdx.push_back(uid);
			}	fout << seqLines[i].second << "\n";
		} else {
			fout << format("\n>contig_{}\n", ++uid);
			cidx++;
		}
	}	fout.close();
}
void outputAnswer(string file, auto& contig_num, auto& ansIdx) {
	int uid = 0, aid = 0;
	ofstream fout(file);
	fout << format("> Scaffold_1\n");
	for (int i = 1; i <= contig_num; i++) {
		if (uid < ansIdx.size() && i >= ansIdx[uid]) {
			fout << format("\n> Scaffold_{}\n", ++uid+1);
		}	fout << format("contig_{} 0\n", i);
	}	fout.close();
}
int main(int argc, char* argv[]) {
	if (argc < 4) {
		print("[error] Usage:\n>>> cut_fna <comp FNA> <contig num> <output dir>\n");
		return 0;
	}
	string cmpFNA(argv[1]);
	string outDir(argv[3]);
	int contig_num(std::stoi(argv[2]));

	int scaff_num;
	vector<int> cutIdx, ansIdx;
	vector<pair<int, string>> seqLines;

	readFNA(cmpFNA, seqLines, scaff_num);
	rand_pick(cutIdx, contig_num-scaff_num, seqLines.size());
	outputFNA(format("{}/draft_{}.fna", outDir, contig_num), seqLines, cutIdx, ansIdx);
	outputAnswer(outDir+"/answerToAll", contig_num, ansIdx);

	return 0;
}
