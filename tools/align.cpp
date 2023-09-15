#include <set>
#include <map>
#include <vector>
#include <utility>

#include <sstream>
#include <fstream>
#include <fmt/core.h>

using fmt::print;
using fmt::format;

using std::set;
using std::map;
using std::pair;
using std::vector;
using std::string;

using std::getline;
using std::ifstream;
using std::ofstream;
using std::stringstream;

#define s first
#define e second

struct Ans {
	string tar;
	bool direct;
};
struct Info {
	string ref;
	pair<int, int> pos;
	double cov;
	bool direct;
};
bool overlap(auto& a, auto& b) {
	return !(a.s < b.s && a.e < b.s || a.s > b.e && a.e > b.e);
}
void readCoords(string coordsFile, auto& alignInfo, double idyThres, double covThres) {
	int uid = 0;
	string line, tmp;
	ifstream fin(coordsFile);

	string ref, tar;
	double identity, refCov, tarCov;
	int refStart, refEnd, tarStart, tarEnd, refLen, tarLen;

	// skip header
	while (getline(fin, line))
		if (line[0] == '=')
			break;
	// parse content
	while (getline(fin, line)) {
		stringstream ss(line);
		ss >> refStart >> refEnd >> tmp >> tarStart >> tarEnd >> tmp >> tmp >> tmp >> tmp
		   >> identity >> tmp >> refLen >> tarLen >> tmp >> refCov >> tarCov >> tmp >> ref >> tar;
		// threshold
		if (identity < idyThres || tarCov < covThres)
			continue;
		// real position
		bool direct = tarEnd-tarStart > 0;
		int left  = direct ? tarStart-1 : tarLen-tarStart;
		int right = direct ? tarLen-tarEnd : tarEnd-1;
		pair<int, int> pos = {refStart - left, refEnd + right};

		// store info
		if (alignInfo.count(tar) != 0) {
			alignInfo[tar].push_back({ref, pos, tarCov, direct});
		} else {
			alignInfo[tar] = vector<Info>(1, {ref, pos, tarCov, direct});
		}
	}	fin.close();
	/*for (auto& [tar, vec]: alignInfo) {
		print("{}:\n", tar);
		for (Info& info: vec)
			print("    {}, ({}, {}), {}, {}\n", info.ref, info.pos.s, info.pos.e, info.cov, info.direct);
	}*/
}
void intervalSchedule(auto& alignInfo, auto& scaffs, auto& remain) {
	// filter intervals
	for (auto& [tar, vec]: alignInfo) {
		// add single intervals only
		if (vec.size() == 1) {
			Info& info = vec[0];
			scaffs[info.ref][pair(info.pos.s, info.pos.e)] = {tar, info.direct};
		} else {
			// skip group intervals for simplicity
			//print("group intervals: {}\n", tar);
			/*auto selectMaxCovInfo = [](auto& vec) {
				double maxCov = vec[0].cov;
				Info& maxInfo = vec[0];
				for (Info& info: vec)
					if (info.cov > maxCov)
						maxCov = info.cov, maxInfo = info;
				return maxInfo;
			};
			Info info = selectMaxCovInfo(vec);
			scaffs[info.ref][pair(info.pos.s, info.pos.e)] = {tar, info.direct};*/
		}
	}

	// deal with overlap
	auto prevInterval = pair(-2, -1);
	for (auto& [ref, intervals]: scaffs)
		for (auto it = intervals.begin(); it != intervals.end(); ) {
			auto& [pos, ans] = *it;
			if (overlap(prevInterval, pos)) {
				intervals.erase(it++);
			} else {
				it++;
				remain.insert(ans.tar);
				prevInterval = pair(pos.s, pos.e);
			}
		}

	// print debug
	for (auto& [ref, intervals]: scaffs) {
		print("{}:\n", ref);
		for (auto& [pos, ans]: intervals)
			print("    ({}, {}): {}, {}\n", pos.s, pos.e, ans.tar, ans.direct ? 0 : 1);
	}
	/*print("\nremain:\n");
	for (auto& tar: remain)
		print("{}\n", tar);*/
	print("tar contig num: {}\n", remain.size());
}
void outputAnswer(string file, auto& scaffs) {
	int uid = 1;
	ofstream fout(file);
	for (auto& [ref, intervals]: scaffs) {
		fout << format("> Scaffold_{}\n", uid++);
		for (auto& [pos, ans]: intervals)
			fout << format("{} {}\n", ans.tar, ans.direct ? 0 : 1);
	}	fout.close();
}
void outputFNA(string inFile, string outFile, auto& remain) {
	ifstream fin(inFile);
	ofstream fout(outFile);

	bool flag;
	string line, contig;
	while (getline(fin, line)) {
		if (line[0] == '>') {
			stringstream ss(string(line.begin()+1, line.end()));
			ss >> contig;

			flag = remain.count(contig) != 0;
		}
		if (flag)
			fout << line << "\n";
	}
	fin.close();
	fout.close();
}
int main(int argc, char* argv[]) {
	if (argc < 6) {
		print("[error] Usage:\n>>> align <comp FNA> <draft FNA> <identity threshold> <coverage threshold> <output dir>\n");
		return 1;
	}
	string refFNA(argv[1]);		// complete fna
	string tarFNA(argv[2]);		// draft fna
	string outDir(argv[5]);
	double idyThres = std::stod(argv[3]);
	double covThres = std::stod(argv[4]);

	// MUMmer stuff
	string command = format("nucmer {} {} --delta={}/delta", refFNA, tarFNA, outDir);
	if (system(command.c_str()) != 0) {
		print("[error] Cannot execute nucmer commands!!\n");
		exit(1);
	}
	command = format("show-coords {}/delta -cql > {}/coords", outDir, outDir);
	if (system(command.c_str()) != 0) {
		print("[error] Cannot execute show-coords commands!!\n");
		exit(1);
	}

	// generate answer with alignment result
	set<string> remain;
	map<string, vector<Info>> alignInfo;
	map<string, map<pair<int, int>, Ans, 
		// sort with end pos
		decltype([](auto& a, auto& b) {
			return a.e < b.e;
		})>> scaffs;

	readCoords(outDir+"/coords", alignInfo, idyThres, covThres);
	intervalSchedule(alignInfo, scaffs, remain);

	outputAnswer(outDir+"/answerToAll", scaffs);
	outputFNA(tarFNA, outDir+"/contigMerged.v3.randOrd", remain);

	return 0;
}
