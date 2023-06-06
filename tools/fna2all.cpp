#include <map>
#include <vector>
#include <utility>
#include <fstream>
#include <cstdlib>
#include <fmt/core.h>

using fmt::print;
using fmt::format;

using std::vector;
using std::string;

using std::map;
using std::pair;
using std::make_pair;

using std::ifstream;
using std::ofstream;

void readFNA(string refFNA, string tarFNA, auto& refContigs, auto& tarContigs) {
	string token;
	std::ifstream fin(refFNA);
	while (fin >> token) {
		if (token[0] == '>')
			refContigs[string(token.begin()+1, token.end())] = 0;
	}	fin.close();

	fin.open(tarFNA);
	while (fin >> token) {
		if (token[0] == '>')
			tarContigs[string(token.begin()+1, token.end())] = 0;
	}	fin.close();
}
void readBlock(string blockFile, auto& refContigs, auto& tarContigs, auto& reorder) {
	// remove blocks only located on ref/tar
	// and make the new marker mapping, "reorder"
	string token, block, marker, seqId, strand, start, end, length;
	ifstream fin(blockFile);
	while (fin >> token)
		if (token[0] == '-')
			break;

	int uid = 1;
	while (fin >> block >> marker) {
		int family = stoi(string(marker.begin()+1, marker.end()));
		//print("marker: {}\n", family);
		fin >> seqId >> strand >> start >> end >> length;

		bool onRef = false, onTar = false;
		while (fin >> seqId) {
			if (seqId[0] == '-')
				break;
			fin >> strand >> start >> end >> length;
			// check if the marker appears on both ref & tar
			stoi(seqId) <= refContigs.size() ? onRef = true : onTar = true ;
		}	reorder.push_back(onRef && onTar ? uid++ : -1);
	}
}
void readPermu(string permuFile, auto& ref, auto& tar, auto& refContigs, auto& tarContigs) {
	ifstream fin(permuFile);
	string token, marker;
	while (fin >> token) {
		string contig(token.begin()+1, token.end());
		auto genomeInfo = refContigs.count(contig) > 0 ? &ref : &tar ;
		auto contigInfo = make_pair(contig, vector<int>());
		while (fin >> marker) {
			if (marker == "$") break;
			contigInfo.second.push_back(std::stoi(marker));
			genomeInfo == &ref ?
				refContigs[contig] += 1 :
				tarContigs[contig] += 1 ;
		}	genomeInfo->push_back(contigInfo);
	}
}
void outputAll(string refFile, string tarFile, auto& ref, auto& tar, auto& reorder) {
	int uid = 1;
	ofstream fout(refFile);
	for (auto& [contig, vec]: ref)
		for (int marker: vec) {
			int newMarker = reorder[abs(marker)];
			if (newMarker != -1)
				fout << format("{} {} {} {}\n", uid++, (marker > 0 ? 1 : -1) * newMarker, contig, 1);
		}
	fout.close();

	uid = 1;
	fout.open(tarFile);
	for (auto& [contig, vec]: tar)
		for (int marker: vec) {
			int newMarker = reorder[abs(marker)];
			if (newMarker != -1)
				fout << format("{} {} {} {}\n", uid++, (marker > 0 ? 1 : -1) * newMarker, contig, 1);
		}
	fout.close();
}
void outputReducedContigs(string rmFile, auto& refContigs, auto& tarContigs) {
	int uid = 1;
	ofstream fout(rmFile);
	fout << format("ref:\n");
	for (auto& [contig, num]: refContigs)
		if (num == 0) fout << format("{} {}\n", uid++, contig);
	uid = 1;
	fout << format("\ntar:\n");
	for (auto& [contig, num]: tarContigs)
		if (num == 0) fout << format("{} {}\n", uid++, contig);
	fout.close();
}
int main(int argc, char* argv[]) {
	string refFNA(argv[1]);
	string tarFNA(argv[2]);
	string outDir(argv[3]);

	string command =
		format("Sibelia -k tools/parameter/paraset_bacterial -m 70 {} {} -o {}", refFNA, tarFNA, outDir);
	//	format("Sibelia -s far -m 20 {} {} -o {}", refFNA, tarFNA, outDir);
	if (system(command.c_str()) != 0) {
		print("[error] Cannot execute Sibelia commands!!\n");
		exit(1);
	}

	map<string, int> refContigs, tarContigs;
	vector<pair<string, vector<int>>> ref, tar;
	vector<int> reorder(1, -1);

	// obtain ref/tar contig num
	readFNA(refFNA, tarFNA, refContigs, tarContigs);

	// read permutations & blocks
	readPermu(outDir+"/genomes_permutations.txt", ref, tar, refContigs, tarContigs);
	readBlock(outDir+"/blocks_coords.txt", refContigs, tarContigs, reorder);

	outputAll(outDir+"/reference.all", outDir+"/target.all", ref, tar, reorder);
	outputReducedContigs(outDir+"/reduced.txt", refContigs, tarContigs);
	return 0;
}
