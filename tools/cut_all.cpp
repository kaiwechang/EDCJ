#include <vector>
#include <unordered_map>

#include <fstream>
#include <sstream>
#include <fmt/core.h>

using fmt::print;
using fmt::format;

using std::vector;
using std::string;
using std::unordered_map;

using std::getline;
using std::ifstream;
using std::ofstream;
using std::stringstream;

inline int abs(int a) {
	return a > 0 ? a : -a ;
}
inline int max(int a, int b) {
	return a > b ? a : b ;
}
struct Marker {
	int id, family, absFamily;
	string contig;

	Marker(void):
		id{0}, family{0}, absFamily{0}, contig{""} {}
	Marker(int id, int family, string contig):
		id{id}, family{family}, absFamily{abs(family)}, contig{contig} {}
	void setMarker(int newId, int newFamily, string newContig) {
		this->id = newId;
		this->family = newFamily;
		this->contig = newContig;
		this->absFamily = abs(newFamily);
	}
	void setFamily(int newFamily) {
		this->family = newFamily;
		this->absFamily = abs(newFamily);
	}
	void show() {
		print("Marker({}, {}, {})\n", this->id, this->family, this->contig);
	}
};
void readALL(string refFile, string tarFile, auto& ref, auto& tar) {
	string contig;
	int id, family, temp;

	ifstream fin(refFile);
	while (fin >> id >> family >> contig >> temp) {
		ref.push_back(Marker(id, family, contig));
	}	fin.close();

	fin.open(tarFile);
	while (fin >> id >> family >> contig >> temp) {
		tar.push_back(Marker(id, family, contig));
	}	fin.close();
}
void cutAnswer(string inFile, string outFile, int start, int end, auto& keep) {
	int uid = 0;
	string line, token;
	ifstream fin(inFile);
	ofstream fout(outFile);

	while (getline(fin, line)) {
		if (line[0] == '>') {
			fout << line << "\n";
		} else {
			stringstream ss(line); ss >> token;
			//if (uid >= start && uid < end) {
			if (uid%10 == 0) {
				keep[token] = 1;
				fout << line << "\n";
			}	uid++;
		}
	}	fin.close();
	fout.close();
}
void outputAll(string refFile, string tarFile, auto& ref, auto& tar, auto& keep, auto& reorder) {
	int uid = 1;
	for (auto& m: tar)
		if (keep.count(m.contig) && !reorder.count(m.absFamily))
				reorder[m.absFamily] = uid++;

	uid = 1;
	ofstream fout(refFile);
	for (auto& m: ref) {
		if (reorder.count(m.absFamily))
			fout << format("{} {} {} 1\n", uid++,
					(m.family > 0 ? 1 : -1) * reorder[m.absFamily], m.contig, 1);
	}	fout.close();

	uid = 1;
	fout.open(tarFile);
	for (auto& m: tar) {
		if (keep.count(m.contig) && reorder.count(m.absFamily))
			fout << format("{} {} {} 1\n", uid++, 
					(m.family > 0 ? 1 : -1) * reorder[m.absFamily], m.contig, 1);
	}	fout.close();
}
int main(int argc, char* argv[]) {
	if (argc < 3) {
		print("[error] Usage:\n>>> cut_all <target dir> <output dir>\n");
		print("[error] target dir: directory containing reference.all / target.all / answerToAll\n");
		return 0;
	}
	string tarDir(argv[1]);
	string outDir(argv[2]);

	vector<Marker> ref, tar;
	unordered_map<int, int> reorder;
	unordered_map<string, int> keep;

	readALL(tarDir+"/reference.all", tarDir+"/target.all", ref, tar);
	cutAnswer(tarDir+"/answerToAll", outDir+"/answerToAll", 18000, 20000, keep);
	outputAll(outDir+"/reference.all", outDir+"/target.all", ref, tar, keep, reorder);

	// around 20000 contigs on HS_data
	print("ref size: {}\n", ref.size());
	print("tar size: {}\n", tar.size());

	return 0;
}
