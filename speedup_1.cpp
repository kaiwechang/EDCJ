#include "markerReorder.cpp"

#include <fmt/core.h>
#include <fstream>
#include <vector>

using fmt::print;
using fmt::format;
using std::vector;
using std::string;

#define max(a, b) (a > b ? a : b)

vector<Marker> refGenome, tarGenome;
vector<int> refFamilySize, tarFamilySize;
int refFamilyMax = 0, tarFamilyMax = 0, newFamily;

void rename(int i, int j, int idx, int jdx, int reversed, bool& done)
{
	if (// check index in range
		i   > 0 && i   < refGenome.size() &&
		idx > 0 && idx < refGenome.size() &&
		j   > 0 && j   < tarGenome.size() &&
		jdx > 0 && jdx < tarGenome.size() &&
		// same contig
		refGenome[i].contig == refGenome[idx].contig &&
		tarGenome[j].contig == tarGenome[jdx].contig &&
		// same family & (not) same sign
		refGenome[idx].family == reversed*tarGenome[jdx].family &&
		// at least one non-singleton
		(refFamilySize[refGenome[idx].absFamily] > 1 || tarFamilySize[tarGenome[jdx].absFamily] > 1)) {

		done = 0;
		// remove original
		refFamilySize[refGenome[idx].absFamily]--;
		tarFamilySize[tarGenome[jdx].absFamily]--;
		newFamily++;

		refGenome[idx].setFamily(newFamily);
		tarGenome[jdx].setFamily(
			refGenome[i].family == -tarGenome[j].family ?
			-newFamily : newFamily);

		if (refFamilySize.size() < newFamily+1)
			refFamilySize.resize(newFamily+1, 0);
		refFamilySize[newFamily] = 1;

		if (tarFamilySize.size() < newFamily+1)
			tarFamilySize.resize(newFamily+1, 0);
		tarFamilySize[newFamily] = 1;

		/*print("newFamily: {}\n", newFamily);
		print("ref[{}]: {}\n", idx, refGenome[idx].family);
		print("tar[{}]: {}\n", idx, tarGenome[jdx].family);*/
	} 
}
int main(int argc, char *argv[])
{
	if (argc < 4) {
		print("[error] no ref/tar genome.\n");
		return 0;
	}
	string contig;
	int id, family, tmp;
	std::ifstream fin(argv[1]);
	while (fin >> id >> family >> contig >> tmp) {
		refFamilyMax = max(refFamilyMax, abs(family));
		refGenome.push_back(Marker(id, family, contig));
	}	fin.close();
	fin.open(argv[2]);
	while (fin >> id >> family >> contig >> tmp) {
		tarFamilyMax = max(tarFamilyMax, abs(family));
		tarGenome.push_back(Marker(id, family, contig));
	}	fin.close();

	// count family size
	newFamily = max(refFamilyMax, tarFamilyMax);
	refFamilySize.assign(newFamily+1, 0);
	tarFamilySize.assign(newFamily+1, 0);
	for (Marker& m: refGenome)
		refFamilySize[m.absFamily]++;
	for (Marker& m: tarGenome)
		tarFamilySize[m.absFamily]++;

	bool done = 0;
	while (!done) {
		done = 1;
		for (int i = 0; i < refGenome.size(); i++)
			for (int j = 0; j < tarGenome.size(); j++) {
				// ref[i] and tar[j] not both singlton and same family
				if (refFamilySize[refGenome[i].absFamily] > 1 ||
					tarFamilySize[tarGenome[j].absFamily] > 1 ||
					refGenome[i].absFamily != tarGenome[j].absFamily)
					continue;
				// ref[i] and tar[j] same sign
				if (refGenome[i].family == tarGenome[j].family) {
					rename(i, j, i+1, j+1, 1, done);
					rename(i, j, i-1, j-1, 1, done);
				} else {
					rename(i, j, i+1, j-1, -1, done);
					rename(i, j, i-1, j+1, -1, done);
				}
			}
		if (!done)
			print("rename\n");
	}
	string out_dir(argc < 4 ? "output" : argv[3]);
	std::ofstream fout(out_dir+"/ref_spd1.all");

	int uid = 1;
	for (int i = 0; i < refGenome.size(); i++) {
		if (tarFamilySize[refGenome[i].absFamily] != 0)
			fout << format("{} {} {} 1\n", uid++, refGenome[i].family, refGenome[i].contig);
	}	fout.close();
	uid = 1;
	fout.open(out_dir+"/tar_spd1.all");
	for (int i = 0; i < tarGenome.size(); i++) {
		if (refFamilySize[tarGenome[i].absFamily] != 0)
			fout << format("{} {} {} 1\n", uid++, tarGenome[i].family, tarGenome[i].contig);
	}	fout.close();

	markerReorder(out_dir+"/ref_spd1.all", out_dir+"/tar_spd1.all");
	return 0;
}
