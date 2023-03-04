#include <iostream>
#include <fstream>
#include <vector>
using namespace std;

#define max(a, b) (a > b ? a : b)
#define abs(a) (a > 0 ? a : -a)

struct Marker
{
	int id, family, absFamily;
	string contig;
	Marker(int id_, int family_, string contig_) {
		id = id_;
		family = family_;
		absFamily = abs(family_);
		contig = contig_;
	}
	void show() {
		printf("Marker(%d, %d, %s)\n", id, family, contig.c_str());
	}
};
vector<Marker> refGenome, tarGenome;
vector<int> refKeep, tarKeep;
vector<int> refFamilySize, tarFamilySize;
int refFamilyMax = 0, tarFamilyMax = 0;

void setKeep(int i, int j, int idx, int jdx, int reversed)
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
		refGenome[idx].family == reversed*tarGenome[jdx].family) {

		refKeep[refGenome[idx].absFamily] = idx;
		tarKeep[tarGenome[jdx].absFamily] = jdx;
	}
}
int main(int argc, char *argv[])
{
	if (argc < 4) {
		printf("[error] no ref/tar genome.\n");
		return 0;
	}
	string contig;
	int id, family, tmp;
	ifstream fin(argv[1]);
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
	refFamilySize.assign(refFamilyMax+1, 0);
	tarFamilySize.assign(tarFamilyMax+1, 0);
	for (Marker& m: refGenome)
		refFamilySize[m.absFamily]++;
	for (Marker& m: tarGenome)
		tarFamilySize[m.absFamily]++;

	refKeep.assign(refFamilyMax+1, -1);
	tarKeep.assign(tarFamilyMax+1, -1);
	for (int i = 0; i < refGenome.size(); i++)
		for (int j = 0; j < tarGenome.size(); j++) {
			// ref[i] and tar[j] not both singlton and same family
			if (refFamilySize[refGenome[i].absFamily] > 1 ||
				tarFamilySize[tarGenome[j].absFamily] > 1 ||
				refGenome[i].absFamily != tarGenome[j].absFamily)
				continue;
			// ref[i] and tar[j] same sign
			if (refGenome[i].family == tarGenome[j].family) {
				setKeep(i, j, i+1, j+1, 1);
				setKeep(i, j, i-1, j-1, 1);
			} else {
				setKeep(i, j, i+1, j-1, -1);
				setKeep(i, j, i-1, j+1, -1);
			}
		}
	string out_dir(argc < 4 ? "output" : argv[3]);
	ofstream fout(out_dir+"/ref_spd1.all");

	int uid = 1;
	for (int i = 0; i < refGenome.size(); i++) {
		int kidx = refKeep[refGenome[i].absFamily];
		if (kidx > 0 && kidx != i) {
			printf("del %d: %d\n", i+1, refGenome[i].absFamily);
			continue;
		}
		fout << uid++ << " " << refGenome[i].family << " " << refGenome[i].contig << " " << 1 << endl;
	}	fout.close();
	uid = 1;
	fout.open(out_dir+"/tar_spd1.all");
	for (int i = 0; i < tarGenome.size(); i++) {
		int kidx = tarKeep[tarGenome[i].absFamily];
		if (kidx > 0 && kidx != i) {
			printf("del %d: %d\n", i+1, tarGenome[i].absFamily);
			continue;
		}
		fout << uid++ << " " << tarGenome[i].family << " " << tarGenome[i].contig << " " << 1 << endl;
	}	fout.close();
	return 0;
}
