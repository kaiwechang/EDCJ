#include "utils.cpp"
#include <fstream>

bool done = false;
int maxFamily = 0;
vector<Marker> ref, tar;
vector<int> refKeep, tarKeep;
vector<int> refFamilySize, tarFamilySize;

void setKeep(int i, int j, int idx, int jdx, int reversed)
{
	if (// check index in range
		i   > 0 && i   < ref.size() &&
		idx > 0 && idx < ref.size() &&
		j   > 0 && j   < tar.size() &&
		jdx > 0 && jdx < tar.size() &&
		// same contig
		ref[i].contig == ref[idx].contig &&
		tar[j].contig == tar[jdx].contig &&
		// same family & (not) same sign
		ref[idx].family == reversed*tar[jdx].family &&
		(refKeep[ref[idx].absFamily] != idx || tarKeep[tar[jdx].absFamily] != jdx)) {

		done = false;
		refKeep[ref[idx].absFamily] = idx;
		tarKeep[tar[jdx].absFamily] = jdx;
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
		maxFamily = max(maxFamily, abs(family));
		ref.push_back(Marker(id, family, contig));
	}	fin.close();
	fin.open(argv[2]);
	while (fin >> id >> family >> contig >> tmp) {
		maxFamily = max(maxFamily, abs(family));
		tar.push_back(Marker(id, family, contig));
	}	fin.close();

	// count family size
	refFamilySize.assign(maxFamily+1, 0);
	tarFamilySize.assign(maxFamily+1, 0);
	for (Marker& m: ref)
		refFamilySize[m.absFamily]++;
	for (Marker& m: tar)
		tarFamilySize[m.absFamily]++;

	// speedup_1 for exemplar
	refKeep.assign(maxFamily+1, -1);
	tarKeep.assign(maxFamily+1, -1);
	while (!done) {
		done = true;
		for (int i = 0; i < ref.size(); i++)
			for (int j = 0; j < tar.size(); j++) {
				// ref[i] and tar[j] not both singlton and same family
				if (refFamilySize[ref[i].absFamily] > 1 ||
					tarFamilySize[tar[j].absFamily] > 1 ||
					ref[i].absFamily != tar[j].absFamily)
					continue;
				// ref[i] and tar[j] same sign
				if (ref[i].family == tar[j].family) {
					setKeep(i, j, i+1, j+1, 1);
					setKeep(i, j, i-1, j-1, 1);
				} else {
					setKeep(i, j, i+1, j-1, -1);
					setKeep(i, j, i-1, j+1, -1);
				}
			}
		// update genome vectors
		int idx = 0;
		for (int i = 0; i < ref.size(); i++) {
			int kidx = refKeep[ref[i].absFamily];
			if (kidx > 0 && kidx != i || tarFamilySize[ref[i].absFamily] == 0) {
				print("del {}: {}\n", i+1, ref[i].absFamily);
				refFamilySize[ref[i].absFamily]--;
				continue;
			}
			if (kidx == i)
				refKeep[ref[i].absFamily] = idx;
			ref[idx++].setMarker(idx, ref[i].family, ref[i].contig);
		}	ref.resize(idx);
		idx = 0;
		for (int i = 0; i < tar.size(); i++) {
			int kidx = tarKeep[tar[i].absFamily];
			if (kidx > 0 && kidx != i || refFamilySize[tar[i].absFamily] == 0) {
				print("del {}: {}\n", i+1, tar[i].absFamily);
				tarFamilySize[tar[i].absFamily]--;
				continue;
			}
			if (kidx == i)
				tarKeep[tar[i].absFamily] = idx;
			tar[idx++].setMarker(idx, tar[i].family, tar[i].contig);
		}	tar.resize(idx);
		if (!done)
			print("loop\n");
	}

	// file output
	string out_dir(argc < 4 ? "output" : argv[3]);
	std::ofstream fout(out_dir+"/ref_spd1.all");
	for (Marker& m: ref) {
		fout << format("{} {} {} 1\n", m.id, m.family, m.contig);
	}	fout.close();
	fout.open(out_dir+"/tar_spd1.all");
	for (Marker& m: tar) {
		fout << format("{} {} {} 1\n", m.id, m.family, m.contig);
	}	fout.close();
	return 0;
}
