#include "utils.h"

void readGenome(string refFile, string tarFile, auto& ref, auto& tar, auto& contigSize, auto& refFamilySize, auto& tarFamilySize, int& maxFamily) {
	string contig;
	int id, family, tmp;

	ifstream fin(refFile);
	while (fin >> id >> family >> contig >> tmp) {
		ref.push_back(Marker(id, family, contig));
		maxFamily = max(maxFamily, abs(family));
	}	fin.close();
	fin.open(tarFile);
	while (fin >> id >> family >> contig >> tmp) {
		tar.push_back(Marker(id, family, contig));
		maxFamily = max(maxFamily, abs(family));
		contigSize[contig]++;
	}	fin.close();

	// count family size
	refFamilySize.assign(maxFamily+1, 0);
	tarFamilySize.assign(maxFamily+1, 0);
	for (Marker& m: ref)
		refFamilySize[m.absFamily]++;
	for (Marker& m: tar)
		tarFamilySize[m.absFamily]++;
}
void speedup(auto& ref, auto& tar, auto& refFamilySize, auto& tarFamilySize, int& maxFamily) {
	bool done = false;
	auto rename = [&](int i, int j, int idx, int jdx, int reversed) {
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
			// at least one non-singleton
			(refFamilySize[ref[idx].absFamily] > 1 || tarFamilySize[tar[jdx].absFamily] > 1)) {

			done = false;
			// remove original
			refFamilySize[ref[idx].absFamily]--;
			tarFamilySize[tar[jdx].absFamily]--;
			maxFamily++;

			ref[idx].setFamily(maxFamily);
			tar[jdx].setFamily(reversed*maxFamily);

			if (refFamilySize.size() < maxFamily+1)
				refFamilySize.resize(maxFamily+1, 0);
			refFamilySize[maxFamily] = 1;

			if (tarFamilySize.size() < maxFamily+1)
				tarFamilySize.resize(maxFamily+1, 0);
			tarFamilySize[maxFamily] = 1;

			/*logging("maxFamily: {}\n", maxFamily);
			logging("ref[{}]: {}, \n", idx, ref[idx].family);
			logging("tar[{}]: {}\n", jdx, tar[jdx].family);*/
		} 
	};

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
					rename(i, j, i+1, j+1, 1);
					rename(i, j, i-1, j-1, 1);
				} else {
					rename(i, j, i+1, j-1, -1);
					rename(i, j, i-1, j+1, -1);
				}
			}
		if (!done)
			logging("loop\n");
	}
}
auto markerReorder(auto& ref, auto& tarFamilySize, int maxFamily) {
	logging("maxFamily: {}\n", maxFamily);
	set<int> refFamily, tarFamily;
	vector<int> reorder(maxFamily+1, 0);
	for (Marker& m: ref)
		if (tarFamilySize[m.absFamily] != 0)
			refFamily.insert(m.absFamily);
	int uid = 1;
	for (int f: refFamily)
		reorder[f] = uid++;
	return reorder;
}
void outputNewGenome(string refFile, string tarFile, auto& ref, auto& tar, auto& reorder) {
	int uid = 1;
	ofstream fout(refFile);
	for (Marker& m: ref) {
		if (reorder[m.absFamily] != 0)
			fout << format("{} {} {} 1\n", uid++,
				(m.family > 0 ? 1 : -1)*reorder[m.absFamily], m.contig);
	}	fout.close();
	uid = 1;
	fout.open(tarFile);
	for (Marker& m: tar) {
		if (reorder[m.absFamily] != 0)
			fout << format("{} {} {} 1\n", uid++,
				(m.family > 0 ? 1 : -1)*reorder[m.absFamily], m.contig);
	}	fout.close();
}
void outputReducedContigs(string filename, auto& tar, auto& contigSize, auto& refFamilySize) {
	ofstream fout(filename);
	logging("contigSize:\n");
	for (Marker& m: tar)
		if (refFamilySize[m.absFamily] == 0)
			contigSize[m.contig]--;
	for (auto& p: contigSize) {
		if (p.second == 0)
			fout << format("{}\n", p.first);
		logging("{}: {}\n", p.first, p.second);
	}	fout.close();
}
int main(int argc, char *argv[]) {
	if (argc < 4) {
		fmt::print("[error] Usage:\n>>> speedup_1 <ref genome> <tar genome> <output_dir>\n");
		return 0;
	}	string out_dir(argv[3]);
	logFile.open(out_dir+"/speedup_1.log");

	int maxFamily = 0;
	vector<Marker> ref, tar;
	vector<int> refFamilySize, tarFamilySize, reorder;
	map<string, int> contigSize;

	// read ref/tar
	readGenome(argv[1], argv[2], ref, tar, contigSize, refFamilySize, tarFamilySize, maxFamily);

	// speedup 1
	speedup(ref, tar, refFamilySize, tarFamilySize, maxFamily);
	reorder = markerReorder(ref, tarFamilySize, maxFamily);

	// output new ref/tar
	outputNewGenome(out_dir+"/ref_spd1.all", out_dir+"/tar_spd1.all", ref, tar, reorder);
	outputReducedContigs(out_dir+"/removed_spd1.txt", tar, contigSize, refFamilySize);

	logFile.close();
	return 0;
}
