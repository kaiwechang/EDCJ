#include "utils.h"

void readGenome(string refFile, string tarFile, auto& ref, auto& tar, auto& contigSize, auto& refFamilySize, auto& tarFamilySize, int& maxFamily) {
	string contig;
	int id, family, tmp;

	ifstream fin(refFile);
	while (fin >> id >> family >> contig >> tmp) {
		maxFamily = max(maxFamily, abs(family));
		ref.push_back(Marker(id, family, contig));
	}	fin.close();
	fin.open(tarFile);
	while (fin >> id >> family >> contig >> tmp) {
		maxFamily = max(maxFamily, abs(family));
		tar.push_back(Marker(id, family, contig));
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
void speedup(auto& ref, auto& tar, auto& contigSize, auto& refFamilySize, auto& tarFamilySize, int maxFamily) {
	bool done = false;
	vector<int> refSurvived(maxFamily+1, -1);
	vector<int> tarSurvived(maxFamily+1, -1);
	auto setSurvived = [&](int i, int j, int idx, int jdx, int reversed) {
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
			(refSurvived[ref[idx].absFamily] != idx || tarSurvived[tar[jdx].absFamily] != jdx)) {

			done = false;
			refSurvived[ref[idx].absFamily] = idx;
			tarSurvived[tar[jdx].absFamily] = jdx;
		}
	};
	auto updateGenome = [&](void) {
		// update genome vectors
		int idx = 0;
		for (int i = 0; i < ref.size(); i++) {
			int kidx = refSurvived[ref[i].absFamily];
			if (kidx > 0 && kidx != i || tarFamilySize[ref[i].absFamily] == 0) {
				logging("del {}: {}\n", i+1, ref[i].absFamily);
				refFamilySize[ref[i].absFamily]--;
				continue;
			}
			if (kidx == i)
				refSurvived[ref[i].absFamily] = idx;
			ref[idx++].setMarker(idx, ref[i].family, ref[i].contig);
		}	ref.resize(idx);
		idx = 0;
		for (int i = 0; i < tar.size(); i++) {
			int kidx = tarSurvived[tar[i].absFamily];
			if (kidx > 0 && kidx != i || refFamilySize[tar[i].absFamily] == 0) {
				logging("del {}: {}\n", i+1, tar[i].absFamily);
				tarFamilySize[tar[i].absFamily]--;
				contigSize[tar[i].contig]--;
				continue;
			}
			if (kidx == i)
				tarSurvived[tar[i].absFamily] = idx;
			tar[idx++].setMarker(idx, tar[i].family, tar[i].contig);
		}	tar.resize(idx);
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
					setSurvived(i, j, i+1, j+1, 1);
					setSurvived(i, j, i-1, j-1, 1);
				} else {
					setSurvived(i, j, i+1, j-1, -1);
					setSurvived(i, j, i-1, j+1, -1);
				}
			}
		updateGenome();
		if (!done)
			logging("loop\n");
	}
}
void markerReorder(auto& ref, auto& tarFamilySize, int maxFamily, auto& reorder) {
	logging("maxFamily: {}\n", maxFamily);
	set<int> refFamily, tarFamily;
	reorder.assign(maxFamily+1, 0);
	for (Marker& m: ref)
		if (tarFamilySize[m.absFamily] != 0)
			refFamily.insert(m.absFamily);
	int uid = 1;
	for (int f: refFamily)
		reorder[f] = uid++;
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
void outputReducedContigs(string filename, auto& contigSize) {
	ofstream fout(filename);
	logging("contigSize:\n");
	for (auto& p: contigSize) {
		if (p.second == 0)
			fout << format("{}\n", p.first);
		logging("{}: {}\n", p.first, p.second);
	}	fout.close();
}
int main(int argc, char *argv[]) {
	if (argc < 4) {
		fmt::print("[error] Usage:\n>>> speedup_1E <ref genome> <tar genome> <output_dir>\n");
		return 0;
	}	string out_dir(argv[3]);
	logFile.open(out_dir+"/speedup_1.log");

	int maxFamily = 0;
	vector<Marker> ref, tar;
	vector<int> refFamilySize, tarFamilySize, reorder;
	map<string, int> contigSize;

	// read ref/tar
	readGenome(argv[1], argv[2], ref, tar, contigSize, refFamilySize, tarFamilySize, maxFamily);

	// speedup 1 for exemplar model
	speedup(ref, tar, contigSize, refFamilySize, tarFamilySize, maxFamily);
	markerReorder(ref, tarFamilySize, maxFamily, reorder);

	// output new ref/tar
	outputNewGenome(out_dir+"/ref_spd1.all", out_dir+"/tar_spd1.all", ref, tar, reorder);
	outputReducedContigs(out_dir+"/removed_spd1.txt", contigSize);

	logFile.close();
	return 0;
}
