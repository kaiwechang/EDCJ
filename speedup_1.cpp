#include "utils.cpp"
#include <fstream>
#include <set>

using std::set;

bool done = false;
int maxFamily = 0;
vector<Marker> ref, tar;
vector<int> refFamilySize, tarFamilySize;

void rename(int i, int j, int idx, int jdx, int reversed)
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

		/*print("maxFamily: {}\n", maxFamily);
		print("ref[{}]: {}, \n", idx, ref[idx].family);
		print("tar[{}]: {}\n", jdx, tar[jdx].family);*/
	} 
}
int main(int argc, char *argv[])
{
	if (argc < 4) {
		print("[error] Usage:\n>>> speedup_1 <ref genome> <tar genome> <output_dir>\n");
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

	// speedup 1
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
			print("loop\n");
	}

	// marker reorder
	print("maxFamily: {}\n", maxFamily);
	set<int> refFamily, tarFamily;
	vector<int> reorder(maxFamily+1, 0);
	for (Marker& m: ref)
		if (tarFamilySize[m.absFamily] != 0)
			refFamily.insert(m.absFamily);
	int uid = 1;
	for (auto f: refFamily)
		reorder[f] = uid++;

	// file output
	string out_dir(argc < 4 ? "output" : argv[3]);
	std::ofstream fout(out_dir+"/ref_spd1.all");

	uid = 1;
	for (Marker& m: ref) {
		if (reorder[m.absFamily] != 0)
			fout << format("{} {} {} 1\n", uid++,
				(m.family > 0 ? 1 : -1)*reorder[m.absFamily], m.contig);
	}	fout.close();
	uid = 1;
	fout.open(out_dir+"/tar_spd1.all");
	for (Marker& m: tar) {
		if (reorder[m.absFamily] != 0)
			fout << format("{} {} {} 1\n", uid++,
				(m.family > 0 ? 1 : -1)*reorder[m.absFamily], m.contig);
	}	fout.close();
	return 0;
}
