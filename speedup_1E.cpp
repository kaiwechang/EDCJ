#include "markerReorder.cpp"
#include <iostream>
#include <fstream>
#include <set>
#include <vector>
using namespace std;
vector<int> refFamily, tarFamily; //start with index 0(input marker id = vector index+1)
vector<string> refContig, tarContig;
vector<int> refMarkerNum, tarMarkerNum;
vector<int> ref_keep, tar_keep;
int refMarkerMax = 0, tarMarkerMax = 0;
void showSingletonPercentage()
{
	int cnt = 0, tt = 0;
	for (int i = 0; i < refFamily.size(); i++)
	{
		if (refMarkerNum[abs(refFamily[i])] == 1 && tarMarkerNum[abs(refFamily[i])] == 1)
			cnt += 2;
	}
	for (int i = 0; i < refFamily.size(); i++)
	{
		if (tarMarkerNum[abs(refFamily[i])] != 0)
			tt++;
	}
	for (int i = 0; i < tarFamily.size(); i++)
	{
		if (refMarkerNum[abs(tarFamily[i])] != 0)
			tt++;
	}
	cout << cnt << " " << tt << endl;
	cout << cnt / (double)tt << endl;
}
int main(int argc, char *argv[])
{
	if (argc < 4) {
		printf("[error] no ref/tar genome.\n");
		exit(1);
	}
	ifstream fin(argv[1]);
	int id, family, tmp;
	string contig;
	while (fin >> id >> family >> contig >> tmp) {
		refFamily.push_back(family);
		refContig.push_back(contig);
	}
	fin.close();
	fin.open(argv[2]);
	while (fin >> id >> family >> contig >> tmp) {
		tarFamily.push_back(family);
		tarContig.push_back(contig);
	}

	//count markernum
	refMarkerNum.resize(refFamily.size() + 1);
	for (int i : refFamily) {
		i = abs(i);
		refMarkerMax = refMarkerMax > i ? refMarkerMax : i;
		refMarkerNum[i]++;
	}
	tarMarkerNum.resize(tarFamily.size() + 1);
	for (int i : tarFamily) {
		i = abs(i);
		tarMarkerMax = tarMarkerMax > i ? tarMarkerMax : i;
		tarMarkerNum[i]++;
	}	showSingletonPercentage();

	ref_keep.assign(refMarkerMax+1, -1);
	tar_keep.assign(tarMarkerMax+1, -1);
	for (int i = 0; i < refFamily.size() - 1; i++) {
		for (int j = 0; j < tarFamily.size() - 1; j++) {
			// ref[i] and tar[j] not both singlton and same family
			if (refMarkerNum[abs(refFamily[i])] > 1 || tarMarkerNum[abs(tarFamily[j])] > 1 ||
				abs(refFamily[i]) != abs(tarFamily[j]))
				continue;
			// ref[i] and tar[j] same sign
			if (refFamily[i] == tarFamily[j]) {
				if (refContig[i] == refContig[i+1] && tarContig[j] == tarContig[j+1] &&
					refFamily[i+1] == tarFamily[j+1]) {
					int id1 = i+1, id2 = j+1;
					ref_keep[abs(refFamily[id1])] = id1;
					tar_keep[abs(tarFamily[id2])] = id2;
				}
				if (i > 0 && j > 0 &&
					refContig[i-1] == refContig[i] && tarContig[j-1] == tarContig[j] &&
					refFamily[i-1] == tarFamily[j-1]) {
					int id1 = i-1, id2 = j-1;
					ref_keep[abs(refFamily[id1])] = id1;
					tar_keep[abs(tarFamily[id2])] = id2;
				}
			} else {
				if (j > 0 &&
					refContig[i] == refContig[i+1] && tarContig[j] == tarContig[j-1] &&
					refFamily[i+1] == -tarFamily[j-1]) {
					int id1 = i+1, id2 = j-1;
					ref_keep[abs(refFamily[id1])] = id1;
					tar_keep[abs(tarFamily[id2])] = id2;
				}
				if (i > 0 &&
					refContig[i-1] == refContig[i] && tarContig[j] == tarContig[j+1] &&
					refFamily[i-1] == -tarFamily[j+1]) {
					int id1 = i-1, id2 = j+1;
					ref_keep[abs(refFamily[id1])] = id1;
					tar_keep[abs(tarFamily[id2])] = id2;
				}
			}
		}
	}
	string out_dir(argc < 4 ? "output" : argv[3]);
	ofstream fout(out_dir+"/ref_spd1.all");

	int uid = 1;
	for (int i = 0; i < refFamily.size(); i++) {
		int kp_idx = ref_keep[abs(refFamily[i])];
		if (kp_idx > 0 && kp_idx != i) {
			printf("%d: %d\n", i+1, abs(refFamily[i]));
			continue;
		}
		fout << uid++ << " " << refFamily[i] << " " << refContig[i] << " " << 1 << endl;
	}	fout.close();

	uid = 1;
	fout.open(out_dir+"/tar_spd1.all");
	for (int i = 0; i < tarFamily.size(); i++) {
		int kp_idx = tar_keep[abs(tarFamily[i])];
		if (kp_idx > 0 && kp_idx != i) {
			printf("%d: %d\n", i+1, abs(tarFamily[i]));
			continue;
		}
		fout << uid++ << " " << tarFamily[i] << " " << tarContig[i] << " " << 1 << endl;
	}
	showSingletonPercentage();
	markerReorder(out_dir+"/ref_spd1.all", out_dir+"/tar_spd1.all");
	return 0;
}
