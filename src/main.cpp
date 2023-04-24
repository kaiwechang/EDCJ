#include <filesystem>
#include <unistd.h>
#include "utils.h"

namespace fs = std::filesystem;

ofstream logFile;

void readOptions(int argc, char* argv[], string& refPath, string& tarPath, string& outDir, Mode& mode, bool& extended, double& gap, int& procs, int& timelimit, bool& spd1, bool& spd3, bool& rewrite, bool& align, bool& capping) {
	char option;
	while ((option = getopt(argc, argv, ":g:l:o:p:r:t:s:miexhacn")) != -1)
		switch (option) {
			case 'e':	mode = EDCJ;					break;
			case 'i':	mode = IDCJ;					break;
			case 'm':	mode = MMDCJ;					break;
			case 'a':	align = true;					break;	// for debugging spd3R
			case 'c':	capping = true;					break;	// for debugging capping
			case 'n':	rewrite = true;					break;	// for debugging ilp_R
			case 'x':	extended = true;				break;
			case 'o':	outDir = optarg;				break;
			case 'r':	refPath = optarg;				break;
			case 't':	tarPath = optarg;				break;
			case 'g':	gap = std::stod(optarg);		break;
			case 'p':	procs = std::stoi(optarg);		break;
			case 'l':	timelimit = std::stoi(optarg);	break;
			case 's':	spd1 = optarg == "1" || optarg == "13" ? true : false;
						spd3 = optarg == "3" || optarg == "13" ? true : false;
						break;
			case 'h':	print(
				"Usage:\n>>> ./DCJ_Scaffolder -r <ref genome> -t <tar genome> [optional options]\n\n"
				"required options:\n\n"
				" -r <ref genome> : reference genome path (\".all\" files from Sibelia)\n\n"
				" -t <tar genome> : target genome path (\".all\" files from Sibelia)\n\n"
				"optional options:\n\n"
				" -o <output dir> : output directory (default: output)\n\n"
				" -m / -i / -e    : maximum matching / intermediate / exemplar model \n"
				"                   (default: exemplar, prioity: e > i > m)\n"
				" -x              : extended version of speedup 3 (default: off)\n\n"
				" -g <gap>        : gap tolerance for Gurobi (default: 0.0001)\n\n"
				" -p <# of procs> : number of processors (default: all processors)\n\n"
				" -l <time limit> : time limit for Gurobi (default: 1800)\n\n"
				" -h              : usage instructions\n\n");	break;
			case '?':	error("Incorrect usage. Please use -h to see instructions.\n");
		}
	// check required options
	if (!fs::exists(refPath))
		error("reference genome \"{}\" not found !!!\n", refPath);
	if (!fs::exists(tarPath))
		error("target genome \"{}\" not found !!!\n", tarPath);
	if (!fs::exists(outDir))
		fs::create_directory(outDir);
}
void preprocess(string refFile, string tarFile, string outDir) {
	int maxFamily = 0;
	vector<Marker> ref, tar;
	vector<string> refOrder, tarOrder;

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
	}	fin.close();

	// contig order
	for (int i = 0; i < ref.size()-1; i++) {
		if (i == 0)
			refOrder.push_back(ref[i].contig);
		if (ref[i].contig != ref[i+1].contig)
			refOrder.push_back(ref[i+1].contig);
	}
	for (int i = 0; i < tar.size()-1; i++) {
		if (i == 0)
			tarOrder.push_back(tar[i].contig);
		if (tar[i].contig != tar[i+1].contig)
			tarOrder.push_back(tar[i+1].contig);
	}

	ofstream fout(outDir+"/info.txt");
	fout << format("          ref   tar\n");
	fout << format("marker: {:5d} {:5d}\n", ref.size(), tar.size());
	fout << format("contig: {:5d} {:5d}\n", refOrder.size(), tarOrder.size());
}
int main(int argc, char* argv[]) {
	double gap = -1;
	int procs = -1, timelimit = -1;
	string refPath, tarPath, outDir = "output";
	bool extended = false, spd1 = true, spd3 = true;
	bool rewrite = false, align = false, capping = false;
	Mode mode = EDCJ;

	readOptions(argc, argv, refPath, tarPath, outDir, mode, extended, gap, procs, timelimit,
				spd1, spd3, rewrite, align, capping);

	preprocess(refPath, tarPath, outDir);

	if (spd1)
		mode != EDCJ ? 
		speedup_1	(refPath, tarPath, outDir):
		speedup_1E	(refPath, tarPath, outDir);

	if (spd3)
		align == false ? 
		speedup_3E	(outDir+"/ref_spd1.all", outDir+"/tar_spd1.all", outDir, extended):
		speedup_3ER	(outDir+"/ref_spd1.all", outDir+"/tar_spd1.all", outDir, extended);

	rewrite == false ? capping == true ? 
		ilp_cap	(outDir+"/ref_spd3.all", outDir+"/tar_spd3.all", outDir, mode, gap, procs, timelimit):
		ilp_old	(outDir+"/ref_spd3.all", outDir+"/tar_spd3.all", outDir, mode, gap, procs, timelimit):
		ilp_new	(outDir+"/ref_spd3.all", outDir+"/tar_spd3.all", outDir, mode, gap, procs, timelimit);

	postprocess	(outDir+"/ref_spd3.all", outDir+"/tar_spd3.all", outDir);

	return 0;
}
