#include <filesystem>
#include <unistd.h>
#include "utils.h"

namespace fs = std::filesystem;

void readOptions(int argc, char* argv[], string& refPath, string& tarPath, string& outDir, Mode& mode, bool& extended, double& gap, int& procs, int& timelimit) {
	char option;
	while ((option = getopt(argc, argv, ":g:l:m:o:p:r:t:miexh")) != -1)
		switch (option) {
			case 'm':	mode = MMDCJ;					break;
			case 'i':	mode = IDCJ;					break;
			case 'e':	mode = EDCJ;					break;
			case 'x':	extended = true;				break;
			case 'o':	outDir = optarg;				break;
			case 'r':	refPath = optarg;				break;
			case 't':	tarPath = optarg;				break;
			case 'g':	gap = std::stod(optarg);		break;
			case 'p':	procs = std::stoi(optarg);		break;
			case 'l':	timelimit = std::stoi(optarg);	break;
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
	if (outDir == "")
		fs::create_directory("output");
	else if (!fs::exists(outDir))
		fs::create_directory(outDir);
}
int main(int argc, char* argv[]) {
	double gap = -1;
	int procs = -1, timelimit = -1;
	string refPath, tarPath, outDir;
	bool extended = false;
	Mode mode = EDCJ;

	readOptions(argc, argv, refPath, tarPath, outDir, mode, extended, gap, procs, timelimit);

	return 0;
}
