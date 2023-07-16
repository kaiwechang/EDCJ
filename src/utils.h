#ifndef DCJ_SCAFFOLDER_UTILS_H_
#define DCJ_SCAFFOLDER_UTILS_H_

// common data structures, macros, etc.

#include <fmt/core.h>
#include <fstream>
#include <vector>
#include <map>
#include <set>

using fmt::print;
using fmt::format;

using std::ifstream;
using std::ofstream;

using std::string;
using std::vector;
using std::map;
using std::set;

extern ofstream logFile;
enum Mode { MMDCJ, IDCJ, EDCJ };

inline int abs(int a) {
	return a > 0 ? a : -a ;
}
inline int max(int a, int b) {
	return a > b ? a : b ;
}
inline int min(int a, int b) {
	return a < b ? a : b ;
}
inline int sign(int a) {
	return (a > 0) - (a < 0);
}
// telos (signed id) to index
inline int idx(int t) {
	return abs(t)-1;
}
template<typename... Args>
void logging(std::string_view fstr, Args... args) {
	logFile << fmt::vformat(fstr, fmt::make_format_args(std::forward<Args>(args)...));
}
template<typename... Args>
void error(std::string_view fstr, Args... args) {
	print("[error] {}", fmt::vformat(fstr, fmt::make_format_args(std::forward<Args>(args)...)));
	exit(1);
}
struct Telos {
	// a telo is stored as its signed id ((+) for head, (-) for tail)
	int lhs, rhs;	// telo on the left/right hand side of a contig
	int ljs, rjs;	// telo on other contig joined with lhs/rhs, 0 when no join
	Telos(void): lhs{0}, rhs{0}, ljs{0}, rjs{0} {}
};
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
		logging("Marker({}, {}, {})\n", this->id, this->family, this->contig);
	}
};

// functions
int speedup_1	(string, string, string);
int speedup_1E	(string, string, string);
int speedup_2E	(string, string, string, bool);
int speedup_2ER	(string, string, string, bool);
int ilp_cap		(string, string, string, Mode, int, int, int);
int ilp_old		(string, string, string, Mode, int, int, int);
int ilp_new		(string, string, string, Mode, int, int, int);
int postprocess	(string, string, string);

#endif	// DCJ_SCAFFOLDER_UTILS_H_
