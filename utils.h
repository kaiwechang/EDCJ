#ifndef DCJ_SCAFFOLDER_UTILS_H_
#define DCJ_SCAFFOLDER_UTILS_H_

// common data structures, macros, etc.

#include <fmt/core.h>
#include <fstream>
#include <vector>
#include <map>
#include <set>

using fmt::format;

using std::ifstream;
using std::ofstream;

using std::string;
using std::vector;
using std::map;
using std::set;

#define abs(a) (a > 0 ? a : -a)
#define max(a, b) (a > b ? a : b)

#define sign(a) ((a > 0) - (a < 0))
#define idx(t) (abs(t)-1)	// telos (signed id) to index

ofstream logFile;

template<typename... Args>
void logging(std::string_view fstr, Args... args) {
	logFile << fmt::vformat(fstr, fmt::make_format_args(std::forward<Args>(args)...));
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

#endif	// DCJ_SCAFFOLDER_UTILS_H_
