#ifndef DCJ_SCAFFOLDER_UTILS_H_
#define DCJ_SCAFFOLDER_UTILS_H_

// common data structures, macros, etc.

#include <fmt/core.h>
#include <fstream>
#include <vector>
#include <map>

using fmt::print;
using fmt::format;

using std::ifstream;
using std::ofstream;

using std::string;
using std::vector;
using std::map;

#define max(a, b) (a > b ? a : b)
#define abs(a) (a > 0 ? a : -a)

ofstream logFile;

template<typename... Args>
void logging(std::string_view fstr, Args... args) {
	logFile << fmt::vformat(fstr, fmt::make_format_args(std::forward<Args>(args)...));
}
struct Telos
{
	int lhs, rhs;
	Telos(void): lhs{0}, rhs{0} {}
	Telos(int lhs, int rhs): lhs{lhs}, rhs{rhs} {}
};
struct Marker
{
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
