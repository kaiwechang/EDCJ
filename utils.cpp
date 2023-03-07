#ifndef DCJ_SCAFFOLDER_UTILS_CPP_
#define DCJ_SCAFFOLDER_UTILS_CPP_

// common data structures, macros, etc.

#include <fmt/core.h>
#include <vector>

using fmt::print;
using fmt::format;
using std::string;
using std::vector;

#define max(a, b) (a > b ? a : b)
#define abs(a) (a > 0 ? a : -a)

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
		print("Marker({}, {}, {})\n", this->id, this->family, this->contig);
	}
};

#endif	// DCJ_SCAFFOLDER_UTILS_CPP_
