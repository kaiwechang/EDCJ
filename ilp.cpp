#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <set>
#include <map>
#include "string.h"
#include "gurobi_c++.h"

#define pb push_back
#define capFamily 999999
using namespace std;

struct Edge
{
	int etype; //0=adj(ref) 1=adj(tar) 2=homo
	bool is_potential;
	int node1, node2;
	Edge(bool is_potential_, int node1_, int node2_, int etype_)
	{
		is_potential = is_potential_;
		node1 = node1_;
		node2 = node2_;
		etype = etype_;
	}
	void show()
	{
		cout << "Edge(" << is_potential << "," << node1 << "," << node2 << "," << etype << ")" << endl;
	}
};
struct Marker
{
	int id, family;
	string contig;
	Marker(int id_, int family_, string contig_)
	{
		id = id_;
		family = family_;
		contig = contig_;
	}
	void show()
	{
		cout << "Marker(" << id << "," << family << "," << contig << ")" << endl;
	}
};
struct Node
{
	int nodeid, markerid;
	int family;
	bool is_head;
	string contig;
	Node(int nodeid_, int markerid_, int family_, string contig_, bool is_head_)
	{
		nodeid = nodeid_;
		markerid = markerid_;
		if (family_ < 0)
			family = -family_;
		else
			family = family_;
		contig = contig_;
		is_head = is_head_;
	}
	void show()
	{
		cout << "Node(" << nodeid << "," << markerid << "," << family << "," << contig << "," << is_head << ")" << endl;
	}
	void showSpaceVersion()
	{
		cout << nodeid << " " << markerid << " " << family << " " << contig << " " << is_head << endl;
	}
};
struct Adjacency
{
	bool is_potential;
	int id1, id2;
	Adjacency(bool is_potential_, int id1_, int id2_)
	{
		is_potential = is_potential_;
		id1 = id1_;
		id2 = id2_;
	}
};

string out_dir;
int refAdjSize = 0, tarAdjSize = 0;
bool isdebug = 0, isdebugadj = 0, isdebugall = 0;
vector<Marker> refMarkers, tarMarkers;
vector<GRBVar> markerVars; // ref--tar--capRef--capTar
vector<GRBVar> homoVars;
vector<GRBVar> cycleVars;
vector<GRBVar> nodeVars; // ref--capRef--tar--capTar
vector<GRBVar> adjVars;
vector<Node> nodeVecs;	// ref--capRef--tar--capTar
vector<Edge> homoEdges; // edge (,markeri,markerj,)
vector<Edge> adjEdges;	// edge (,nodei,nodej,)
vector<int> refFamSize, tarFamSize;
vector<int> refTelomeres, tarTelomeres;
vector<string> refContig, tarContig;
set<int> familySet;
map<pair<int, int>, bool> refPotentialAdj, tarPotentialAdj; //is there a contig with marker pair int,int
map<pair<string, string>, bool> refAdj, tarAdj;				//is contig pair string,string available
GRBLinExpr teloWithCapSum = 0;
GRBLinExpr refRemoveMarkerNum = 0;
int refCapSize = 0, tarCapSize = 0, rc = 0, tc = 0;
int refNodeSize;
int refMarkerSize = 0, tarMarkerSize = 0; //without capping gene
int refMarkerNum = 0, tarContigNum = 0, refContigNum = 0;
int refRealAdjNumc2c = 0, tarRealAdjNumc2c = 0;
int refPotAdjNumt2c = 0, refPotAdjNumt2t = 0, refPotAdjNumc2c = 0, tarPotAdjNumt2c = 0, tarPotAdjNumt2t = 0, tarPotAdjNumc2c = 0;
int jabs(int n)
{
	if (n > 0)
		return n;
	return -n;
}
void getInput(string refFile, string tarFile)
{
	ifstream fin;
	fin.open(refFile);
	int id, family, tmp;
	string contig;
	string prevContig = "NotAContig";
	vector<int> families;
	while (fin >> id >> family >> contig >> tmp)
	{
		Marker marker(id, family, contig);
		if (contig != prevContig)
		{
			refContigNum++;
			for (auto it1 = families.begin(); it1 != families.end(); it1++)
			{
				for (auto it2 = it1 + 1; it2 != families.end(); it2++)
				{
					pair<int, int> famPair;
					if (*it1 < *it2)
						famPair = make_pair(*it1, *it2);
					else
						famPair = make_pair(*it2, *it1);
					refPotentialAdj[famPair] = 1;
				}
			}
			prevContig = contig;
			refContig.pb(contig);
			families.clear();
		}
		families.pb(jabs(family));

		if (isdebug)
			marker.show();
		refMarkers.pb(marker);
		refMarkerNum++;
	}
	for (auto it1 = families.begin(); it1 != families.end(); it1++)
	{
		for (auto it2 = it1 + 1; it2 != families.end(); it2++)
		{
			pair<int, int> famPair;
			if (*it1 < *it2)
				famPair = make_pair(*it1, *it2);
			else
				famPair = make_pair(*it2, *it1);
			refPotentialAdj[famPair] = 1;
		}
	}
	refMarkerSize = refMarkers.size();
	fin.close();
	fin.open(tarFile);
	prevContig = "NotAContig";
	families.clear();
	while (fin >> id >> family >> contig >> tmp)
	{
		Marker marker(id, family, contig);
		if (contig != prevContig)
		{
			tarContigNum++;
			for (auto it1 = families.begin(); it1 != families.end(); it1++)
			{
				for (auto it2 = it1 + 1; it2 != families.end(); it2++)
				{
					pair<int, int> famPair;
					if (*it1 < *it2)
						famPair = make_pair(*it1, *it2);
					else
						famPair = make_pair(*it2, *it1);
					tarPotentialAdj[famPair] = 1;
				}
			}
			prevContig = contig;
			tarContig.pb(contig);
			families.clear();
		}
		families.pb(jabs(family));
		if (isdebug)
			marker.show();
		tarMarkers.pb(marker);
	}
	for (auto it1 = families.begin(); it1 != families.end(); it1++)
	{
		for (auto it2 = it1 + 1; it2 != families.end(); it2++)
		{
			pair<int, int> famPair;
			if (*it1 < *it2)
				famPair = make_pair(*it1, *it2);
			else
				famPair = make_pair(*it2, *it1);
			tarPotentialAdj[famPair] = 1;
		}
	}
	tarMarkerSize = tarMarkers.size();
	fin.close();
}
void addDuplicatedVar(GRBModel &model)
{
	for (auto m : refMarkers)
		familySet.insert(jabs(m.family));

	auto it = familySet.end();
	it--;
	refFamSize.assign(*it + 1, 0), tarFamSize.assign(*it + 1, 0);
	for (auto m : refMarkers)
		refFamSize[jabs(m.family)]++;
	for (auto m : tarMarkers)
		tarFamSize[jabs(m.family)]++;

	for (auto m : refMarkers)
	{
		if (refFamSize[jabs(m.family)] == 1)
		{
			GRBVar x = model.addVar(0, 0, 0, GRB_BINARY);
			refRemoveMarkerNum += x;
			markerVars.pb(x);
		}
		else
		{
			GRBVar x = model.addVar(0, 1, 0, GRB_BINARY);
			refRemoveMarkerNum += x;
			markerVars.pb(x);
		}
	}
	for (auto m : tarMarkers)
	{
		if (tarFamSize[jabs(m.family)] == 1)
		{
			GRBVar x = model.addVar(0, 0, 0, GRB_BINARY);
			markerVars.pb(x);
		}
		else
		{
			GRBVar x = model.addVar(0, 1, 0, GRB_BINARY);
			markerVars.pb(x);
		}
	}
}

void addMemberConstraint(GRBModel &model)
{

	for (auto it = familySet.begin(); it != familySet.end(); ++it)
	{
		GRBLinExpr refMembers = 0, tarMembers = 0;
		for (int i = 0; i < refMarkerSize; ++i)
		{
			if (jabs(refMarkers[i].family) == *it)
				refMembers += (1 - markerVars[i]);
		}
		for (int i = 0; i < tarMarkerSize; ++i)
		{
			if (jabs(tarMarkers[i].family) == *it)
				tarMembers += (1 - markerVars[i + refMarkerSize]);
		}
		model.addConstr(refMembers == tarMembers);	// same gene content
		model.addConstr(refMembers >= 1);			// intermediate
		//model.addConstr(refMembers == min(refFamSize[*it], tarFamSize[*it])); // maximum matching
		//model.addConstr(refMembers == 1);			// exemplar

		if (isdebug)
			cout << "add family " << *it << " constraint" << endl;
	}
}
void addHomoVar(GRBModel &model)
{
	string prevContig = "NotAContig";
	for (int i = 0; i < refMarkerSize; i++) //add ref capping marker
	{
		//telomeres
		//node on a new contig || node on the end of an old contig
		if (prevContig != refMarkers[i].contig || i == refMarkerSize - 1 || refMarkers[i].contig != refMarkers[i + 1].contig)
		{
			if (prevContig != refMarkers[i].contig && (i == refMarkerSize - 1 || refMarkers[i].contig != refMarkers[i + 1].contig))
			{
				prevContig = refMarkers[i].contig;
				//add two capping genes
				//refTelomeres.pb(i);
				Marker marker(refMarkerSize + refCapSize, capFamily, "CAP" + to_string(refCapSize++));
				if (isdebug)
					marker.show();
				refMarkers.pb(marker);
				//refTelomeres.pb(i);
				Marker marker2(refMarkerSize + refCapSize, capFamily, "CAP" + to_string(refCapSize++));
				if (isdebug)
					marker2.show();
				refMarkers.pb(marker2);
			}
			else
			{
				prevContig = refMarkers[i].contig;
				//refTelomeres.pb(i);
				//add a capping gene(marker)
				Marker marker(refMarkerSize + refCapSize, capFamily, "CAP" + to_string(refCapSize++));
				if (isdebug)
					marker.show();
				refMarkers.pb(marker);
			}
		}
	}
	prevContig = "NotAContig";
	for (int i = 0; i < tarMarkerSize; i++) //add tar capping marker
	{
		//telomeres
		//node on a new contig || node on the end of an old contig
		if (prevContig != tarMarkers[i].contig || i == tarMarkerSize - 1 || tarMarkers[i].contig != tarMarkers[i + 1].contig)
		{
			//single marker on a contig (add two cap)
			if (prevContig != tarMarkers[i].contig && (i == tarMarkerSize - 1 || tarMarkers[i].contig != tarMarkers[i + 1].contig))
			{
				prevContig = tarMarkers[i].contig;
				//tarTelomeres.pb(i);
				//add two capping genes(marker)
				Marker marker(tarMarkerSize + tarCapSize, capFamily, "CAP" + to_string(tarCapSize++));
				if (isdebug)
					marker.show();
				tarMarkers.pb(marker);
				//tarTelomeres.pb(i);
				Marker marker2(tarMarkerSize + tarCapSize, capFamily, "CAP" + to_string(tarCapSize++));
				if (isdebug)
					marker2.show();
				tarMarkers.pb(marker2);
			}
			else
			{
				prevContig = tarMarkers[i].contig;
				//tarTelomeres.pb(i);
				//add a capping gene(marker)
				Marker marker(tarMarkerSize + tarCapSize, capFamily, "CAP" + to_string(tarCapSize++));
				if (isdebug)
					marker.show();
				tarMarkers.pb(marker);
			}
		}
	}
	rc = refCapSize;
	tc = tarCapSize;
	if (refCapSize > tarCapSize)
	{
		while (refCapSize - tarCapSize) //1207 debug
		{
			Marker marker(tarMarkerSize + tarCapSize, capFamily, "CAP" + to_string(tarCapSize++));
			if (isdebug)
				marker.show();
			tarMarkers.pb(marker);
		}
	}
	else
	{
		while (tarCapSize - refCapSize) //1207 debug
		{
			Marker marker(refMarkerSize + refCapSize, capFamily, "CAP" + to_string(refCapSize++));
			if (isdebug)
				marker.show();
			refMarkers.pb(marker);
		}
	}
	for (int i = 0; i < refMarkerSize + refCapSize; i++) //find same family marker pairs
	{
		for (int j = 0; j < tarMarkerSize + tarCapSize; j++)
		{
			if (refMarkers[i].family == tarMarkers[j].family || -refMarkers[i].family == tarMarkers[j].family)
			{
				Edge edge(0, i, j, 2);
				if (isdebug)
					edge.show();
				if (i < refMarkerSize && j < tarMarkerSize)
				{
					//cout << i << " " << j << " " << refMarkers[i].family << " " << tarMarkers[j].family << endl;
					if (refFamSize[jabs(refMarkers[i].family)] == 1 && tarFamSize[jabs(tarMarkers[j].family)] == 1)
					{
						GRBVar x = model.addVar(1, 1, 0, GRB_BINARY);
						homoEdges.pb(edge);
						homoVars.pb(x);
					}
					else
					{
						GRBVar x = model.addVar(0, 1, 0, GRB_BINARY);
						homoEdges.pb(edge);
						homoVars.pb(x);
					}
				}
				else
				{
					GRBVar x = model.addVar(0, 1, 0, GRB_BINARY);
					homoEdges.pb(edge);
					homoVars.pb(x);
				}
			}
		}
	}
}
void addHomoConstraint(GRBModel &model)
{
	for (int i = 0; i < refMarkerSize; i++) //without cap
	{
		GRBLinExpr currNodeSum = 0;
		for (int j = 0; j < homoEdges.size(); j++)
		{
			if (homoEdges[j].node1 == i)
				currNodeSum += homoVars[j];
		}
		model.addConstr(currNodeSum == 1 - markerVars[i]);
	}
	for (int i = 0; i < tarMarkerSize; i++) //without cap
	{
		GRBLinExpr currNodeSum = 0;
		for (int j = 0; j < homoEdges.size(); j++)
		{
			if (homoEdges[j].node2 == i)
				currNodeSum += homoVars[j];
		}
		model.addConstr(currNodeSum == 1 - markerVars[i + refMarkerSize]);
	}

	for (int i = refMarkerSize; i < refMarkerSize + refCapSize; i++) //deal cap
	{
		GRBLinExpr currNodeSum = 0;
		for (int j = 0; j < homoEdges.size(); j++)
		{
			if (homoEdges[j].node1 == i)
				currNodeSum += homoVars[j];
		}
		model.addConstr(currNodeSum == 1);
	}
	for (int i = tarMarkerSize; i < tarMarkerSize + tarCapSize; i++) //deal cap
	{
		GRBLinExpr currNodeSum = 0;
		for (int j = 0; j < homoEdges.size(); j++)
		{
			if (homoEdges[j].node2 == i)
				currNodeSum += homoVars[j];
		}
		model.addConstr(currNodeSum == 1);
	}
}
void addNodeVar(GRBModel &model)
{
	int nodeid = 1;
	bool is_head;
	for (int i = 0; i < refMarkerSize; i++) //two node for a ref marker
	{
		if (refMarkers[i].family > 0)
			is_head = 0;
		else
			is_head = 1;

		Node node1(nodeid, refMarkers[i].id, refMarkers[i].family, refMarkers[i].contig, is_head);
		GRBVar y1 = model.addVar(1, nodeid++, 0, GRB_INTEGER);
		Node node2(nodeid, refMarkers[i].id, refMarkers[i].family, refMarkers[i].contig, !is_head);
		GRBVar y2 = model.addVar(1, nodeid++, 0, GRB_INTEGER);
		nodeVecs.pb(node1);
		nodeVars.pb(y1);
		nodeVecs.pb(node2);
		nodeVars.pb(y2);
	}
	int refCapMarkerID = refMarkers.size();
	for (int i = refMarkerSize; i < refMarkerSize + refCapSize; i++) //one node for cap
	{
		Node node(nodeid, refCapMarkerID++, capFamily, refMarkers[i].contig, 0);
		GRBVar y = model.addVar(1, nodeid++, 0, GRB_INTEGER);
		nodeVecs.pb(node);
		nodeVars.pb(y);
	}
	refNodeSize = nodeVars.size();
	for (int i = 0; i < tarMarkerSize; i++) //two node for a tar marker
	{
		if (tarMarkers[i].family > 0)
			is_head = 0;
		else
			is_head = 1;

		Node node1(nodeid, tarMarkers[i].id, tarMarkers[i].family, tarMarkers[i].contig, is_head);
		GRBVar y1 = model.addVar(1, nodeid++, 0, GRB_INTEGER);
		Node node2(nodeid, tarMarkers[i].id, tarMarkers[i].family, tarMarkers[i].contig, !is_head);
		GRBVar y2 = model.addVar(1, nodeid++, 0, GRB_INTEGER);
		nodeVecs.pb(node1);
		nodeVars.pb(y1);
		nodeVecs.pb(node2);
		nodeVars.pb(y2);
	}
	int tarCapMarkerID = tarMarkers.size();
	for (int i = tarMarkerSize; i < tarMarkerSize + tarCapSize; i++) //one node for cap
	{
		Node node(nodeid, tarCapMarkerID++, capFamily, tarMarkers[i].contig, 0);
		GRBVar y = model.addVar(1, nodeid++, 0, GRB_INTEGER);
		nodeVecs.pb(node);
		nodeVars.pb(y);
	}
}
void addAdjCheck()
{
	for (int i = 0; i < refMarkerSize; i++)
	{
		for (int j = i + 1; j < refMarkerSize; j++)
		{
			if (refMarkers[i].contig != refMarkers[j].contig)
			{
				int fam1 = jabs(refMarkers[i].family), fam2 = jabs(refMarkers[j].family);
				if (fam1 > fam2)
					swap(fam1, fam2);
				pair<int, int> mypair(fam1, fam2);
				if (tarPotentialAdj[mypair] == 1)
				{
					pair<string, string> spair(refMarkers[i].contig, refMarkers[j].contig);
					refAdj[spair] = 1;
				}
			}
		}
	}
	for (int i = 0; i < tarMarkerSize; i++)
	{
		for (int j = i + 1; j < tarMarkerSize; j++)
		{
			if (tarMarkers[i].contig != tarMarkers[j].contig)
			{
				int fam1 = jabs(tarMarkers[i].family), fam2 = jabs(tarMarkers[j].family);
				if (fam1 > fam2)
					swap(fam1, fam2);
				pair<int, int> mypair(fam1, fam2);
				if (refPotentialAdj[mypair] == 1)
				{
					pair<string, string> spair(tarMarkers[i].contig, tarMarkers[j].contig);
					tarAdj[spair] = 1;
				}
			}
		}
	}
}
void addAdjVar(GRBModel &model)
{
	//ref
	for (int i = 0; i < 2 * refMarkerSize - 1; i++) //true adj
	{
		if (nodeVecs[i].markerid != nodeVecs[i + 1].markerid && nodeVecs[i].contig == nodeVecs[i + 1].contig)
		{
			Edge edge(0, i, i + 1, 0);
			if (isdebugadj)
				edge.show();
			GRBVar PA = model.addVar(1, 1, 0, GRB_BINARY);
			adjEdges.pb(edge);
			adjVars.pb(PA);
		}
	}
	string prevContig = "NotAContig";
	int currCapMarkerID = 2 * refMarkerSize;
	for (int i = 0; i < 2 * refMarkerSize; i++) //telo with ori cap
	{
		if (prevContig != nodeVecs[i].contig || i == 2 * refMarkerSize - 1 || nodeVecs[i].contig != nodeVecs[i + 1].contig)
		{
			prevContig = nodeVecs[i].contig;
			refTelomeres.pb(i);
			Edge edge(1, i, currCapMarkerID++, 0);
			if (isdebugadj)
				edge.show();
			GRBVar PA = model.addVar(0, 1, 0, GRB_BINARY);
			adjEdges.pb(edge);
			//teloWithCapSum += PA; //0103 mod
			adjVars.pb(PA);
			refPotAdjNumt2c++;
		}
	}

	for (int i = 0; i < refTelomeres.size(); i++) //telo with telo
	{
		for (int j = i + 1; j < refTelomeres.size(); j++)
		{
			if (nodeVecs[refTelomeres[i]].contig == nodeVecs[refTelomeres[j]].contig) //telos on the same contig
				continue;
			pair<string, string> mypair(nodeVecs[refTelomeres[i]].contig, nodeVecs[refTelomeres[j]].contig);
			if (refAdj[mypair] != 1)
				continue;
			Edge edge(1, refTelomeres[i], refTelomeres[j], 0);
			if (isdebugadj)
				edge.show();
			GRBVar PA = model.addVar(0, 1, 0, GRB_BINARY);
			adjEdges.pb(edge);
			adjVars.pb(PA);
			refPotAdjNumt2t++;
		}
	}
	for (int i = 2 * refMarkerSize; i < 2 * refMarkerSize + rc; i++) //cap with cap
	{
		for (int j = i + 1; j < 2 * refMarkerSize + rc; j++)
		{
			Edge edge(1, i, j, 0);
			if (isdebugadj)
				edge.show();
			GRBVar PA = model.addVar(0, 1, 0, GRB_BINARY);
			adjEdges.pb(edge);
			adjVars.pb(PA);
			refPotAdjNumc2c++;
		}
	}
	if(tc > rc)
	{
		for (int i = 2 * refMarkerSize + rc; i < 2 * refMarkerSize + refCapSize; i+=2) //true cap with cap
		{
			Edge edge(0, i, i+1, 0);
			if (isdebugadj)
				edge.show();
			GRBVar PA = model.addVar(1, 1, 0, GRB_BINARY);
			adjEdges.pb(edge);
			adjVars.pb(PA);
			refRealAdjNumc2c++;
		}
	}
	/***********************************************************/
	/*
	for (int i = 2 * refMarkerSize; i < 2 * refMarkerSize + refCapSize; i++) //true cap with cap
	{
		for (int j = i + 1; j < 2 * refMarkerSize + refCapSize; j++)
		{
			Edge edge(0, i, j, 0);
			if (isdebugadj)
				edge.show();
			GRBVar PA = model.addVar(0, 1, 0, GRB_BINARY);
			adjEdges.pb(edge);
			adjVars.pb(PA);
		}
	}
	/***********************************************************/

	refAdjSize = adjEdges.size();
	//tar
	int offset = 2 * refMarkerSize + refCapSize;
	for (int i = 0; i < 2 * tarMarkerSize - 1; i++) //true adj
	{
		if (nodeVecs[offset + i].markerid != nodeVecs[offset + i + 1].markerid && nodeVecs[offset + i].contig == nodeVecs[offset + i + 1].contig)
		{
			Edge edge(0, i, i + 1, 0);
			if (isdebugadj)
				edge.show();
			GRBVar PA = model.addVar(1, 1, 0, GRB_BINARY);
			adjEdges.pb(edge);
			adjVars.pb(PA);
		}
	}
	prevContig = "NotAContig";
	currCapMarkerID = 2 * tarMarkerSize;
	for (int i = 0; i < 2 * tarMarkerSize; i++) //telo with ori cap
	{
		if (prevContig != nodeVecs[offset + i].contig || i == 2 * tarMarkerSize - 1 || nodeVecs[offset + i].contig != nodeVecs[offset + i + 1].contig)
		{
			prevContig = nodeVecs[offset + i].contig;
			tarTelomeres.pb(i);
			Edge edge(0, i, currCapMarkerID++, 0);
			if (isdebugadj)
				edge.show();
			GRBVar PA = model.addVar(0, 1, 0, GRB_BINARY);
			adjEdges.pb(edge);
			//teloWithCapSum += PA; //0103 mod
			adjVars.pb(PA);
			tarPotAdjNumt2c++;
		}
	}
	for (int i = 0; i < tarTelomeres.size(); i++) //telo with telo
	{
		for (int j = i + 1; j < tarTelomeres.size(); j++)
		{
			if (nodeVecs[offset + tarTelomeres[i]].contig == nodeVecs[offset + tarTelomeres[j]].contig) //telos on the same contig
				continue;
			pair<string, string> mypair(nodeVecs[offset + tarTelomeres[i]].contig, nodeVecs[offset + tarTelomeres[j]].contig);
			if (tarAdj[mypair] != 1)
				continue;
			Edge edge(1, tarTelomeres[i], tarTelomeres[j], 0);
			if (isdebugadj)
				edge.show();
			GRBVar PA = model.addVar(0, 1, 0, GRB_BINARY);
			adjEdges.pb(edge);
			adjVars.pb(PA);
			tarPotAdjNumt2t++;
		}
	}
	for (int i = 2 * tarMarkerSize; i < 2 * tarMarkerSize + tc; i++) //cap with cap
	{
		for (int j = i + 1; j < 2 * tarMarkerSize + tc; j++)
		{
			Edge edge(1, i, j, 0);
			if (isdebugadj)
				edge.show();
			GRBVar PA = model.addVar(0, 1, 0, GRB_BINARY);
			adjEdges.pb(edge);
			adjVars.pb(PA);
			tarPotAdjNumc2c++;
		}
	}
	if(rc > tc)
	{
		for (int i = 2 * tarMarkerSize + tc; i < 2 * tarMarkerSize + tarCapSize; i+=2) //true cap with cap
		{
			Edge edge(0, i, i+1, 0);
			if (isdebugadj)
				edge.show();
			GRBVar PA = model.addVar(1, 1, 0, GRB_BINARY);
			adjEdges.pb(edge);
			adjVars.pb(PA);
			tarRealAdjNumc2c++;
		}
	}
	tarAdjSize = adjEdges.size() - refAdjSize;
}
void addAdjConstraint(GRBModel &model)
{
	for (int i = 0; i < 2 * refMarkerSize + refCapSize; i++)
	{
		GRBLinExpr currNodeEdgeSum = 0;
		for (int j = 0; j < refAdjSize; j++)
		{
			if (adjEdges[j].node1 == i || adjEdges[j].node2 == i)
			{
				currNodeEdgeSum += adjVars[j];
			}
		}
		model.addConstr(currNodeEdgeSum == 1);
	}
	for (int i = 0; i < 2 * tarMarkerSize + tarCapSize; i++)
	{
		GRBLinExpr currNodeEdgeSum = 0;
		for (int j = 0; j < tarAdjSize; j++)
		{
			if (adjEdges[refAdjSize + j].node1 == i || adjEdges[refAdjSize + j].node2 == i)
			{
				currNodeEdgeSum += adjVars[refAdjSize + j]; //1207 find bug
			}
		}
		model.addConstr(currNodeEdgeSum == 1);
	}
}
void addCycleConstraint(GRBModel &model)
{
	for (int i = 0; i < homoEdges.size(); i++) //make nodes of homo edge have same labels
	{
		Edge edge = homoEdges[i];
		int ref = edge.node1, tar = edge.node2;
		int refhead, tarhead, reftail, tartail;
		int offset = 2 * refMarkerSize + refCapSize;
		if (ref < refMarkerSize) //not cap
		{

			if (nodeVecs[2 * ref].is_head)
			{
				refhead = 2 * ref;
				reftail = 2 * ref + 1;
			}
			else
			{
				refhead = 2 * ref + 1;
				reftail = 2 * ref;
			}
			if (nodeVecs[offset + 2 * tar].is_head)
			{
				tarhead = offset + 2 * tar;
				tartail = offset + 2 * tar + 1;
			}
			else
			{
				tarhead = offset + 2 * tar + 1;
				tartail = offset + 2 * tar;
			}
			if (isdebugall)
			{
				cout << "homo node " << refhead << " <= node " << tarhead << endl;
				cout << "homo node " << reftail << " <= node " << tartail << endl;
			}
			model.addConstr(nodeVars[refhead] <= (nodeVars[tarhead] + (1 - homoVars[i]) * (refhead + 1)));
			model.addConstr(nodeVars[reftail] <= (nodeVars[tartail] + (1 - homoVars[i]) * (reftail + 1)));
			model.addConstr(nodeVars[tarhead] <= (nodeVars[refhead] + (1 - homoVars[i]) * (tarhead + 1)));
			model.addConstr(nodeVars[tartail] <= (nodeVars[reftail] + (1 - homoVars[i]) * (tartail + 1))); //1207 debug
		}
		else //cap
		{
			if (isdebugall)
			{
				cout << "homo node " << ref + refMarkerSize << " <= node " << offset + tar + tarMarkerSize << endl;
			}

			model.addConstr(nodeVars[ref + refMarkerSize] <= (nodeVars[offset + tar + tarMarkerSize] + (1 - homoVars[i]) * (ref + refMarkerSize + 1)));
			model.addConstr(nodeVars[offset + tar + tarMarkerSize] <= (nodeVars[ref + refMarkerSize] + (1 - homoVars[i]) * (offset + tar + tarMarkerSize + 1)));
		}
	}
	for (int i = 0; i < refAdjSize; i++) //ref adj
	{
		Edge edge = adjEdges[i];
		int ref1 = edge.node1, ref2 = edge.node2;
		if (isdebugall)
		{
			cout << "ref adj node " << ref1 << " <= node " << ref2 << endl;
		}
		model.addConstr(nodeVars[ref1] <= (nodeVars[ref2] + (1 - adjVars[i]) * (ref1 + 1)));
		model.addConstr(nodeVars[ref2] <= (nodeVars[ref1] + (1 - adjVars[i]) * (ref2 + 1)));
	}
	int offset = 2 * refMarkerSize + refCapSize;
	for (int i = refAdjSize; i < adjEdges.size(); i++) //tar adj
	{
		Edge edge = adjEdges[i];
		int tar1 = edge.node1, tar2 = edge.node2;
		if (isdebugall)
		{
			cout << "tar adj node " << offset + tar1 << " <= node " << offset + tar2 << endl;
		}

		model.addConstr(nodeVars[offset + tar1] <= (nodeVars[offset + tar2] + (1 - adjVars[i]) * (offset + tar1 + 1)));
		model.addConstr(nodeVars[offset + tar2] <= (nodeVars[offset + tar1] + (1 - adjVars[i]) * (offset + tar2 + 1)));
	}
	return;
	for (int i = 0; i < markerVars.size(); i++)	// int adj same label
	{
		int node1 = 2*i   + i < refMarkerSize ? 0 : rc ;
		int node2 = 2*i+1 + i < refMarkerSize ? 0 : rc ;
		if (isdebugall)
		{
			cout << "int adj: node " << node1 << " <= node " << node2 << endl;
		}
		model.addConstr(nodeVars[node1] <= (nodeVars[node2] + (1 - markerVars[i]) * (node1 + 1)));
		model.addConstr(nodeVars[node2] <= (nodeVars[node1] + (1 - markerVars[i]) * (node2 + 1)));
	}
}
void addInternalConstraint(GRBModel &model)
{
	for (int i = 0; i < refMarkerSize; i++)
	{
		if (isdebugall)
		{
			cout << "internal node " << 2 * i << " <= node " << 2 * i + 1 << endl;
		}
		model.addConstr(nodeVars[2 * i] <= (nodeVars[2 * i + 1] + (1 - markerVars[i]) * (2 * i + 1)));
		model.addConstr(nodeVars[2 * i + 1] <= (nodeVars[2 * i] + (1 - markerVars[i]) * (2 * i + 1 + 1)));
	}
	int offset = 2 * refMarkerSize + refCapSize;
	for (int i = 0; i < tarMarkerSize; i++)
	{
		if (isdebugall)
		{
			cout << "internal node " << offset + 2 * i << " <= node " << offset + 2 * i + 1 << endl;
		}
		model.addConstr(nodeVars[offset + 2 * i] <= (nodeVars[offset + 2 * i + 1] + (1 - markerVars[refMarkerSize + i]) * (offset + 2 * i + 1)));
		model.addConstr(nodeVars[offset + 2 * i + 1] <= (nodeVars[offset + 2 * i] + (1 - markerVars[refMarkerSize + i]) * (offset + 2 * i + 1 + 1)));
	}
}
void addCycleVarAndConstraint(GRBModel &model)
{
	//fix cycle of length 2
	/*
	int offset = 2 * refMarkerSize + refCapSize;
	for (int i = 0; i < refMarkerSize - 1; i++)
	{
		if (refFamSize[jabs(refMarkers[i].family)] == 1 && refFamSize[jabs(refMarkers[i + 1].family)] == 1) //singleton on ref
		{
			if (refMarkers[i].contig == refMarkers[i + 1].contig) //iso on ref singlton's right
			{
				for (int j = 0; j < tarMarkerSize - 1; j++)
				{
					if (tarMarkers[j].family == refMarkers[i].family && tarFamSize[jabs(tarMarkers[j].family)] == 1 && tarFamSize[jabs(tarMarkers[j + 1].family)] == 1) //shared + singleton on tar
					{
						//cout<<"ref"<<i<<" tar"<<j<<endl;
						if (tarMarkers[j + 1].family == refMarkers[i + 1].family && tarMarkers[j].contig == tarMarkers[j + 1].contig) //shared + iso on tar
						{
							model.addConstr(nodeVars[i * 2 + 1] == i * 2 + 2);
							model.addConstr(nodeVars[i * 2 + 2] == i * 2 + 2);
							model.addConstr(nodeVars[offset + j * 2 + 1] == i * 2 + 2);
							model.addConstr(nodeVars[offset + j * 2 + 2] == i * 2 + 2);
						}
					}
					//neg
					else if (tarMarkers[j].family == -refMarkers[i].family && tarFamSize[jabs(tarMarkers[j].family)] == 1 && j - 1 >= 0 && tarFamSize[jabs(tarMarkers[j - 1].family)] == 1) //singleton on tar
					{
						if (tarMarkers[j - 1].family == -refMarkers[i + 1].family && tarMarkers[j].contig == tarMarkers[j - 1].contig) //shared + iso on tar
						{
							model.addConstr(nodeVars[i * 2 + 1] == i * 2 + 2);
							model.addConstr(nodeVars[i * 2 + 2] == i * 2 + 2);
							model.addConstr(nodeVars[offset + j * 2] == i * 2 + 2);
							model.addConstr(nodeVars[offset + j * 2 - 1] == i * 2 + 2);
						}
					}
				}
			}
			if (i - 1 >= 0 && refMarkers[i].contig == refMarkers[i - 1].contig) //iso on ref singleton's left
			{
				for (int j = 0; j < tarMarkerSize - 1; j++)
				{
					if (tarMarkers[j].family == refMarkers[i].family && tarFamSize[jabs(tarMarkers[j].family)] == 1 && j - 1 >= 0 && tarFamSize[jabs(tarMarkers[j - 1].family)] == 1) //shared + singleton on tar
					{
						//cout<<"p2ref"<<i<<" tar"<<j<<endl;
						//cout<<tarIso[j-1]<<endl;
						if (tarMarkers[j - 1].family == refMarkers[i - 1].family && tarMarkers[j].contig == tarMarkers[j - 1].contig) //shared + iso on tar
						{
							model.addConstr(nodeVars[i * 2] == i * 2);
							model.addConstr(nodeVars[i * 2 - 1] == i * 2);
							model.addConstr(nodeVars[offset + j * 2] == i * 2);
							model.addConstr(nodeVars[offset + j * 2 - 1] == i * 2);
						}
					}
					//neg
					else if (tarMarkers[j].family == -refMarkers[i].family && tarFamSize[jabs(tarMarkers[j].family)] == 1) //singleton on tar
					{
						if (tarMarkers[j + 1].family == -refMarkers[i - 1].family && tarMarkers[j].contig == tarMarkers[j + 1].contig) //shared + iso on tar
						{
							model.addConstr(nodeVars[i * 2] == i * 2);
							model.addConstr(nodeVars[i * 2 - 1] == i * 2);
							model.addConstr(nodeVars[offset + j * 2 + 1] == i * 2);
							model.addConstr(nodeVars[offset + j * 2 + 2] == i * 2);
						}
					}
				}
			}
		}
	}
	*/
	for (int i = 0; i < nodeVars.size(); i++) //using refNodeSize instead of nodeVars.size() since upper bound of ref < tar
	{
		GRBVar z = model.addVar(0, 1, 0, GRB_BINARY);
		cycleVars.pb(z);
		model.addConstr(z * (i + 1) <= nodeVars[i]);
		model.addConstr(1 - z <= (i + 1) - nodeVars[i]);
	}
}
void findMax(GRBModel &model)
{
	GRBLinExpr cycleSum = 0;
	for (int i = 0; i < cycleVars.size(); i++) //using refNodeSize instead of cycleVars.size() since upper bound of ref < tar
		cycleSum += cycleVars[i];
	//model.setObjective(2 * refRemoveMarkerNum - cycleSum); //removal cost = 3
	model.setObjective(-refRemoveMarkerNum - cycleSum); //no removal cost
	//model.setObjective(-cycleSum);
	model.getEnv().set("TimeLimit", "1800");

	//model.setObjectiveN(-refRemoveMarkerNum - cycleSum, 0, 1);
	//model.setObjectiveN(-teloWithCapSum,1,0);// teloWithCapSum statements had been commented out

	//auto env0 = model.getMultiobjEnv(0);
	//auto env1 = model.getMultiobjEnv(1);

	/* shao's config*/
	//env0.set(GRB_DoubleParam_Heuristics, 0.5);
	//env0.set(GRB_IntParam_MIPFocus, 1);

	//env0.set("TimeLimit", "3600");
	//env1.set("TimeLimit","3600");
	model.update();
	//model.write(out_dir+"/debug.lp");
	model.optimize();
}

void showResult()
{
	int cycles = 0;
	for (int i = 0; i < cycleVars.size(); i++)
		if (cycleVars[i].get(GRB_DoubleAttr_X) == 1)
			cycles++;
	cout << "Cycles: " << cycles << endl;
	if (isdebug)
	{
		for (int i = 0; i < markerVars.size(); i++)
		{
			cout << "marker " << (i < refMarkerSize ? refMarkers[i].family : tarMarkers[i-refMarkerSize].family) << ": " << markerVars[i].get(GRB_DoubleAttr_X) << endl;
		}
		cout << "-----------------------Homo---------------------" << endl;
		for (int i = 0; i < homoEdges.size(); i++)
		{
			if (homoVars[i].get(GRB_DoubleAttr_X) == 1)
				cout << "homoEdges " << homoEdges[i].node1 << " " << homoEdges[i].node2 << endl;
		}
		cout << "-----------------------Adj---------------------" << endl;
		for (int i = 0; i < adjEdges.size(); i++)
		{
			if (adjVars[i].get(GRB_DoubleAttr_X) == 1)
				cout << "adjEdges " << adjEdges[i].node1 << " " << adjEdges[i].node2 << endl;
		}
		for (int i = 0; i < nodeVars.size(); i++)
		{
			cout << "Cycle" << i << ": " << nodeVars[i].get(GRB_DoubleAttr_X) << endl;
		}
	}
	int offset = 2 * refMarkerSize + refCapSize;
	ofstream fout(out_dir+"/joins.txt");
	for (int i = refAdjSize; i < adjEdges.size(); i++)
	{
		if (adjVars[i].get(GRB_DoubleAttr_X) == 1 && adjEdges[i].is_potential == 1)
		{
			int sign1 = -1, sign2 = -1;
			//cout<<"adjEdges "<<adjEdges[i].node1<<" "<<adjEdges[i].node2<<endl;
			//nodeVecs[adjEdges[i].node1+offset].showSpaceVersion();
			//nodeVecs[adjEdges[i].node2+offset].showSpaceVersion();
			if (nodeVecs[adjEdges[i].node1 + offset].is_head)
				sign1 = 1;
			if (nodeVecs[adjEdges[i].node2 + offset].is_head)
				sign2 = 1;
			if (jabs(nodeVecs[adjEdges[i].node1 + offset].markerid) <= tarMarkerSize && jabs(nodeVecs[adjEdges[i].node2 + offset].markerid) <= tarMarkerSize)
				fout << "marker " << sign1 * nodeVecs[adjEdges[i].node1 + offset].markerid << " " << sign2 * nodeVecs[adjEdges[i].node2 + offset].markerid << endl;
		}
	}
	fout.close();
	fout.open(out_dir+"/duplicated.txt");
	int refRemoveMarkerNums = 0;
	for (int i = 0; i < tarMarkerSize; i++)
	{
		if (markerVars[refMarkerSize + i].get(GRB_DoubleAttr_X) == 1)
		{
			fout << tarMarkers[i].contig << endl;
			refRemoveMarkerNums++;
		}
	}
	// cout << "refMarkerSize: " << refMarkerSize << " " << refMarkerNum << endl;
	// cout << "refRemoveMarkerNum: " << refRemoveMarkerNum.getValue() << endl;
	// cout << "tarContigNum: " << tarContigNum << endl;
	// cout << "refContigNum: " << refContigNum << endl;
	//cout<<"----------------------node---------------------"<<endl;
}
void showVars()
{
	cout << "markerVar size: " << markerVars.size() << endl;
	cout << "homoVar size: " << homoVars.size() << endl;
	cout << "cycleVar size: " << cycleVars.size() << endl;
	cout << "nodeVar size: " << nodeVars.size() << endl;
	cout << "adjVar size: " << adjVars.size() << endl;
}
void showDebugInfo(bool isdebug)
{
	if (!isdebugall)
		return;
	cout << endl
		 << endl;
	cout << "------------------Homo-------------------" << endl;
	for (int i = 0; i < homoEdges.size(); i++)
		homoEdges[i].show();
	cout << "------------------Adj--------------------" << endl;
	for (int i = 0; i < adjEdges.size(); i++)
		adjEdges[i].show();
	cout << "--------------Ref Marker-----------------" << endl;
	for (int i = 0; i < refMarkers.size(); i++)
		refMarkers[i].show();
	cout << "--------------Tar Marker-----------------" << endl;
	for (int i = 0; i < tarMarkers.size(); i++)
		tarMarkers[i].show();
	cout << "--------------Ref Node-------------------" << endl;
	for (int i = 0; i < 2 * refMarkerSize + refCapSize; i++)
	{
		cout << "Node " << i << " :";
		nodeVecs[i].show();
	}
	cout << "--------------Tar Node-------------------" << endl;
	for (int i = 2 * refMarkerSize + refCapSize; i < nodeVecs.size(); i++)
	{
		cout << "Node " << i << " :";
		nodeVecs[i].show();
	}
	cout << "------------------------------------------" << endl
		 << endl;
}
void graphDebug(bool debug)
{
	if (!debug)
		return;
	printf("===== ===== ===== graphDebug ===== ===== =====\n");
	printf("refMarkerSize: %d, refCapSize: %d\n", refMarkerSize, refCapSize);
	printf("tarMarkerSize: %d, tarCapSize: %d\n", tarMarkerSize, tarCapSize);
	// ref
	for (int i = 0; i < refMarkerSize*2; i++) {
		Marker& temp = refMarkers[i/2];
		printf("[%d] marker: (%d, %d, %s), removed: %d, label: %d\n", i, temp.id, temp.family, temp.contig.c_str(), (int)markerVars[i/2].get(GRB_DoubleAttr_X), (int)nodeVars[i].get(GRB_DoubleAttr_X));
	}
	// ref cap
	for (int i = refMarkerSize*2; i < refMarkerSize*2+refCapSize; i++) {
		Marker& temp = refMarkers[i-refMarkerSize];
		printf("[%d] marker: (%d, %d, %s), removed: %d, label: %d\n", i, temp.id, temp.family, temp.contig.c_str(), 999, (int)nodeVars[i].get(GRB_DoubleAttr_X));
	}
	// tar
	int offset = refMarkerSize*2 + refCapSize;
	for (int i = 0; i < tarMarkerSize*2; i++) {
		Marker& temp = tarMarkers[i/2];
		printf("[%d] marker: (%d, %d, %s), removed: %d, label: %d\n", i, temp.id, temp.family, temp.contig.c_str(), (int)markerVars[i/2+refMarkerSize].get(GRB_DoubleAttr_X), (int)nodeVars[i+offset].get(GRB_DoubleAttr_X));
	}
	// tar cap
	for (int i = tarMarkerSize*2; i < tarMarkerSize*2+tarCapSize; i++) {
		Marker& temp = tarMarkers[i-tarMarkerSize];
		printf("[%d] marker: (%d, %d, %s), removed: %d, label: %d\n", i, temp.id, temp.family, temp.contig.c_str(), 999, (int)nodeVars[i+offset].get(GRB_DoubleAttr_X));
	}
	printf("===== ===== ===== homoEdges ===== ===== =====\n");
	for (int i = 0; i < homoEdges.size(); i++) {
		Edge& temp = homoEdges[i];
		printf("edge: (%d, %d), kept: %d\n", temp.node1, temp.node2, (int)homoVars[i].get(GRB_DoubleAttr_X));
	}
	printf("===== ===== ===== adjEdges ===== ===== =====\n");
	for (int i = 0; i < adjEdges.size(); i++) {
		Edge& temp = adjEdges[i];
		printf("edge: (%d, %d), potential: %d, kept: %d\n", temp.node1, temp.node2, temp.is_potential, (int)adjVars[i].get(GRB_DoubleAttr_X));
	}
}
int main(int argc, char *argv[])
{
	out_dir = string(argc < 4 ? "output" : argv[3]);
	try
	{
		GRBEnv env = GRBEnv(true);
		env.set("LogToConsole", "0");
		env.set("LogFile", out_dir+"/gurobi.log");
		env.start();
		GRBModel model = GRBModel(env);
		//model.set("TimeLimit","120");
		getInput(argv[1], argv[2]);
		addDuplicatedVar(model);	//todo: no need to determine singletons
		addMemberConstraint(model); //1.same gene content 2.at least one in each fam
		if (isdebug)
			cout << "p1" << endl;
		addHomoVar(model);
		addHomoConstraint(model);
		addNodeVar(model);
		addAdjCheck();
		addAdjVar(model);
		if (isdebug)
			cout << "p2" << endl;
		//20210515 code review pass
		addAdjConstraint(model);
		addCycleConstraint(model);
		addInternalConstraint(model);
		addCycleVarAndConstraint(model);
		if (isdebug)
			cout << "p3" << endl;
		showDebugInfo(isdebug);
		showVars();
		// cout << "refRealAdjNumc2c: " << refRealAdjNumc2c << endl;
		// cout << "refPotentialAdjNumt2t: " << refPotAdjNumt2t << endl;
		// cout << "refPotentialAdjNumt2c: " << refPotAdjNumt2c << endl;
		// cout << "refPotentialAdjNumc2c: " << refPotAdjNumc2c << endl;
		// cout << "tarRealAdjNumc2c: " << tarRealAdjNumc2c << endl;
		// cout << "tarPotentialAdjNumt2t: " << tarPotAdjNumt2t << endl;
		// cout << "tarPotentialAdjNumt2c: " << tarPotAdjNumt2c << endl;
		// cout << "tarPotentialAdjNumc2c: " << tarPotAdjNumc2c << endl;
		findMax(model);
		showResult();
		graphDebug(true);
	}
	catch (GRBException e)
	{
		cout << "Error code = " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
	}
	catch (...)
	{
		cout << "Exception during optimization" << endl;
	}
	return 0;
}
