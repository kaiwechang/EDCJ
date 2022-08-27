#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>
#include <string>
#include <map>
#include <utility>

#define pb push_back
using namespace std;

struct Telomeres
{
	int lhs, rhs;
};
struct Edge
{
	int node1, node2;
};
struct Marker
{
	int id, family;
	string contig;
};

map<string, Telomeres> contigMap;
map<string, int> visited;
map<string, string> scontig;
int joinSide; //0 for lhs, 1 rhs
vector<Edge> edges;
vector<vector<pair<string, int>>> scaffolds;
vector<pair<string, int>> currScaffold;
string findContig(int node)
{
	for (auto c : contigMap)
	{
		if (c.second.lhs == node)
		{
			joinSide = 0;
			return c.first;
		}
		if (c.second.rhs == node)
		{
			joinSide = 1;
			return c.first;
		}
	}
	cout << "NOTFOUND!!";
}
void traverse(int node, bool isright)
{
	string nextContig;
	for (int i = 0; i < edges.size(); i++)
	{
		if (node == edges[i].node1 || node == edges[i].node2)
		{
			int nextNode = edges[i].node1;
			pair<string, int> ContigWithDir;
			if (node == nextNode)
				nextNode = edges[i].node2;
			//cout<<"NN:"<<nextNode<<endl;
			nextContig = findContig(nextNode);
			//cout<<"NC:"<<nextContig<<endl;
			if (visited[nextContig] == 1)
				return;
			visited[nextContig] = 1;
			ContigWithDir.first = nextContig;
			if (joinSide == isright) //join two contigs with same side => one is reversed
			{
				ContigWithDir.second = 1;
			}
			else
			{
				ContigWithDir.second = 0;
			}
			//cout<<ContigWithDir.first<<" "<<ContigWithDir.second<<endl;
			currScaffold.push_back(ContigWithDir);
			if (joinSide == 0)
			{
				//cout<<"tra("<<contigMap[nextContig].rhs<<",1)"<<endl;
				traverse(contigMap[nextContig].rhs, 1);
			}
			else
			{
				//cout<<"tra("<<contigMap[nextContig].lhs<<",0)"<<endl;
				traverse(contigMap[nextContig].lhs, 0);
			}
			break;
		}
	}
}
void scaffold(pair<string, Telomeres> contig)
{
	if (visited[contig.first] == 1)
		return;
	//cout << "start" << contig.first << endl;
	visited[contig.first] = 1;

	currScaffold.clear();
	//cout<<"tr"<<endl;
	pair<string, int> firstContig(contig.first, 0);
	//cout<<contig.first<<" "<<0<<endl;
	currScaffold.push_back(firstContig);

	traverse(contig.second.lhs, 0);
	//cout<<"tr2"<<endl;
	reverse(currScaffold.begin(), currScaffold.end());
	traverse(contig.second.rhs, 1);
	//cout<<"tr3"<<endl;
	scaffolds.push_back(currScaffold);
}
int main(int argc, char *argv[])
{
	string target_dir(argc < 4 ? "output" : argv[3]);

	ifstream fin;
	fin.open(argv[1]);

	int id, family, tmp;
	string contig, pcontig="";
	int count = 0; 
	vector<Marker> markers;
	map<string, int> oriMarkerNum, removedMarkerNum;
	while (fin >> id >> family >> contig >> tmp)
	{
		if(contig != pcontig)
		{
			count++;
			scontig[to_string(count)] = contig;
		}
		pcontig = contig;
		Marker m;
		m.id = id;
		m.family = family;
		m.contig = to_string(count);
		markers.push_back(m);
		oriMarkerNum[contig]++;
	}
	fin.close();
	fin.open(target_dir+"/duplicated.txt");
    while (fin >> contig)
        removedMarkerNum[contig]++;
    fin.close();
	string prevContig = "NotAContig";
	Telomeres t;
	for (int i = 0; i < markers.size() - 1; i++)
	{
		Marker m;
		m = markers[i];
		if (m.contig != prevContig) //begin of a contig
		{
			prevContig = m.contig;
			int startIsHead = -1;
			if (m.family < 0)
				startIsHead = 1;
			t.lhs = m.id * startIsHead;
		}
		if (m.contig != markers[i + 1].contig) //end of a contig
		{
			int endIsHead = -1;
			if (m.family > 0)
				endIsHead = 1;
			t.rhs = m.id * endIsHead;
			contigMap[m.contig] = t;
			visited[m.contig] = 0;
		}
	}
	Marker m;
	m = markers[markers.size() - 1];
	if (m.contig == prevContig)
	{
		int endIsHead = -1;
		if (m.family > 0)
			endIsHead = 1;
		t.rhs = m.id * endIsHead;
		contigMap[m.contig] = t;
		visited[m.contig] = 0;
	}
	else
	{
		int startIsHead = -1;
		if (m.family < 0)
			startIsHead = 1;
		t.lhs = m.id * startIsHead;
		t.rhs = -t.lhs;
		contigMap[m.contig] = t;
		visited[m.contig] = 0;
	}
	// for(auto c:contigMap)
	// {
	// 	cout<<c.first<<" "<<c.second.lhs<<" "<<c.second.rhs<<endl;
		
	// }

	string joinFile = argv[2];
	ifstream fin2(joinFile);
	fin.open(argv[2]);
	//fin2.open(joinFile);
	string stmp;
	int node1, node2;
	while (fin >> stmp >> node1 >> node2)
	{
		Edge edge;
		edge.node1 = node1;
		edge.node2 = node2;
		edges.push_back(edge);
	}
	fin.close();
	for (auto c : contigMap)
	{
		scaffold(c);
		//cout<<"g0"<<endl;
	}

	
	for (auto currScaffold = scaffolds.begin(); currScaffold != scaffolds.end(); currScaffold++)
	{
		for (auto contig = (*currScaffold).begin(); contig != (*currScaffold).end(); contig++)
		{
			(*contig).first = scontig[(*contig).first];
		}
	}
	// for (auto currScaffold : scaffolds)
	// {
	// 	for (auto contig : currScaffold)
	// 	{
	// 		cout << contig.first << " " << contig.second << endl;
	// 	}
	// }
	// cout << "finish" << endl;
	bool sign = 0;
	for (auto currScaffold = scaffolds.begin(); currScaffold != scaffolds.end(); currScaffold++)
	{
		for (auto contig = (*currScaffold).begin(); contig != (*currScaffold).end(); contig++)
		{
			if ((*contig).second == 1)
				sign = !sign;
			(*contig).second = sign;
		}
	}

//	for (auto currScaffold : scaffolds)
//	{
//		for (auto contig : currScaffold)
//		{
//			cout << contig.first << " " << contig.second << endl;
//		}
//	}
	for (auto it = scaffolds.begin(); it != scaffolds.end(); it++)
	{
		for (auto it2 = (*it).begin(); it2 != (*it).end(); it2++)
		{
			vector<vector<pair<string, int>>> scaffolds;
			if(oriMarkerNum[(*it2).first] == removedMarkerNum[(*it2).first])
			{
				(*it).erase(it2);
				it2--;
			}
		}
	} 
	map<string, vector<pair<string, int>>> tar;
	fin.open(target_dir+"/tarContigMerge.all");
	string key;
	while (fin >> key)
	{
		vector<pair<string, int>> temp;
		while(fin >> contig)
		{
			if(contig == "end")
			{
				break;
			}
			fin >> tmp;
			temp.pb(make_pair(contig, tmp));
		}
		tar[key] = temp;
	}
	fin.close();
//	for(auto a : tar)
//	{
//		cout << "[" << a.first << "]" << endl;
//		for(auto b : a.second)
//		{
//			cout << b.first << " " << b.second << endl;
//		}
//	}
//	for (auto currScaffold : scaffolds)
//	{
//		for (auto contig : currScaffold)
//		{
//			cout << contig.first << " " << contig.second << endl;
//		}
//		cout << endl;
//	}
	while(!tar.empty())
	{
		for (auto it = scaffolds.begin(); it != scaffolds.end(); it++)
		{
			for (auto it2 = (*it).begin(); it2 != (*it).end(); it2++)
	 		{
				if (tar.find((*it2).first) != tar.end())
				{
					//cout << "find " << (*it2).first << endl;
//					vector<vector<pair<string, int>>> scaffolds;
//					map<string, vector<pair<string, int>>> tar;
					vector<pair<string, int>> temp;
					for(auto c : tar[(*it2).first])
					{
						//cout << c.first << " " << c.second << endl;
						if(c.second == 1 && (*it2).second == 1)
							temp.pb(make_pair(c.first, 0));
						else if(c.second == 0 && (*it2).second == 1)
							temp.pb(make_pair(c.first, 1));
						else
							temp.pb(make_pair(c.first, c.second));
					}
					if((*it2).second == 1)
					{
						//cout << "reverse" << endl; 
						reverse(temp.begin(), temp.end());
					} 
					tar.erase(tar.find((*it2).first));
					(*it).erase(it2);
					(*it).insert(it2, temp.begin(), temp.end());
					break;
				}
			}
		}
	}
	int scaffoldID = 1;
	ofstream fout(target_dir+"/myScaffold.txt");
	for (auto currScaffold : scaffolds)
	{
		fout << ">scaffold_" << scaffoldID++ << endl;
		for (auto contig : currScaffold)
		{
			fout << contig.first << " " << contig.second << endl;
		}
		fout << endl;
	}

	return 0;
}
