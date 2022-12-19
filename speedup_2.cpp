#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <string>
#include <algorithm> 

#define pb push_back
#define abs(x) (x > 0 ? x : -x)

using namespace std;

struct Marker{
	int id,family;
	string contig;
	Marker(int id_,int family_,string contig_)
	{
		id=id_;
		family=family_;
		contig=contig_;
	}
	void show()
	{
		cout<<"Marker("<<id<<","<<family<<","<<contig<<")"<<endl;
	}
};

struct SharedSegment{//index start from 0
	vector<int> refMarkerID,tarMarkerID;
};

vector<Marker> refMarkers,tarMarkers;//Markerid = index+1, start from 1
vector<int> refFamilyNum,tarFamilyNum;
int refMarkerSize=0,tarMarkerSize=0;
int refFamilyMax=0,tarFamilyMax=0;
vector<SharedSegment> sharedsegs,PHFSSs;

void getInput(string refFile,string tarFile)
{
	ifstream fin;
	fin.open(refFile);
	int id,family,tmp;
	string contig;
	while(fin>>id>>family>>contig>>tmp)
	{
		Marker marker(id,family,contig);
		refMarkers.pb(marker);
	}
	refMarkerSize=refMarkers.size();
	fin.close();
	fin.open(tarFile);
	while(fin>>id>>family>>contig>>tmp)
	{
		Marker marker(id,family,contig);
		tarMarkers.pb(marker);
	}
	tarMarkerSize=tarMarkers.size();
	fin.close();
}

void calFamilySize()
{
	refFamilyNum.assign(refMarkers.size()+1,0);// at most every one is singleton
	for(auto m:refMarkers)//cal sizes of ref families
	{
		int fam;
		fam = abs(m.family);
		refFamilyMax = refFamilyMax>fam?refFamilyMax:fam;
		refFamilyNum[fam]++;
	}
	
	tarFamilyNum.assign(tarMarkers.size()+1,0);// at most every one is singleton
	for(auto m:tarMarkers)//cal sizes of tar families
	{
		int fam;
		fam = abs(m.family);
		tarFamilyMax = tarFamilyMax>fam?tarFamilyMax:fam;
		tarFamilyNum[fam]++;
	}
}

void findSS()
{
	SharedSegment ss;
	sharedsegs.clear();
	int refRHS,tarRHS;
	string refContig,tarContig;
	int shift = 0;
	for(int refLHS=0;refLHS<refMarkers.size();refLHS++)//ref tar same dir
	{
		refContig=refMarkers[refLHS].contig;
		for(int tarLHS=0;tarLHS<tarMarkers.size();tarLHS++)
		{
			tarContig=tarMarkers[tarLHS].contig;
			refRHS=refLHS;
			tarRHS=tarLHS;
			int cnt=0;
			while(refRHS<refMarkers.size()&&tarRHS<tarMarkers.size()&&
				refMarkers[refRHS].family==tarMarkers[tarRHS].family&&
				refMarkers[refRHS].contig==refContig&&tarMarkers[tarRHS].contig==tarContig)
			{
				refRHS++;
				tarRHS++;
				cnt++;
			}
			if(cnt>0)
			{
				vector<int> refs,tars;
				for(int i=0;i<cnt;i++)
				{
					refs.pb(refLHS+i);
					tars.pb(tarLHS+i);
					
				}
				ss.refMarkerID=refs;
				ss.tarMarkerID=tars;
				sharedsegs.pb(ss);
			}
		}
	}
	for(int refLHS=0;refLHS<refMarkers.size();refLHS++)//ref tar diff dir
	{
		refContig=refMarkers[refLHS].contig;
		for(int tarLHS=tarMarkers.size()-1;tarLHS>=0;tarLHS--)
		{
			tarContig=tarMarkers[tarLHS].contig;
			refRHS=refLHS;
			tarRHS=tarLHS;
			int cnt=0;
			while(refRHS<refMarkers.size()&&tarRHS>=0&&
				refMarkers[refRHS].family==-tarMarkers[tarRHS].family&&
				refMarkers[refRHS].contig==refContig&&tarMarkers[tarRHS].contig==tarContig)
			{
				refRHS++;
				tarRHS--;
				cnt++;
			}
			if(cnt>0)
			{
				vector<int> refs,tars;
				for(int i=0;i<cnt;i++)
				{
					refs.pb(refLHS+i);
					tars.pb(tarLHS-i);
					ss.refMarkerID=refs;
					ss.tarMarkerID=tars;
					
				}
				if(ss.tarMarkerID.size() != 1 && ss.tarMarkerID[0] > ss.tarMarkerID[1])
					reverse(ss.tarMarkerID.begin(), ss.tarMarkerID.end());
				sharedsegs.pb(ss);
			}
		}
	}
}
void findPHFSS()
{
	PHFSSs.clear();
	for(auto ss:PHFSSs)
	{
		bool wanted=1;
		//ref singleton tar not singleton
		for(int ref:ss.refMarkerID)
		{
			if(refFamilyNum[abs(refMarkers[ref].family)]!=1)
			{
				wanted=0;
				break;
			}
		}
		if(wanted)
		{
			for(int tar:ss.tarMarkerID)
			{
				if(tarFamilyNum[abs(tarMarkers[tar].family)]==1)//tar singleton
				{
					wanted=0;
					break;
				}
			}

		}
		if(wanted)
		{
			PHFSSs.pb(ss);
			continue;
		}
		//tar singleton ref not singleton
		wanted=1;
		for(int tar:ss.tarMarkerID)
		{
			if(tarFamilyNum[abs(tarMarkers[tar].family)]!=1)
			{
				wanted=0;
				break;
			}
		}
		if(wanted)
		{
			for(int ref:ss.refMarkerID)
			{
				if(refFamilyNum[abs(refMarkers[ref].family)]==1)//ref singleton
				{
					wanted=0;
					break;
				}
			}
		}
		if(wanted)
			PHFSSs.pb(ss);
	}
}
int main(int argc, char *argv[])
{
	getInput(argv[1],argv[2]);
	calFamilySize();
	bool isend = 0;
	while(!isend)
	{
		isend = 1;
		findSS();
		findPHFSS();
		int refLHS, refRHS, tarLHS, tarRHS;
		for(int i=0;i<PHFSSs.size();i++)
		{
			cout << "ref: " << endl;
			for(int id : PHFSSs[i].refMarkerID)
				cout << id << " " << refMarkers[id].family << endl;
			cout << "tar: " << endl;
			for(int id : PHFSSs[i].tarMarkerID)
				cout << id << " " << tarMarkers[id].family << endl;
			cout << endl;
			refLHS = PHFSSs[i].refMarkerID.front();
			refRHS = PHFSSs[i].refMarkerID.back();
			tarLHS = PHFSSs[i].tarMarkerID.front();
			tarRHS = PHFSSs[i].tarMarkerID.back();
			if(refLHS == 0 || refRHS == refMarkers.size()-1 || tarLHS == 0 || tarRHS == tarMarkers.size()-1)
			{
				continue;
			}
			if(refMarkers[refLHS].contig != refMarkers[refLHS-1].contig || refMarkers[refRHS].contig != refMarkers[refRHS+1].contig ||
				tarMarkers[tarLHS].contig != tarMarkers[tarLHS-1].contig || tarMarkers[tarRHS].contig != tarMarkers[tarRHS+1].contig)
			{
				continue;
			}
			if(refMarkers[refLHS].family == tarMarkers[tarLHS].family)
			{
				if(refMarkers[refLHS-1].family == -tarMarkers[tarRHS+1].family && refMarkers[refRHS+1].family == -tarMarkers[tarLHS-1].family)
				{
					if(refFamilyNum[abs(refMarkers[refLHS-1].family)] == tarFamilyNum[abs(tarMarkers[tarRHS+1].family)] ==
						refFamilyNum[abs(refMarkers[refRHS+1].family)] == tarFamilyNum[abs(tarMarkers[tarLHS-1].family)] == 1)
					{
						isend = 0;
						for(int id : PHFSSs[i].refMarkerID)
						{
							refFamilyNum[abs(refMarkers[id].family)]--;
							refMarkers[id].family = ++refFamilyMax;
							if (refFamilyNum.size() == refFamilyMax)
								refFamilyNum.pb(1);
							else
								refFamilyNum[refFamilyMax] = 1;
						}
						for(int id : PHFSSs[i].tarMarkerID)
						{
							tarFamilyNum[abs(tarMarkers[id].family)]--;
							tarMarkers[id].family = ++tarFamilyMax;
							if (tarFamilyNum.size() == tarFamilyMax)
								tarFamilyNum.pb(1);
							else
								tarFamilyNum[tarFamilyMax] = 1;
						}
						cout << "find" << endl;
						// int n;
						// cin >> n;
					}
				}
			}
			else
			{
				if(refMarkers[refLHS-1].family == tarMarkers[tarLHS-1].family && refMarkers[refRHS+1].family == tarMarkers[tarRHS+1].family)
				{
					if(refFamilyNum[abs(refMarkers[refLHS-1].family)] == tarFamilyNum[abs(tarMarkers[tarRHS+1].family)] ==
						refFamilyNum[abs(refMarkers[refRHS+1].family)] == tarFamilyNum[abs(tarMarkers[tarLHS-1].family)] == 1)
					{
						isend = 0;
						for(int id : PHFSSs[i].refMarkerID)
						{
							refFamilyNum[abs(refMarkers[id].family)]--;
							refMarkers[id].family = ++refFamilyMax;
							if (refFamilyNum.size() == refFamilyMax)
								refFamilyNum.pb(1);
							else
								refFamilyNum[refFamilyMax] = 1;
						}
						vector<int> temp = PHFSSs[i].tarMarkerID;
						reverse(temp.begin(), temp.end());
						for(int id : temp)
						{
							tarFamilyNum[abs(tarMarkers[id].family)]--;
							tarMarkers[id].family = -(++tarFamilyMax);
							if (tarFamilyNum.size() == tarFamilyMax)
								tarFamilyNum.pb(1);
							else
								tarFamilyNum[tarFamilyMax] = 1;
						}
						cout << "find" << endl;
						// int n;
						// cin >> n;
					}
				}
			}
		}
	}

	string out_dir(argc < 4 ? "output" : argv[3]);
	ofstream fout(out_dir+"/ref_spd2.all");
	int tmpid=1;
	for(int i=0;i<refMarkers.size();i++)
	{
		if(tarFamilyNum[abs(refMarkers[i].family)]!=0)
			fout<<tmpid++<<" "<<refMarkers[i].family<<" "<<refMarkers[i].contig<<" "<<1<<endl;
	}
	fout.close();
	
	fout.open(out_dir+"/tar_spd2.all");
	tmpid=1;
	for(int i=0;i<tarMarkers.size();i++)
	{
		if(refFamilyNum[abs(tarMarkers[i].family)]!=0)
			fout<<tmpid++<<" "<<tarMarkers[i].family<<" "<<tarMarkers[i].contig<<" "<<1<<endl;
	}
	fout.close();
	return 0;
}
