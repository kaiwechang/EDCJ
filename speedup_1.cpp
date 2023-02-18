#include "markerReorder.cpp"
#include <iostream>
#include <fstream>
#include <set>
#include <vector>
#define pb push_back
using namespace std;
vector<int> refFamily, tarFamily; //start with index 0(input marker id = vector index+1)
vector<string> refContig, tarContig;
vector<int> refMarkerNum, tarMarkerNum;
vector<bool> refIso, tarIso; //a marker is iso
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
	ifstream fin(argv[1]);
	int id, family, tmp;
	string contig;
	while (fin >> id >> family >> contig >> tmp)
	{
		refFamily.pb(family);
		refContig.pb(contig);
	}
	fin.close();
	fin.open(argv[2]);
	while (fin >> id >> family >> contig >> tmp)
	{
		tarFamily.pb(family);
		tarContig.pb(contig);
	}
	refMarkerNum.resize(refFamily.size() + 1);
	refIso.resize(refFamily.size(), 1);
	for (int i : refFamily) //count markernum
	{
		i = abs(i);
		refMarkerMax = refMarkerMax > i ? refMarkerMax : i;
		refMarkerNum[i]++;
	}
	for (int i = 0; i < refFamily.size(); i++) // find iso
	{
		int absFam = abs(refFamily[i]);
		if ((refMarkerNum[absFam] <= 1)) //must be iso
			continue;
		for (int j = i + 2; j < refFamily.size(); j++)
		{
			if (refContig[i] != refContig[j])
				break;

			if (refFamily[j] == refFamily[i]) //same sign
			{
				if (j + 1 < refFamily.size() && refFamily[j + 1] == refFamily[i + 1])
				{
					if (refContig[j + 1] == refContig[i + 1] && refContig[i] == refContig[i + 1])
					{
						//cout<<i<<" "<<i+1<<" "<<j<<" "<<j+1<<endl;
						refIso[i] = refIso[j] = refIso[i + 1] = refIso[j + 1] = 0;
					}
				}
			}
			else if (refFamily[j] == -refFamily[i])
			{
				if (refFamily[j - 1] == -refFamily[i + 1])
				{
					if (refContig[j - 1] == refContig[i + 1] && refContig[i] == refContig[i + 1])
					{
						//cout<<i<<" "<<j<<" "<<i+1<<" "<<j-1<<endl;
						refIso[i] = refIso[j] = refIso[i + 1] = refIso[j - 1] = 0;
					}
				}
			}
		}
	}
	tarMarkerNum.resize(tarFamily.size() + 1);
	tarIso.resize(tarFamily.size(), 1);
	for (int i : tarFamily)
	{
		i = abs(i);
		tarMarkerMax = tarMarkerMax > i ? tarMarkerMax : i;
		tarMarkerNum[i]++;
	}
	for (int i = 0; i < tarFamily.size(); i++) // find iso
	{
		int absFam = abs(tarFamily[i]);
		if ((tarMarkerNum[absFam] <= 1)) //must be iso
			continue;
		for (int j = i + 2; j < tarFamily.size(); j++)
		{
			if (tarContig[i] != tarContig[j])
				break;

			if (tarFamily[j] == tarFamily[i]) //same sign
			{
				if (j + 1 < tarFamily.size() && tarFamily[j + 1] == tarFamily[i + 1])
				{
					if (tarContig[j + 1] == tarContig[i + 1] && tarContig[i] == tarContig[i + 1])
					{
						//cout<<i<<" "<<i+1<<" "<<j<<" "<<j+1<<endl;
						tarIso[i] = tarIso[j] = tarIso[i + 1] = tarIso[j + 1] = 0;
					}
				}
			}
			else if (tarFamily[j] == -tarFamily[i])
			{
				if (tarFamily[j - 1] == -tarFamily[i + 1])
				{
					if (tarContig[j - 1] == tarContig[i + 1] && tarContig[i] == tarContig[i + 1])
					{
						//cout<<i<<" "<<j<<" "<<i+1<<" "<<j-1<<endl;
						tarIso[i] = tarIso[j] = tarIso[i + 1] = tarIso[j - 1] = 0;
					}
				}
			}
		}
	}
	showSingletonPercentage();

	bool isend = 0;
	while (!isend)
	{
		cout << "+" << endl;
		isend = 1;
		for (int i = 0; i < refFamily.size() - 1; i++)
		{
			int refAbsFam = abs(refFamily[i]);
			if (refMarkerNum[refAbsFam] == 1) //singleton on ref
			{
				if (refContig[i] == refContig[i + 1]) //iso on ref singlton's right
				{
					for (int j = 0; j < tarFamily.size() - 1; j++)
					{
						if (tarFamily[j] == refFamily[i] && tarMarkerNum[abs(tarFamily[j])] == 1) //shared + singleton on tar
						{
							//cout<<"ref"<<i<<" tar"<<j<<endl;
							if (tarFamily[j + 1] == refFamily[i + 1] && tarContig[j] == tarContig[j + 1]) //shared + iso on tar
							{
								if (refMarkerNum[abs(refFamily[i + 1])] > 1 || tarMarkerNum[abs(tarFamily[j + 1])] > 1) //not all singleton
								{
									isend = 0;
									//cout<<"ref:"<<i<<" "<<i+1<<" tar:"<<j<<" "<<j+1<<endl;
									refMarkerNum[abs(refFamily[i + 1])]--;
									tarMarkerNum[abs(tarFamily[j + 1])]--;

									refFamily[i + 1] = ++refMarkerMax;
									if (refMarkerNum.size() == refMarkerMax) //warning: Fool Proof
										refMarkerNum.pb(1);
									else
										refMarkerNum[refMarkerMax] = 1;

									tarFamily[j + 1] = ++tarMarkerMax;
									if (tarMarkerNum.size() == tarMarkerMax) //warning: Fool Proof
										tarMarkerNum.pb(1);
									else
										tarMarkerNum[tarMarkerMax] = 1;
								}
							}
						}
						//neg
						else if (tarFamily[j] == -refFamily[i] && tarMarkerNum[abs(tarFamily[j])] == 1) //singleton on tar
						{
							if (j - 1 >= 0 && tarFamily[j - 1] == -refFamily[i + 1] && tarContig[j] == tarContig[j - 1]) //shared + iso on tar
							{
								if (refMarkerNum[abs(refFamily[i + 1])] > 1 || tarMarkerNum[abs(tarFamily[j - 1])] > 1)
								{
									isend = 0;
									//cout<<"ref:"<<i<<" "<<i+1<<" tar:"<<j<<" "<<j-1<<endl;
									refMarkerNum[abs(refFamily[i + 1])]--;
									tarMarkerNum[abs(tarFamily[j - 1])]--;

									refFamily[i + 1] = ++refMarkerMax;
									if (refMarkerNum.size() == refMarkerMax) //warning: Fool Proof
										refMarkerNum.pb(1);
									else
										refMarkerNum[refMarkerMax] = 1;
									tarFamily[j - 1] = -(++tarMarkerMax);
									if (tarMarkerNum.size() == tarMarkerMax) //warning: Fool Proof
										tarMarkerNum.pb(1);
									else
										tarMarkerNum[tarMarkerMax] = 1;
								}
							}
						}
					}
				}
				if (i - 1 >= 0 && refContig[i] == refContig[i - 1]) //iso on ref singleton's left
				{
					for (int j = 0; j < tarFamily.size() - 1; j++)
					{
						if (tarFamily[j] == refFamily[i] && tarMarkerNum[abs(tarFamily[j])] == 1) //shared + singleton on tar
						{
							//cout<<"p2ref"<<i<<" tar"<<j<<endl;
							//cout<<tarIso[j-1]<<endl;
							if (j - 1 >= 0 && tarFamily[j - 1] == refFamily[i - 1] && tarContig[j] == tarContig[j - 1]) //shared + iso on tar
							{
								//cout<<"p2+ref"<<i<<" tar"<<j<<endl;
								if (refMarkerNum[abs(refFamily[i - 1])] > 1 || tarMarkerNum[abs(tarFamily[j - 1])] > 1) //not all singleton
								{
									//cout<<"p2++ref"<<i<<" tar"<<j<<endl;
									isend = 0;
									//cout<<"ref:"<<i<<" "<<i+1<<" tar:"<<j<<" "<<j+1<<endl;
									refMarkerNum[abs(refFamily[i - 1])]--;
									tarMarkerNum[abs(tarFamily[j - 1])]--;

									refFamily[i - 1] = ++refMarkerMax;
									if (refMarkerNum.size() == refMarkerMax) //warning: Fool Proof
										refMarkerNum.pb(1);
									else
										refMarkerNum[refMarkerMax] = 1;

									tarFamily[j - 1] = ++tarMarkerMax;
									if (tarMarkerNum.size() == tarMarkerMax) //warning: Fool Proof
										tarMarkerNum.pb(1);
									else
										tarMarkerNum[tarMarkerMax] = 1;
								}
							}
						}
						//neg
						else if (tarFamily[j] == -refFamily[i] && tarMarkerNum[abs(tarFamily[j])] == 1) //singleton on tar
						{
							if (j + 1 < tarFamily.size() && tarFamily[j + 1] == -refFamily[i - 1] && tarContig[j + 1] == tarContig[j]) //shared + iso on tar
							{
								if (refMarkerNum[abs(refFamily[i - 1])] > 1 || tarMarkerNum[abs(tarFamily[j + 1])] > 1)
								{
									isend = 0;
									//cout<<"ref:"<<i<<" "<<i+1<<" tar:"<<j<<" "<<j-1<<endl;
									refMarkerNum[abs(refFamily[i - 1])]--;
									tarMarkerNum[abs(tarFamily[j + 1])]--;

									refFamily[i - 1] = ++refMarkerMax;
									if (refMarkerNum.size() == refMarkerMax)
										refMarkerNum.pb(1);
									else
										refMarkerNum[refMarkerMax] = 1;

									tarFamily[j + 1] = -(++tarMarkerMax);
									if (tarMarkerNum.size() == tarMarkerMax) //warning: Fool Proof
										tarMarkerNum.pb(1);
									else
										tarMarkerNum[tarMarkerMax] = 1;
								}
							}
						}
					}
				}
			}
		}
	}

	string out_dir(argc < 4 ? "output" : argv[3]);
	ofstream fout(out_dir+"/ref_spd1.all");
	int tmpid = 1;
	for (int i = 0; i < refFamily.size(); i++)
	{
		if (tarMarkerNum[abs(refFamily[i])] != 0)
			fout << tmpid++ << " " << refFamily[i] << " " << refContig[i] << " " << 1 << endl;
	}
	fout.close();

	fout.open(out_dir+"/tar_spd1.all");
	tmpid = 1;
	for (int i = 0; i < tarFamily.size(); i++)
	{
		if (refMarkerNum[abs(tarFamily[i])] != 0)
			fout << tmpid++ << " " << tarFamily[i] << " " << tarContig[i] << " " << 1 << endl;
	}
	showSingletonPercentage();

	markerReorder(out_dir+"/ref_spd1.all", out_dir+"/tar_spd1.all");
	return 0;
}
