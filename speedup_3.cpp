#include <iostream>
#include <fstream>
#include <set>
#include <vector>
#include <map>
#include <algorithm> 
#define pb push_back
#define abs(x) (x > 0 ? x : -x)
using namespace std;
vector<int> refFamily, tarFamily; //start with index 0(input marker id = vector index+1)
vector<string> refContig, tarContig;
vector<int> refMarkerNum, tarMarkerNum;
map<string, vector<pair<string, int>>> tarContigMerge, refContigMerge;
int refMarkerMax = 0, tarMarkerMax = 0, refContigSize = 0, tarContigSize = 0;
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
void showNumberOfMarkers()
{
	cout << "ref: " << refFamily.size() << endl;
	cout << "tar: " << tarFamily.size() << endl;	
}
int main(int argc, char *argv[])
{
	ifstream fin(argv[1]);
	int id, family, tmp, preFamily = 0;
	string contig, preContig = "";
	while (fin >> id >> family >> contig >> tmp)
	{
		if(preContig != contig)
			refContigSize++;
		preFamily = family;
		preContig = contig;
		refFamily.pb(family);
		refContig.pb(contig);
	}
	fin.close();
	preFamily = 0, preContig = "";
	fin.open(argv[2]);
	while (fin >> id >> family >> contig >> tmp)
	{
		if(preContig != contig)
			tarContigSize++;
		preFamily = family;
		preContig = contig;
		tarFamily.pb(family);
		tarContig.pb(contig);
	}
	refMarkerNum.resize(refFamily.size() + 1);

	for (int i : refFamily) //count markernum
	{
		i = abs(i);
		refMarkerMax = refMarkerMax > i ? refMarkerMax : i;
		refMarkerNum[i]++;
	}

	tarMarkerNum.resize(tarFamily.size() + 1);
	for (int i : tarFamily) //count markernum
	{
		i = abs(i);
		tarMarkerMax = tarMarkerMax > i ? tarMarkerMax : i;
		tarMarkerNum[i]++;
	}

//	for(int i=0; i < refFamily.size(); i++)
//	{
//		if(refTelomere[i]==1)
//		{
//			cout << "telo:" << refFamily[i] << " ";
//		}
//	}
//	cout << endl;
//	for(int i=0; i < tarFamily.size(); i++)
//	{
//		if(tarTelomere[i]==1)
//		{
//			cout << "telo:" << tarFamily[i] << " ";
//		}
//	}
//	cout << endl;
	showNumberOfMarkers();
	auto rfit = refFamily.begin();
	auto rcit = refContig.begin();
	auto tfit = tarFamily.begin(); 
	auto tcit = tarContig.begin();
	int tarNumberMerge = 0, refNumberMerge = 0;
	// bool isend = 0;
	//speedup6
	//merge ref contig
	mergeRefContig:
	for(rfit = refFamily.begin(), rcit = refContig.begin(); rfit+1 != refFamily.end(); rfit++, rcit++)//find first telo on ref
	{
		if(refMarkerNum[abs(*rfit)] == 1)
		{
			if(rfit == refFamily.begin() || *rcit != *(rcit+1) || *rcit != *(rcit-1))
			{
				auto rfit2 = rfit+1;
				auto rcit2 = rcit+1;
				for(;rfit2 != refFamily.end(); rfit2++, rcit2++)//find second telo on ref
				{
					if(refMarkerNum[abs(*rfit2)] == 1)
					{
						if(*rcit != *rcit2)
						{
							if(rfit2+1 == refFamily.end() || *rcit2 != *(rcit2+1) || *rcit2 != *(rcit2-1))//2 telo found
							{
								//cout << *tfit << " " << *tfit2 << endl;
								//L+L									
								if((rfit == refFamily.begin() || *rcit != *(rcit-1)) && *rcit2 != *(rcit2-1))
								{						
									tfit = tarFamily.begin();
									tcit = tarContig.begin();
									for(; tfit+1 != tarFamily.end(); tfit++, tcit++)
									{
										if(tarMarkerNum[abs(*tfit)] == 1 && tarMarkerNum[abs(*(tfit+1))] == 1 && *tcit == *(tcit+1))
										{
											if((*rfit == -*tfit && *rfit2 == *(tfit+1)) || (*rfit == *(tfit+1) && *rfit2 == -*tfit) )
											{
												cout << "*" << endl;
												refNumberMerge++; 
//													for(int i : tarFamily)
//														cout << i << "\t ";
//													cout << endl;
//													for(string s : tarContig)
//														cout << s << " ";
//													cout << endl
												// cout << "*r1*" << endl;
												// cout << *rcit2 << " inverse & insert to left of " << *rcit << endl;
												if(refContigMerge.find(*rcit) != refContigMerge.end())
												{
													vector<pair<string, int>> temp = refContigMerge[*rcit];
													temp.insert(temp.begin(), make_pair(*rcit2, 1));
													refContigMerge[*rcit] = temp;
												}
												else
												{
													vector<pair<string, int>> temp({make_pair(*rcit2, 1), make_pair(*rcit, 0)});
													refContigMerge[*rcit] = temp;
												}
												auto pos_1_f = rfit2;
												auto pos_1_c = rcit2;
												auto pos_2_f = rfit2;
												auto pos_2_c = rcit2;
												for(;pos_2_c != refContig.end() && *pos_2_c == *rcit2; pos_2_f++, pos_2_c++)
												{
													*pos_2_f = -*pos_2_f;
												}
												vector<int> inf(pos_1_f, pos_2_f);
												reverse(inf.begin(), inf.end());
												refFamily.erase(pos_1_f, pos_2_f);
												refFamily.insert(rfit, inf.begin(), inf.end());
												vector<string> inc(inf.size(), *rcit);
												refContig.erase(pos_1_c, pos_2_c);
												refContig.insert(rcit, inc.begin(), inc.end());
//													for(int i : tarFamily)
//														cout << i << "\t ";
//													cout << endl;
//													for(string s : tarContig)
//														cout << s << " ";
//													cout << endl;
												goto mergeRefContig;
											}
										}
									}
								}
								//L+R		
								if((rfit == refFamily.begin() || *rcit != *(rcit-1)) && (rfit2+1 == refFamily.end() || *rcit2 != *(rcit2+1)))
								{				
									tfit = tarFamily.begin();
									tcit = tarContig.begin();
									for(; tfit+1 != tarFamily.end(); tfit++, tcit++)
									{
										if(tarMarkerNum[abs(*tfit)] == 1 && tarMarkerNum[abs(*(tfit+1))] == 1 && *tcit == *(tcit+1))
										{
											if((*rfit == *(tfit+1) && *rfit2 == *tfit) || (*rfit == -*tfit && *rfit2 == -*(tfit+1)))
											{
												cout << "*" << endl;
												refNumberMerge++; 
//													for(int i : tarFamily)
//														cout << i << "\t ";
//													cout << endl;
//													for(string s : tarContig)
//														cout << s << " ";
//													cout << endl;
												// cout << "*r2*" << endl;
												// cout << *rcit2 << " insert to left of " << *rcit << endl;
												if(refContigMerge.find(*rcit) != refContigMerge.end())
												{
													vector<pair<string, int>> temp = refContigMerge[*rcit];
													temp.insert(temp.begin(), make_pair(*rcit2, 0));
													refContigMerge[*rcit] = temp;
												}
												else
												{
													vector<pair<string, int>> temp({make_pair(*rcit2, 0), make_pair(*rcit, 0)});
													refContigMerge[*rcit] = temp;
												}
												auto pos_1_f = rfit2;
												auto pos_1_c = rcit2;
												auto pos_2_f = rfit2+1;
												auto pos_2_c = rcit2+1;
												for(;pos_1_c != refContig.begin() && *(pos_1_c-1) == *rcit2; pos_1_f--, pos_1_c--);
												vector<int> inf(pos_1_f, pos_2_f);
												refFamily.erase(pos_1_f, pos_2_f);
												refFamily.insert(rfit, inf.begin(), inf.end());
												vector<string> inc(inf.size(), *rcit);
												refContig.erase(pos_1_c, pos_2_c);
												refContig.insert(rcit, inc.begin(), inc.end());
//													for(int i : tarFamily)
//														cout << i << "\t ";
//													cout << endl;
//													for(string s : tarContig)
//														cout << s << " ";
//													cout << endl;
												goto mergeRefContig;
											}
										}
									}
								}
								//R+L
								if((rfit+1 == refFamily.end() || *rcit != *(rcit+1)) && (rfit2 == refFamily.begin() || *rcit2 != *(rcit2-1)))
								{				
									tfit = tarFamily.begin();
									tcit = tarContig.begin();
									for(; tfit+1 != tarFamily.end(); tfit++, tcit++)
									{
										if(tarMarkerNum[abs(*tfit)] == 1 && tarMarkerNum[abs(*(tfit+1))] == 1 && *tcit == *(tcit+1))
										{
											if((*rfit == *tfit && *rfit2 == *(tfit+1)) || (*rfit == -*(tfit+1) && *rfit2 == -*tfit))
											{
												cout << "*" << endl;
												refNumberMerge++; 
//													for(int i : tarFamily)
//														cout << i << "\t ";
//													cout << endl;
//													for(string s : tarContig)
//														cout << s << " ";
//													cout << endl;
												// cout << "*r3*" << endl;
												// cout << *rcit2 << " insert to right of " << *rcit << endl;
												if(refContigMerge.find(*rcit) != refContigMerge.end())
												{
													vector<pair<string, int>> temp = refContigMerge[*rcit];
													temp.insert(temp.end(), make_pair(*rcit2, 0));
													refContigMerge[*rcit] = temp;
												}
												else
												{
													vector<pair<string, int>> temp({make_pair(*rcit, 0), make_pair(*rcit2, 0)});
													refContigMerge[*rcit] = temp;
												}
												auto pos_1_f = rfit2;
												auto pos_1_c = rcit2;
												auto pos_2_f = rfit2;
												auto pos_2_c = rcit2;
												for(;pos_2_c != refContig.end() && *pos_2_c == *rcit2; pos_2_f++, pos_2_c++);
												vector<int> inf(pos_1_f, pos_2_f);
												refFamily.erase(pos_1_f, pos_2_f);
												refFamily.insert(rfit+1, inf.begin(), inf.end());
												vector<string> inc(inf.size(), *rcit);
												refContig.erase(pos_1_c, pos_2_c);
												refContig.insert(rcit+1, inc.begin(), inc.end());
//													for(int i : tarFamily)
//														cout << i << "\t ";
//													cout << endl;
//													for(string s : tarContig)
//														cout << s << " ";
//													cout << endl;
												goto mergeRefContig;
											}
										}
									}
								}
								//R+R
								if((rfit+1 == refFamily.end() || *rcit != *(rcit+1)) && (rfit2+1 == refFamily.end() || *rcit2 != *(rcit2+1)))
								{
									tfit = tarFamily.begin();
									tcit = tarContig.begin();
									for(; tfit+1 != tarFamily.end(); tfit++, tcit++)
									{
										if(tarMarkerNum[abs(*tfit)] == 1 && tarMarkerNum[abs(*(tfit+1))] == 1 && *tcit == *(tcit+1))
										{
											if((*rfit == -*(tfit+1) && *rfit2 == *tfit) || (*rfit == *tfit && *rfit2 == -*(tfit+1)))
											{
												cout << "*" << endl;
												refNumberMerge++; 
//													for(int i : tarFamily)
//														cout << i << "\t ";
//													cout << endl;
//													for(string s : tarContig)
//														cout << s << " ";
//													cout << endl;
												// cout << "*r4*" << endl;
												// cout << *rcit2 << " inverse & insert to right of " << *rcit << endl;
												if(refContigMerge.find(*rcit) != refContigMerge.end())
												{
													vector<pair<string, int>> temp = refContigMerge[*rcit];
													temp.insert(temp.end(), make_pair(*rcit2, 1));
													refContigMerge[*rcit] = temp;
												}
												else
												{
													vector<pair<string, int>> temp({make_pair(*rcit, 0), make_pair(*rcit2, 1)});
													refContigMerge[*rcit] = temp;
												}
												auto pos_1_f = rfit2;
												auto pos_1_c = rcit2;
												auto pos_2_f = rfit2+1;
												auto pos_2_c = rcit2+1;
												for(;pos_1_c != refContig.begin() && *(pos_1_c) == *rcit2; pos_1_f--, pos_1_c--)
												{
													*pos_1_f = -*pos_1_f;
												}
												pos_1_f++, pos_1_c++;
												vector<int> inf(pos_1_f, pos_2_f);
												reverse(inf.begin(), inf.end());
												refFamily.erase(pos_1_f, pos_2_f);
												refFamily.insert(rfit+1, inf.begin(), inf.end());
												vector<string> inc(inf.size(), *rcit);
												refContig.erase(pos_1_c, pos_2_c);
												refContig.insert(rcit+1, inc.begin(), inc.end());
//													for(int i : tarFamily)
//														cout << i << "\t ";
//													cout << endl;
//													for(string s : tarContig)
//														cout << s << " ";
//													cout << endl;
												goto mergeRefContig;
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
	//merge tar contig
	mergeTarContig:
	for(tfit = tarFamily.begin(), tcit = tarContig.begin();tfit+1 != tarFamily.end(); tfit++, tcit++)//find first telo on tar
	{
		if(tarMarkerNum[abs(*tfit)] == 1)
		{
			if(tfit == tarFamily.begin() || *tcit != *(tcit+1) || *tcit != *(tcit-1))
			{
				auto tfit2 = tfit+1;
				auto tcit2 = tcit+1;
				for(;tfit2 != tarFamily.end(); tfit2++, tcit2++)//find second telo on tar
				{
					if(tarMarkerNum[abs(*tfit2)] == 1)
					{
						if(*tcit != *tcit2)
						{
							if(tfit2+1 == tarFamily.end() || *tcit2 != *(tcit2+1) || *tcit2 != *(tcit2-1))//2 telo found
							{
								//cout << *tfit << " " << *tfit2 << endl;
								//L+L									
								if((tfit == tarFamily.begin() || *tcit != *(tcit-1)) && *tcit2 != *(tcit2-1))
								{						
									rfit = refFamily.begin();
									rcit = refContig.begin();
									for(; rfit+1 != refFamily.end(); rfit++, rcit++)
									{
										if(refMarkerNum[abs(*rfit)] == 1 && refMarkerNum[abs(*(rfit+1))] == 1 && *rcit == *(rcit+1))
										{
											if((*tfit == -*rfit && *tfit2 == *(rfit+1)) || (*tfit == *(rfit+1) && *tfit2 == -*rfit) )
											{
												cout << "+" << endl;
												tarNumberMerge++; 
//													for(int i : tarFamily)
//														cout << i << "\t ";
//													cout << endl;
//													for(string s : tarContig)
//														cout << s << " ";
//													cout << endl
												// cout << "*t1*" << endl;
												// cout << *tcit2 << " inverse & insert to left of " << *tcit << endl;
												if(tarContigMerge.find(*tcit) != tarContigMerge.end())
												{
													vector<pair<string, int>> temp = tarContigMerge[*tcit];
													temp.insert(temp.begin(), make_pair(*tcit2, 1));
													tarContigMerge[*tcit] = temp;
												}
												else
												{
													vector<pair<string, int>> temp({make_pair(*tcit2, 1), make_pair(*tcit, 0)});
													tarContigMerge[*tcit] = temp;
												}
												auto pos_1_f = tfit2;
												auto pos_1_c = tcit2;
												auto pos_2_f = tfit2;
												auto pos_2_c = tcit2;
												for(;pos_2_c != tarContig.end() && *pos_2_c == *tcit2; pos_2_f++, pos_2_c++)
												{
													*pos_2_f = -*pos_2_f;
												}
												vector<int> inf(pos_1_f, pos_2_f);
												reverse(inf.begin(), inf.end());
												tarFamily.erase(pos_1_f, pos_2_f);
												tarFamily.insert(tfit, inf.begin(), inf.end());
												vector<string> inc(inf.size(), *tcit);
												tarContig.erase(pos_1_c, pos_2_c);
												tarContig.insert(tcit, inc.begin(), inc.end());
//													for(int i : tarFamily)
//														cout << i << "\t ";
//													cout << endl;
//													for(string s : tarContig)
//														cout << s << " ";
//													cout << endl;
												goto mergeTarContig;
											}
										}
									}
								}
								//L+R		
								if((tfit == tarFamily.begin() || *tcit != *(tcit-1)) && (tfit2+1 == tarFamily.end() || *tcit2 != *(tcit2+1)))
								{				
									rfit = refFamily.begin();
									rcit = refContig.begin();
									for(; rfit+1 != refFamily.end(); rfit++, rcit++)
									{
										if(refMarkerNum[abs(*rfit)] == 1 && refMarkerNum[abs(*(rfit+1))] == 1 && *rcit == *(rcit+1))
										{
											if((*tfit == *(rfit+1) && *tfit2 == *rfit) || (*tfit == -*rfit && *tfit2 == -*(rfit+1)))
											{
												cout << "+" << endl;
												tarNumberMerge++; 
//													for(int i : tarFamily)
//														cout << i << "\t ";
//													cout << endl;
//													for(string s : tarContig)
//														cout << s << " ";
//													cout << endl;
												// cout << "*t2*" << endl;
												// cout << *tcit2 << " insert to left of " << *tcit << endl;
												if(tarContigMerge.find(*tcit) != tarContigMerge.end())
												{
													vector<pair<string, int>> temp = tarContigMerge[*tcit];
													temp.insert(temp.begin(), make_pair(*tcit2, 0));
													tarContigMerge[*tcit] = temp;
												}
												else
												{
													vector<pair<string, int>> temp({make_pair(*tcit2, 0), make_pair(*tcit, 0)});
													tarContigMerge[*tcit] = temp;
												}
												auto pos_1_f = tfit2;
												auto pos_1_c = tcit2;
												auto pos_2_f = tfit2+1;
												auto pos_2_c = tcit2+1;
												for(;pos_1_c != tarContig.begin() && *(pos_1_c-1) == *tcit2; pos_1_f--, pos_1_c--);
												vector<int> inf(pos_1_f, pos_2_f);
												tarFamily.erase(pos_1_f, pos_2_f);
												tarFamily.insert(tfit, inf.begin(), inf.end());
												vector<string> inc(inf.size(), *tcit);
												tarContig.erase(pos_1_c, pos_2_c);
												tarContig.insert(tcit, inc.begin(), inc.end());
//													for(int i : tarFamily)
//														cout << i << "\t ";
//													cout << endl;
//													for(string s : tarContig)
//														cout << s << " ";
//													cout << endl;
												goto mergeTarContig;
											}
										}
									}
								}
								//R+L
								if((tfit+1 == tarFamily.end() || *tcit != *(tcit+1)) && (tfit2 == tarFamily.begin() || *tcit2 != *(tcit2-1)))
								{				
									rfit = refFamily.begin();
									rcit = refContig.begin();
									for(; rfit+1 != refFamily.end(); rfit++, rcit++)
									{
										if(refMarkerNum[abs(*rfit)] == 1 && refMarkerNum[abs(*(rfit+1))] == 1 && *rcit == *(rcit+1))
										{
											if((*tfit == *rfit && *tfit2 == *(rfit+1)) || (*tfit == -*(rfit+1) && *tfit2 == -*rfit))
											{
												cout << "+" << endl;
												tarNumberMerge++; 
//													for(int i : tarFamily)
//														cout << i << "\t ";
//													cout << endl;
//													for(string s : tarContig)
//														cout << s << " ";
//													cout << endl;
												// cout << "*t3*" << endl;
												// cout << *tcit2 << " insert to right of " << *tcit << endl;
												if(tarContigMerge.find(*tcit) != tarContigMerge.end())
												{
													vector<pair<string, int>> temp = tarContigMerge[*tcit];
													temp.insert(temp.end(), make_pair(*tcit2, 0));
													tarContigMerge[*tcit] = temp;
												}
												else
												{
													vector<pair<string, int>> temp({make_pair(*tcit, 0), make_pair(*tcit2, 0)});
													tarContigMerge[*tcit] = temp;
												}
												auto pos_1_f = tfit2;
												auto pos_1_c = tcit2;
												auto pos_2_f = tfit2;
												auto pos_2_c = tcit2;
												for(;pos_2_c != tarContig.end() && *pos_2_c == *tcit2; pos_2_f++, pos_2_c++);
												vector<int> inf(pos_1_f, pos_2_f);
												tarFamily.erase(pos_1_f, pos_2_f);
												tarFamily.insert(tfit+1, inf.begin(), inf.end());
												vector<string> inc(inf.size(), *tcit);
												tarContig.erase(pos_1_c, pos_2_c);
												tarContig.insert(tcit+1, inc.begin(), inc.end());
//													for(int i : tarFamily)
//														cout << i << "\t ";
//													cout << endl;
//													for(string s : tarContig)
//														cout << s << " ";
//													cout << endl;
												goto mergeTarContig;
											}
										}
									}
								}
								//R+R
								if((tfit+1 == tarFamily.end() || *tcit != *(tcit+1)) && (tfit2+1 == tarFamily.end() || *tcit2 != *(tcit2+1)))
								{
									rfit = refFamily.begin();
									rcit = refContig.begin();
									for(; rfit+1 != refFamily.end(); rfit++, rcit++)
									{
										if(refMarkerNum[abs(*rfit)] == 1 && refMarkerNum[abs(*(rfit+1))] == 1 && *rcit == *(rcit+1))
										{
											if((*tfit == -*(rfit+1) && *tfit2 == *rfit) || (*tfit == *rfit && *tfit2 == -*(rfit+1)))
											{
												cout << "+" << endl;
												tarNumberMerge++; 
//													for(int i : tarFamily)
//														cout << i << "\t ";
//													cout << endl;
//													for(string s : tarContig)
//														cout << s << " ";
//													cout << endl;
												// cout << "*t4*" << endl;
												// cout << *tcit2 << " inverse & insert to right of " << *tcit << endl;
												if(tarContigMerge.find(*tcit) != tarContigMerge.end())
												{
													vector<pair<string, int>> temp = tarContigMerge[*tcit];
													temp.insert(temp.end(), make_pair(*tcit2, 1));
													tarContigMerge[*tcit] = temp;
												}
												else
												{
													vector<pair<string, int>> temp({make_pair(*tcit, 0), make_pair(*tcit2, 1)});
													tarContigMerge[*tcit] = temp;
												}
												auto pos_1_f = tfit2;
												auto pos_1_c = tcit2;
												auto pos_2_f = tfit2+1;
												auto pos_2_c = tcit2+1;
												for(;pos_1_c != tarContig.begin() && *(pos_1_c) == *tcit2; pos_1_f--, pos_1_c--)
												{
													*pos_1_f = -*pos_1_f;
												}
												pos_1_f++, pos_1_c++;
												vector<int> inf(pos_1_f, pos_2_f);
												reverse(inf.begin(), inf.end());
												tarFamily.erase(pos_1_f, pos_2_f);
												tarFamily.insert(tfit+1, inf.begin(), inf.end());
												vector<string> inc(inf.size(), *tcit);
												tarContig.erase(pos_1_c, pos_2_c);
												tarContig.insert(tcit+1, inc.begin(), inc.end());
//													for(int i : tarFamily)
//														cout << i << "\t ";
//													cout << endl;
//													for(string s : tarContig)
//														cout << s << " ";
//													cout << endl;
												goto mergeTarContig;
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
	// while (!isend)
	// {
	// 	isend = 1;
		// while(rfit != refFamily.end())//find consecutive same marker on ref
		// {
		// 	if(abs(*rfit) == abs(preFamily) && *rcit == preContig)
		// 	{
		// 		//cout << "find consecutive same marker on ref:" << *rfit << " Contig:" << *rcit << endl;
		// 		isend = 0;
		// 		refMarkerNum[abs(*rfit)]--;
		// 		refFamily.erase(rfit);
		// 		refContig.erase(rcit);
		// 	}
		// 	else
		// 	{
		// 		preFamily = *rfit, preContig = *rcit; 
		// 		rfit++, rcit++;
		// 	}
		// }
		// preFamily = 0, preContig = "";
		// tfit = tarFamily.begin(), tcit = tarContig.begin();
		// while(tfit != tarFamily.end())//find consecutive same marker on tar
		// {
		// 	if(abs(*tfit) == abs(preFamily) && *tcit == preContig)
		// 	{
		// 		//cout << "find consecutive same marker on tar:" << *tfit << " Contig:" << *tcit << endl;
		// 		isend = 0;
		// 		tarMarkerNum[abs(*rfit)]--;
		// 		tarFamily.erase(tfit);
		// 		tarContig.erase(tcit);
		// 	}
		// 	else
		// 	{
		// 		preFamily = *tfit, preContig = *tcit; 
		// 		tfit++, tcit++;
		// 	}
		// }
		
		
		
//		//fix cycles of length two
//		fix:
//		//case1
//		for (rfit = refFamily.begin(), rcit = refContig.begin(); rfit+1 != refFamily.end(); rfit++, rcit++)//for each marker on ref
//		{
//			if (refMarkerNum[abs(*rfit)] == 1) //singleton
//			{	
//				//cout << "find singleton on ref: " << *rfit << endl;
//				if (*rcit == *(rcit+1)) //right is same contig 
//				{					
//					for (tfit = tarFamily.begin(), tcit = tarContig.begin(); tfit+1 != tarFamily.end(); tfit++, tcit++)//for each marker on tar 
//					{
//						if(*tfit == *rfit && tarMarkerNum[abs(*tfit)] == 1)//same sign & same family & singleton 
//						{
//							//cout << "find singleton on tar: " << *tfit << endl;
//							if(*(tfit+1) == *(rfit+1) && *tcit == *(tcit+1))//right is same sign & family & same contig
//							{	
//								if(refMarkerNum[abs(*(rfit+1))] > 1 || tarMarkerNum[abs(*(tfit+1))] > 1)
//								{ 
//								// cout << "1 remove ref:" << *(rfit) << " " << *(rfit+1) << "*rcit:" << *rcit << " " << "*(rcit+1):" << *(rcit+1) << endl;
//								// cout << "1 remove tar:" << *(tfit) << " " << *(tfit+1) << "*tcit:" << *tcit << " " << "*(tcit+1):" << *(tcit+1) << endl;
//								isend = 0;
//								refMarkerNum[abs(*(rfit+1))]--;
//								tarMarkerNum[abs(*(tfit+1))]--;
//								*(rfit+1) = ++refMarkerMax;
//								if (refMarkerNum.size() == refMarkerMax) //warning: Fool Proof
//									refMarkerNum.pb(1);
//								else
//									refMarkerNum[refMarkerMax] = 1;
//
//								*(tfit+1) = ++tarMarkerMax;
//								if (tarMarkerNum.size() == tarMarkerMax) //warning: Fool Proof
//									tarMarkerNum.pb(1);
//								else
//									tarMarkerNum[tarMarkerMax] = 1;
////								refFamily.erase(rfit+1);
////								refContig.erase(rcit+1);
////								tarFamily.erase(tfit+1);
////								tarContig.erase(tcit+1);
//								count++;
//								//goto fix;
//								}
//							}
//						}
//					}
//				}
//			}
//		}
//		//case2
//		for (rfit = refFamily.begin(), rcit = refContig.begin(); rfit+1 != refFamily.end(); rfit++, rcit++)//for each marker on ref
//		{
//			if (refMarkerNum[abs(*rfit)] == 1) //singleton
//			{	
//				if (*rcit == *(rcit+1)) //same contig 
//				{
//					for (tfit = tarFamily.end()-1, tcit = tarContig.end()-1; tfit != tarFamily.begin(); tfit--, tcit--)//for each marker on tar 
//					{
//						
//						if (*tfit == -*rfit && tarMarkerNum[abs(*tfit)] == 1)//different sign & same family & singleton 
//						{
//							
//							if (*(tfit-1) == -*(rfit+1) && *tcit == *(tcit-1))//left is differnrt sign & same family & same contig
//							{
//								if(refMarkerNum[abs(*(rfit+1))] > 1 || tarMarkerNum[abs(*(tfit-1))] > 1)
//								{
//								// cout << "2 remove ref:" << *(rfit) << " " << *(rfit+1) << "*rcit:" << *rcit << " " << "*(rcit+1):" << *(rcit+1) << endl;
//								// cout << "2 remove tar:" << *(tfit) << " " << *(tfit-1) << "*tcit:" << *tcit << " " << "*(tcit-1):" << *(tcit-1) << endl;
//								isend = 0;								
//								refMarkerNum[abs(*(rfit+1))]--;
//								tarMarkerNum[abs(*(tfit-1))]--;
//								*(rfit+1) = ++refMarkerMax;
//								if (refMarkerNum.size() == refMarkerMax) //warning: Fool Proof
//									refMarkerNum.pb(1);
//								else
//								refMarkerNum[refMarkerMax] = 1;
//								*(tfit-1) = -(++tarMarkerMax);
//								if (tarMarkerNum.size() == tarMarkerMax) //warning: Fool Proof
//									tarMarkerNum.pb(1);
//								else
//									tarMarkerNum[tarMarkerMax] = 1;
////								refFamily.erase(rfit+1);
////								refContig.erase(rcit+1);
////								tarFamily.erase(tfit-1);
////								tarContig.erase(tcit-1);
//								count++;
//								//goto fix;
//								}
//							}
//						}
//					}
//				}
//			}
//		}
//		//case3
//		for (rfit = refFamily.end()-1, rcit = refContig.end()-1; rfit != refFamily.begin(); rfit--, rcit--)//for each marker on ref
//		{
//			if (refMarkerNum[abs(*rfit)] == 1)//singleton
//			{
//				if (*rcit == *(rcit-1))//left is same contig 
//				{
//					for (tfit = tarFamily.end()-1, tcit = tarContig.end()-1; tfit != tarFamily.begin(); tfit--, tcit--)//for each marker on tar
//					{
//						if(*tfit == *rfit && tarMarkerNum[abs(*tfit)] == 1)//same sign & same family & singleton 
//						{	
//							if(*(tfit-1) == *(rfit-1) && *tcit == *(tcit-1))//left is same sign & same family & same contig
//							{	
//								if(refMarkerNum[abs(*(rfit-1))] > 1 || tarMarkerNum[abs(*(tfit-1))] > 1)
//								{
//								// cout << "3 remove ref:" << *(rfit) << " " << *(rfit-1) << "*rcit:" << *rcit << " " << "*(rcit-1):" << *(rcit-1) << endl;
//								// cout << "3 remove tar:" << *(tfit) << " " << *(tfit-1) << "*tcit:" << *tcit << " " << "*(tcit-1):" << *(tcit-1) << endl;
//								isend = 0;
//								refMarkerNum[abs(*(rfit-1))]--;
//								tarMarkerNum[abs(*(tfit-1))]--;
//								*(rfit-1) = ++refMarkerMax;
//								if (refMarkerNum.size() == refMarkerMax) //warning: Fool Proof
//									refMarkerNum.pb(1);
//								else
//									refMarkerNum[refMarkerMax] = 1;
//
//								*(tfit-1) = ++tarMarkerMax;
//								if (tarMarkerNum.size() == tarMarkerMax) //warning: Fool Proof
//									tarMarkerNum.pb(1);
//								else
//									tarMarkerNum[tarMarkerMax] = 1;
////								refFamily.erase(rfit-1);
////								refContig.erase(rcit-1);
////								tarFamily.erase(tfit-1);
////								tarContig.erase(tcit-1);
//								count++;
//								//goto fix;
//								}
//							}						
//						}
//					} 
//				}
//			}
//		}
//		//case4
//		for (rfit = refFamily.end()-1, rcit = refContig.end()-1; rfit != refFamily.begin(); rfit--, rcit--)//for each marker on ref
//		{
//			if (refMarkerNum[abs(*rfit)] == 1)//singleton
//			{
//				if (*rcit == *(rcit-1))//left is same contig 
//				{
//					for (tfit = tarFamily.begin(), tcit = tarContig.begin(); tfit+1 != tarFamily.end(); tfit++, tcit++)//for each marker on tar
//					{
//						if (*tfit == -*rfit && tarMarkerNum[abs(*tfit)] == 1) //different sign & same family & singleton 
//						{
//							if (*(tfit+1) == -*(rfit-1) && *tcit == *(tcit+1))//right is different sign & same family & same contig
//							{
//								if(refMarkerNum[abs(*(rfit-1))] > 1 || tarMarkerNum[abs(*(tfit+1))] > 1)
//								{ 
//								// cout << "4 remove ref:" << *(rfit) << " " << *(rfit-1) << "*rcit:" << *rcit << " " << "*(rcit-1):" << *(rcit-1) << endl;
//								// cout << "4 remove tar:" << *(tfit) << " " << *(tfit+1) << "*tcit:" << *tcit << " " << "*(tcit+1):" << *(tcit+1) << endl;
//								isend = 0;
//								refMarkerNum[abs(*(rfit-1))]--;
//								tarMarkerNum[abs(*(tfit+1))]--;	
//								*(rfit-1) = ++refMarkerMax;
//								if (refMarkerNum.size() == refMarkerMax)
//									refMarkerNum.pb(1);
//								else
//									refMarkerNum[refMarkerMax] = 1;
//
//								*(tfit+1) = -(++tarMarkerMax);
//								if (tarMarkerNum.size() == tarMarkerMax) //warning: Fool Proof
//									tarMarkerNum.pb(1);
//								else
//									tarMarkerNum[tarMarkerMax] = 1;							
////								refFamily.erase(rfit-1);
////								refContig.erase(rcit-1);
////								tarFamily.erase(tfit+1);
////								tarContig.erase(tcit+1);	
//								count++;	
//								//goto fix;
//								} 
//							}
//						}
//					}
//				}
//			}
//		}
		// fix:
		// //case1
		// for (rfit = refFamily.begin(), rcit = refContig.begin(), rtit = refTelomere.begin(); rfit+2 != refFamily.end(); rfit++, rcit++, rtit++)//for each marker on ref
		// {
		// 	if (refMarkerNum[abs(*rfit)] == 1) //singleton
		// 	{	
		// 		//cout << "find singleton on ref: " << *rfit << endl;
		// 		if (*(rcit+1) == *(rcit+2) && *rcit == *(rcit+1)) //right is not telo & same contig 
		// 		{					
		// 			for (tfit = tarFamily.begin(), tcit = tarContig.begin(); tfit+2 != tarFamily.end(); tfit++, tcit++)//for each marker on tar 
		// 			{
		// 				if(*tfit == *rfit && tarMarkerNum[abs(*tfit)] == 1)//same sign & same family & singleton 
		// 				{
		// 					//cout << "find singleton on tar: " << *tfit << endl;
		// 					if(*(tcit+1) == *(tcit+2) && *(tfit+1)== *(rfit+1) && *tcit == *(tcit+1))//rigth is not telo & same sign & family & same contig
		// 					{	
		// 						// cout << "1 remove ref:" << *(rfit) << " " << *(rfit+1) << "*rcit:" << *rcit << " " << "*(rcit+1):" << *(rcit+1) << endl;
		// 						// cout << "1 remove tar:" << *(tfit) << " " << *(tfit+1) << "*tcit:" << *tcit << " " << "*(tcit+1):" << *(tcit+1) << endl;
		// 						isend = 0;
		// 						refMarkerNum[abs(*(rfit+1))]--;
		// 						tarMarkerNum[abs(*(tfit+1))]--;
		// 						refFamily.erase(rfit+1);
		// 						refContig.erase(rcit+1);
		// 						tarFamily.erase(tfit+1);
		// 						tarContig.erase(tcit+1);
		// 						count++;
		// 						goto fix;
		// 					}
		// 				}
		// 			}
		// 		}
		// 	}
		// }
		// //case2
		// for (rfit = refFamily.begin(), rcit = refContig.begin(); rfit+2 != refFamily.end(); rfit++, rcit++)//for each marker on ref
		// {
		// 	if (refMarkerNum[abs(*rfit)] == 1) //singleton
		// 	{	
		// 		if (*(rcit+1) == *(rcit+2) && *rcit == *(rcit+1)) //rigth is not telo & same contig 
		// 		{
		// 			for (tfit = tarFamily.end()-1, tcit = tarContig.end()-1; tfit-1 != tarFamily.begin(); tfit--, tcit--)//for each marker on tar 
		// 			{
						
		// 				if (*tfit == -*rfit && tarMarkerNum[abs(*tfit)] == 1)//different sign & same family & singleton 
		// 				{
							
		// 					if (*(tcit-1) == *(tcit-2) && *(tfit-1) == -*(rfit+1) && *tcit == *(tcit-1))//left is not telo & differnrt sign & same family & same contig
		// 					{
		// 						// cout << "2 remove ref:" << *(rfit) << " " << *(rfit+1) << "*rcit:" << *rcit << " " << "*(rcit+1):" << *(rcit+1) << endl;
		// 						// cout << "2 remove tar:" << *(tfit) << " " << *(tfit-1) << "*tcit:" << *tcit << " " << "*(tcit-1):" << *(tcit-1) << endl;
		// 						isend = 0;								
		// 						refMarkerNum[abs(*(rfit+1))]--;
		// 						tarMarkerNum[abs(*(tfit-1))]--;
		// 						refFamily.erase(rfit+1);
		// 						refContig.erase(rcit+1);
		// 						tarFamily.erase(tfit-1);
		// 						tarContig.erase(tcit-1);
		// 						count++;
		// 						goto fix;
		// 					}
		// 				}
		// 			}
		// 		}
		// 	}
		// }
		// //case3
		// for (rfit = refFamily.end()-1, rcit = refContig.end()-1; rfit-1 != refFamily.begin(); rfit--, rcit--)//for each marker on ref
		// {
		// 	if (refMarkerNum[abs(*rfit)] == 1)//singleton
		// 	{
		// 		if (*(rcit-1) == *(rcit-2) && *rcit == *(rcit-1))//left is not telo & same contig 
		// 		{
		// 			for (tfit = tarFamily.end()-1, tcit = tarContig.end()-1; tfit-1 != tarFamily.begin(); tfit--, tcit--)//for each marker on tar
		// 			{
		// 				if(*tfit == *rfit && tarMarkerNum[abs(*tfit)] == 1)//same sign & same family & singleton 
		// 				{	
		// 					if(*(tcit-1) ==  *(tcit-2) && *(tfit-1) == *(rfit-1) && *tcit == *(tcit-1))//left is not telo & same sign & same family & same contig
		// 					{	
		// 						// cout << "3 remove ref:" << *(rfit) << " " << *(rfit-1) << "*rcit:" << *rcit << " " << "*(rcit-1):" << *(rcit-1) << endl;
		// 						// cout << "3 remove tar:" << *(tfit) << " " << *(tfit-1) << "*tcit:" << *tcit << " " << "*(tcit-1):" << *(tcit-1) << endl;
		// 						isend = 0;
		// 						refMarkerNum[abs(*(rfit-1))]--;
		// 						tarMarkerNum[abs(*(tfit-1))]--;
		// 						refFamily.erase(rfit-1);
		// 						refContig.erase(rcit-1);
		// 						tarFamily.erase(tfit-1);
		// 						tarContig.erase(tcit-1);
		// 						count++;
		// 						goto fix;
		// 					}						
		// 				}
		// 			} 
		// 		}
		// 	}
		// }
		// //case4
		// for (rfit = refFamily.end()-1, rcit = refContig.end()-1; rfit-1 != refFamily.begin(); rfit--, rcit--)//for each marker on ref
		// {
		// 	if (refMarkerNum[abs(*rfit)] == 1)//singleton
		// 	{
		// 		if (*(rcit-1) == *(rcit-2) && *rcit == *(rcit-1))//left is not telo & same contig 
		// 		{
		// 			for (tfit = tarFamily.begin(), tcit = tarContig.begin(); tfit+2 != tarFamily.end(); tfit++, tcit++)//for each marker on tar
		// 			{
		// 				if (*tfit == -*rfit && tarMarkerNum[abs(*tfit)] == 1) //different sign & same family & singleton 
		// 				{
		// 					if (*(tcit+1) == *(tcit+2) && *(tfit+1) == -*(rfit-1) && *tcit == *(tcit+1))//right is not telo & different sign & same family & same contig
		// 					{
		// 						// cout << "4 remove ref:" << *(rfit) << " " << *(rfit-1) << "*rcit:" << *rcit << " " << "*(rcit-1):" << *(rcit-1) << endl;
		// 						// cout << "4 remove tar:" << *(tfit) << " " << *(tfit+1) << "*tcit:" << *tcit << " " << "*(tcit+1):" << *(tcit+1) << endl;
		// 						isend = 0;
		// 						refMarkerNum[abs(*(rfit-1))]--;
		// 						tarMarkerNum[abs(*(tfit+1))]--;								
		// 						refFamily.erase(rfit-1);
		// 						refContig.erase(rcit-1);
		// 						tarFamily.erase(tfit+1);
		// 						tarContig.erase(tcit+1);	
		// 						count++;	
		// 						goto fix;
		// 					}
		// 				}	
		// 			}
		// 		}
		// 	}
		// }
		//find unique marker on ref
	// 	rfit = refFamily.begin(), rcit = refContig.begin();
	// 	while(rfit != refFamily.end())
	// 	{
	// 		if(tarMarkerNum[abs(*rfit)] == 0)
	// 		{	
	// 			cout << tarMarkerNum[abs(*rfit)] << "find unique marker on ref:" << *rfit << endl;
	// 			isend = 0;
	// 			refMarkerNum[abs(*rfit)]--;
	// 			refFamily.erase(rfit);
	// 			refContig.erase(rcit);
	// 		}
	// 		else
	// 		{
	// 			rfit++, rcit++;
	// 		}
	// 	}
	// 	//find unique marker on tar
	// 	tfit = tarFamily.begin(), tcit = tarContig.begin();
	// 	while(tfit != tarFamily.end())
	// 	{
	// 		if(refMarkerNum[abs(*tfit)] == 0)
	// 		{
	// 			cout << refMarkerNum[abs(*tfit)] << "find unique marker on tar:" << *tfit << endl;
	// 			isend = 0;
	// 			tarMarkerNum[abs(*tfit)]--;
	// 			tarFamily.erase(tfit);
	// 			tarContig.erase(tcit);
	// 		}
	// 		else
	// 		{
	// 			tfit++, tcit++;
	// 		}
	// 	}
	// }
	// rfit = refFamily.begin();
	// rcit = refContig.begin();
	// while(rfit != refFamily.end() - 1)
	// {
	// 	if(refMarkerNum[abs(*rfit)] == 1 && refMarkerNum[abs(*(rfit+1))] == 1 && *rcit == *(rcit+1))
	// 	{
	// 		tfit = tarFamily.begin();
	// 		tcit = tarContig.begin();
	// 		for(;tfit != tarFamily.end() - 1; tfit++, tcit++)
	// 		{
	// 			if(*tfit == *rfit && *(tfit+1) == *(rfit+1) && tarMarkerNum[abs(*tfit)] == 1 && tarMarkerNum[abs(*(tfit+1))] == 1 && *tcit == *(tcit+1))
	// 			{
	// 				refMarkerNum[abs(*(rfit+1))]--;
	// 				tarMarkerNum[abs(*(tfit+1))]--;
	// 				refFamily.erase(rfit+1);
	// 				refContig.erase(rcit+1);
	// 				tarFamily.erase(tfit+1);
	// 				tarContig.erase(tcit+1);
	// 				count++;
	// 				break;
	// 			}
	// 		}
	// 		rfit++, rcit++;
	// 	}
	// 	else
	// 	{
	// 		rfit++, rcit++;
	// 	}
	// }
	cout << "# of ref contigs:" << refContigSize << endl;
	cout << "# of tar contigs:" << tarContigSize << endl;
	cout << "# of tar merges:" << tarNumberMerge << endl;
	cout << "# of ref merges:" << refNumberMerge << endl;

	string out_dir(argc < 4 ? "output" : argv[3]);
	ofstream fout(out_dir+"/ref_spd3.all");

	int tmpid = 1;
	for (int i = 0; i < refFamily.size(); i++)
	{
		fout << tmpid++ << " " << refFamily[i] << " " << refContig[i] << " " << 1 << endl;
		// else
		// 	cout << "find unique marker on ref:" << refFamily[i] << endl;
	}
	fout.close();

	fout.open(out_dir+"/tar_spd3.all");
	tmpid = 1;
	for (int i = 0; i < tarFamily.size(); i++)
	{
		fout << tmpid++ << " " << tarFamily[i] << " " << tarContig[i] << " " << 1 << endl;
		// else
		// 	cout << "find unique marker on tar:" << tarFamily[i] << endl;
	}
	fout.close();

	fout.open(out_dir+"/tarContigMerge.all");
	for(auto it : tarContigMerge)
	{
		fout << it.first << endl;
		for(auto it2 : it.second)
		{
			fout << it2.first << " " << it2.second << endl;
		}
		fout << "end" << endl;
	}
	fout.close();

	fout.open(out_dir+"/refContigMerge.all");
	for(auto it : refContigMerge)
	{
		fout << it.first << endl;
		for(auto it2 : it.second)
		{
			fout << it2.first << " " << it2.second << endl;
		}
		fout << "end" << endl;
	}
	fout.close();
	// cout << "ref:" << endl;
	// for(int i=0; i < refFamily.size(); i++)
	// {
	// 	if(tarMarkerNum[abs(refFamily[i])] != 0)
	// 	{
	// 		cout << refFamily[i] << " " << refMarkerNum[abs(refFamily[i])] << " " << tarMarkerNum[abs(refFamily[i])] << " " << refContig[i] << endl;
	// 	}
	// }
	// cout << "tar:" << endl;
	// for(int i=0; i < tarFamily.size(); i++)
	// {
	// 	if(refMarkerNum[abs(tarFamily[i])] != 0)
	// 	{
	// 		cout << tarFamily[i] << " " << tarMarkerNum[abs(tarFamily[i])] << " " << tarMarkerNum[abs(tarFamily[i])] << " " << tarContig[i] << endl;
	// 	}
	// }
	// set<int> ref;
	// for(int i=0; i < refFamily.size(); i++)
	// {
	// 	ref.insert(abs(refFamily[i]));
	// }
	// cout << "ref" << endl;
	// for(auto it = ref.begin(); it != ref.end(); it++)
	// {
	// 	cout << *it << " " << refMarkerNum[abs(*it)] << " " << tarMarkerNum[abs(*it)] << endl;
	// }
	// set<int> tar;
	// for(int i=0; i < tarFamily.size(); i++)
	// {
	// 	tar.insert(abs(tarFamily[i]));
	// }
	// cout << "tar" << endl;
	// for(auto it = tar.begin(); it != tar.end(); it++)
	// {
	// 	cout << *it << " " << tarMarkerNum[abs(*it)] << " " << refMarkerNum[abs(*it)] << endl;
	// }
	// showNumberOfMarkers();
	// int sinOnRef = 0;
	// int dupOnRef = 0;
	// for(int i = 0; i < refFamily.size(); i++)
	// {
	// 	if(tarMarkerNum[abs(refFamily[i])] != 0 && refMarkerNum[abs(refFamily[i])] == 1)
	// 	{
	// 		sinOnRef++;
	// 	}
	// 	if(tarMarkerNum[abs(refFamily[i])] != 0 && refMarkerNum[abs(refFamily[i])] != 1)
	// 	{
	// 		dupOnRef++;
	// 	}
	// }
	// int sinOnTar = 0;
	// int dupOnTar = 0;
	// for(int i = 0; i < tarFamily.size(); i++)
	// {
	// 	if(refMarkerNum[abs(tarFamily[i])] != 0 && tarMarkerNum[abs(tarFamily[i])] == 1)
	// 	{
	// 		sinOnTar++;
	// 	}
	// 	if(refMarkerNum[abs(tarFamily[i])] != 0 && tarMarkerNum[abs(tarFamily[i])] != 1)
	// 	{
	// 		dupOnTar++;
	// 	}
	// }
	// cout << "#sinOnRef: " << sinOnRef << " #dupOnRef: " << dupOnRef << endl;
	// cout << "#sinOnTar: " << sinOnTar << " #dupOnTar: " << dupOnTar << endl;
}
