#include<fstream>
#include<iostream>
#include<vector>
#include<set>
#define pb push_back
#define abs(a) (a > 0 ? a : -a)
using namespace std;
struct Marker{
	int id, family, absFamily;
	string contig;
	Marker(int id_, int family_, string contig_) {
		id = id_;
		family = family_;
		absFamily = abs(family_);
		contig = contig_;
	}
	void setFamily(int newFamily) {
		family = newFamily;
		absFamily = abs(newFamily);
	}
	void show() {
		cout<<"Marker("<<id<<", "<<family<<", "<<contig<<")"<<endl;
	}
};
int markerReorder(string ref_file, string tar_file)
{
	cout<<"===== ===== ===== markerReorder ===== ===== ====="<<endl;
	ifstream fin(ref_file);
	int id,fam,tmp;
	string contig;
	set<int> fa,fb;
	vector<Marker> ref,tar;

	while(fin>>id>>fam>>contig>>tmp)
	{
		fa.insert(abs(fam));
		Marker m(id,fam,contig);
		ref.pb(m);
	}
	for(int i:fa)
		cout<<i<<" ";
	cout<<endl;
	fin.close();
	fin.open(tar_file);

	while(fin>>id>>fam>>contig>>tmp)
	{
		fb.insert(abs(fam));
		Marker m(id,fam,contig);
		tar.pb(m);
	}
	cout<<"-----------------------"<<endl;
	for(auto it=fa.begin(),it2=fb.begin();it!=fa.begin(),it2!=fb.begin();it++,it2++)
	{
		if(*it!=*it2)
			cout<<"g"<<endl;

	}
	auto it = fa.end();
	it--;
	cout<<"it val:"<<*it<<endl;
	vector<int> reorder(*it+1,0);
	cout<<"rsize:"<<reorder.size()<<endl;
	int tid = 1;
	for(int n:fa)
		reorder[n]=tid++;
	cout<<fa.size()<<endl;
	fin.close();

	//ref_file.replace(ref_file.find(".all"), 4, "_mr.all");
	ofstream fout(ref_file);
	int tmpid=1;
	cout<<"p1"<<endl;
	for(auto m:ref)
	{
		if(m.family<0)
			fout<<tmpid++<<" "<<-reorder[-m.family]<<" "<<m.contig<<" "<<1<<endl;
		else
			fout<<tmpid++<<" "<<reorder[m.family]<<" "<<m.contig<<" "<<1<<endl;
	}
	fout.close();

	//tar_file.replace(tar_file.find(".all"), 4, "_mr.all");
	fout.open(tar_file);
	cout<<"p2"<<endl;
	tmpid=1;
	for(auto m:tar)
	{
		if(m.family<0)
			fout<<tmpid++<<" "<<-reorder[-m.family]<<" "<<m.contig<<" "<<1<<endl;
		else
			fout<<tmpid++<<" "<<reorder[m.family]<<" "<<m.contig<<" "<<1<<endl;
	}
	cout<<"p3"<<endl;
	//fout.close();
	return 0;
}
