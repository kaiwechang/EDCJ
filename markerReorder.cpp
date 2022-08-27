#include<fstream>
#include<iostream>
#include<vector>
#include<set>
#define pb push_back
#define abs(x) (x>0?x:-x)
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
int main(int argc, char *argv[])
{
	ifstream fin(argv[1]);
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
	fin.open(argv[2]);

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

	string str(argv[1]);
	str.replace(str.find(".all"), 4, "_mr.all");
	ofstream fout(str);
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

	str = string(argv[2]);
	str.replace(str.find(".all"), 4, "_mr.all");
	fout.open(str);
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
