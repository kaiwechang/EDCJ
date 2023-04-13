#include <iostream>
#include <fstream>
#include <algorithm>
#include <string>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <functional>
#include <random>

using namespace std;

//generate random position within the genome
int Position_rand(int len_gen)
{
	int pos_rand;
	random_device rd;

	mt19937 gen = mt19937(rd());
	uniform_real_distribution<> dis(0, 1000000);
	auto randfun = bind(dis, gen);

	pos_rand = int(randfun()) % len_gen;

	return pos_rand;
}

//build up initial genome
vector<int> getInitialGenome(int num_marker)
{
	vector<int> genome;
	
	for(int i = 0; i < num_marker; i++)
	{
		genome.push_back(i+1);
	}

	return genome;
}

//output .all file
int Output_ALL(vector<int> genome, vector<int> pos_contig, string filename_all)
{
	int counter_contig;
	int len_gen;
	fstream file;

	len_gen = genome.size();
	counter_contig = 0;

	file.open(filename_all, ios::out);

	if(pos_contig.size() == 0)
	{
		for(int i = 0; i < len_gen; i ++)
		{
			file << i + 1 << " " << genome[i] << " contig_" << counter_contig + 1 << " 1" << endl;
		}
	}
	else
	{
		for(int i = 0; i < len_gen; i ++)
		{		
			if((i >= pos_contig[counter_contig]) && (counter_contig < pos_contig.size()))
			{
				counter_contig++;
			}
		
			file << i + 1 << " " << genome[i] << " contig_" << counter_contig + 1 << " 1" << endl;
		}
	}

	file.close();

	return 0;
}

//output genome to log file
int Output_Genome(vector<int> genome, string filename_log)
{	
	int len_gen = genome.size();
	fstream file_log;

	file_log.open(filename_log, ios::out | ios::app);

	for(int i = 0; i < len_gen; i++)
	{
		if(i == len_gen - 1)
		{
			file_log << genome[i] << endl;

			file_log.close();

			return 0;
		}

		file_log << genome[i] << ", ";
	}

	file_log.close();

	return 0;
}

//create answer to scaffold
int CreateAnswer(int num_contig, string filename_ans)
{
	fstream file;

	file.open(filename_ans, ios::out);

	file << "> Scaffold_1" << endl;

	for(int i = 0; i < num_contig; i++)
	{
		file << "contig_" << i + 1 << " 0" << endl;
	}

	file.close();

	return 0;
}

//check if the position is unique in the vector
bool isUnique(vector<int> pos_contig, int pos)
{
	for(int i = 0; i < pos_contig.size(); i++)
	{
		if(pos_contig[i] == pos)
		{
			return false;
		}
	}

	return true;
}

vector<int> getPosition_Contigs(int len_gen, int num_contig)
{
	vector<int> pos_contig;
	int pos_rand;

	if(num_contig == 1)
	{
		return pos_contig;
	}
	else
	{
		for(int i = 0; i < num_contig - 1; i++)  //cut (n - 1) times to create n contigs
		{
			do
			{
				pos_rand = Position_rand(len_gen - 1);  //for n elemnets, there are only (n - 1) positions could be cut
			}while(!isUnique(pos_contig, pos_rand + 1));  //check if pos_rand is unique in pos_contig
		
			pos_contig.push_back(pos_rand + 1);  //deviate by 1, since we cannot cut the position before the first element
		}
	}

	sort(pos_contig.begin(), pos_contig.end());
	
	return pos_contig;
}

int Duplicate(vector<int> &genome, int len_dup, string filename_log)
{
	int len_gen = genome.size();
	int pos_dup = Position_rand(len_gen - len_dup + 1);  //cannot form a duplicate by less than len_dup elements
	int pos_ins = Position_rand(len_gen + 1);  //for n elements, there are (n +1) positions could be inserted
	vector<int> duplicate;
	fstream file_log;

	file_log.open(filename_log, ios::out | ios::app);
	
	file_log << "Start position: " << pos_dup + 1 << "    ";
	file_log << "End position: " << pos_dup + len_dup << "    ";  //(pos_dup + len_dup - 1) + 1
	file_log << "Insert position: " << pos_ins + 1 << endl;
	
	file_log.close();

	//build a copy of the duplicate
	for(int i = 0; i < len_dup; i++)
	{
		duplicate.push_back(genome[pos_dup + i]);
	}

	//insert the duplicate into pos_ins
	for(int i = 0; i < len_dup; i++)
	{
		genome.insert(genome.begin() + pos_ins + i, duplicate[i]);
	}

	Output_Genome(genome, filename_log);

	return 0;
}

int Inverse(vector<int> &genome, string filename_log)
{
	int len_gen = genome.size();
	int pos_1 = Position_rand(len_gen);
	int pos_2 = Position_rand(len_gen);
	fstream file_log;

	file_log.open(filename_log, ios::out | ios::app);

	//to make sure pos_1 is the smaller one
	if(pos_1 > pos_2)
	{
		swap(pos_1, pos_2);
	}

	file_log << "Start position: " << pos_1 + 1 << "    ";
	file_log << "End position: " << pos_2 + 1 << endl;

	file_log.close();

	//reverse the sequence of the markers between pos_1 and pos_2
	while(pos_1 <= pos_2)
	{
		if(pos_1 == pos_2)
		{
			genome[pos_1] = genome[pos_1] * (-1);
			pos_1++;
			pos_2--;
		}else
		{
			genome[pos_1] = genome[pos_1] * (-1);
			genome[pos_2] = genome[pos_2] * (-1);
			swap(genome[pos_1], genome[pos_2]);
			pos_1++;
			pos_2--;
		}
	}

	Output_Genome(genome, filename_log);

	return 0;
}

int Simulator(int num_marker, int prob_event, int len_dup, int num_evol, int num_contig, string filename, string filename_log)
{
	vector<int> genome;
	vector<int> pos_contig;
	int len_gen;
	int counter_inv = 0;
	int counter_dup = 0;
	int num_evol_inv = int(num_evol * prob_event / 100); 
	fstream file_result;
	fstream file_log;	

	//build up initial genome
	genome = getInitialGenome(num_marker);

	file_log.open(filename_log, ios::out);  //without ios::app, clear the content of the original file

	file_log << "initial genome:" << endl;

	file_log.close();

	Output_Genome(genome, filename_log);

	file_log.open(filename_log, ios::out | ios::app);

	file_log << endl;

	file_log.close();

	//start evolution
	file_log.open(filename_log, ios::out | ios::app);

	for(int i = 0; i < num_evol; i++)
	{
		file_log << "Step " << i + 1 << ": ";

		int prob = rand() % 100;
		
		if((prob + 1) <= prob_event)  //inverse occurs with probability of (prob_event)%
		{	
			if(counter_inv >= num_evol_inv)
			{
				file_log << "duplicate" << endl;

				Duplicate(genome, len_dup, filename_log);

				counter_dup ++;
			}
			else
			{
				file_log << "inverse" << endl;

				Inverse(genome, filename_log);

				counter_inv ++;
			}
		}
		else if((prob + 1) > prob_event)  //duplicate occurs with probability of (100 - prob_event)%
		{	
			if(counter_dup >= (num_evol - num_evol_inv))
			{
				file_log << "inverse" << endl;

				Inverse(genome, filename_log);

				counter_inv ++;
			}
			else
			{
				file_log << "duplicate" << endl;

				Duplicate(genome, len_dup, filename_log);

				counter_dup ++;
			}
		}

		file_log << endl;
	}

	//determine where to cut the genome
	len_gen = genome.size();
	pos_contig = getPosition_Contigs(len_gen, num_contig);

	//output to .all file
	Output_ALL(genome, pos_contig, filename);
		
	return 0;
}

int CreateAncestor(int num_marker, int num_contig, string filename_anc)
{
	vector<int> genome;
	vector<int> pos_contig;
	int len_gen;
	
	//build up ancestor genome
	genome = getInitialGenome(num_marker);

	//determine where to cut the genome
	len_gen = genome.size();
	pos_contig = getPosition_Contigs(len_gen, num_contig);

	//output to .all file
	Output_ALL(genome, pos_contig, filename_anc);

	return 0;
}

int main(int argc, char *argv[])
{
	if (argc < 8) {
		cout << "[error] Usage:\n>>> Simulator <# initial markers> <inverse rate> <duplicate length> <# evolutions> <# ref contigs> <# tar contigs> <output_dir>\n";
		return 0;
	}

	int num_marker = stoi(argv[1]);  //number of initial markers
	int prob_event = stoi(argv[2]);  //probability of inverse event, probability of duplicate event is (1 - prob_event)
	int len_dup = stoi(argv[3]);  //length of duplicate segment
	int num_evo = stoi(argv[4]);  //number of evolution
	int num_contig_ref = stoi(argv[5]);  //number of reference contigs
	int num_contig_que = stoi(argv[6]);  //number of target contigs

	string out_dir(argv[7]);
	//string filename_anc = "reference_test.all";			// for using ancestor genome as reference
	string filename_ref		= out_dir + "/reference.all";	// for using evloved genome as reference
	string filename_que		= out_dir + "/target.all";		// used to be called query.all
	string filename_ref_log = out_dir + "/reference_process";
	string filename_que_log = out_dir + "/target_process";
	string filename_ans		= out_dir + "/answerToAll";
	
	//int num_evol = float(num_marker) / (((100 - float(prob_event)) / 100) * float(len_dup));  //number of evolutions = N / ((1 - P) * L)

	//srand(time(NULL));
		
	//CreateAncestor(num_marker, num_contig_ref, filename_anc);  //create ancestor genome as reference

	Simulator(num_marker, prob_event, len_dup, num_evo, num_contig_ref, filename_ref, filename_ref_log);  //create evolved genome as reference
	Simulator(num_marker, prob_event, len_dup, num_evo, num_contig_que, filename_que, filename_que_log);

	//create answer to scaffold
	CreateAnswer(num_contig_que, filename_ans);
		
	return 0;
}
