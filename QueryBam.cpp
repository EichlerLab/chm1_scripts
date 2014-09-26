#include "BitNucVector.h"
#include "sam.h"
#include "stdlib.h"
#include <vector>
#include <string>
#include <set>
#include <iostream>
#include <fstream>
#include <map>
#include <sstream>
#include <numeric>
#include <pthread.h>
#include <semaphore.h>
#include <fcntl.h>
                         //0, 1, 2, 3, 4, 5, 6, 7, 8
const int BamToTwoBit[] = {0, 0, 1, 0, 2, 0, 0, 0, 3};
const char BamToAscii[] = {0,'A','C',0,'G',0,0,0,'T',0,0,0,0,0,0,0,0,0,0,0,'N'};

sem_t *semreader;
sem_t *semcount;



int InitializeTupleFromBam(uint8_t *seq, int pos, int length, BitNucVector &tuple) {
	int i;
	tuple.Clear();

	bool isInitialized = false;
	i = 0;
	while (isInitialized == false and pos < length - tuple.k + 1) {
		for (i = 0; i < tuple.k; i++) {
			int nuc = bam1_seqi(seq, i+pos);
			if (nuc == 16) {
				pos = i;
				tuple.Clear();
				break;
			}
			else {
				tuple.Set(i, BamToTwoBit[nuc]);
			}
		}
		if (i == tuple.k) {
			isInitialized = true;
		}
	}
	if (isInitialized) {
		return pos;
	}
	else {
		return -1;
	}
}

int AdvanceNucleotide(uint8_t* seq, int pos, int length, BitNucVector &tuple) {
	int nuc = bam1_seqi(seq, pos+tuple.k);
	if (pos +tuple.k >= length) {
		return pos + tuple.k;
	}
	if (nuc != 16) {
		tuple.ShiftOneNucleotide(BamToTwoBit[nuc]);
		// 
		// this nucleotide is fine, next iteration will not reset tuple.
    return 0;
	}
	else {
		// signal that must start search past 'N'
		return pos + tuple.k + 1;
	}
}

class Counts {
public:
	samfile_t *in;  	
	map<BitNucVector, int> *tupleToIndex;
	vector<int> *queryCount;
	int *readIndex ;
};

string GetSeq(bam1_t *b) {
	int lSeq = b->core.l_qseq;
	string res;
	res.resize(lSeq);
	int i;
	uint8_t *seq = bam1_seq(b);
	for (i = 0; i < lSeq; i++) {
		res[i] = BamToAscii[bam1_seqi(seq,i)];
	}
	return res;
}

void CountWords(void *data) {
	map<BitNucVector, int> &tupleToIndex =  *((Counts*)data)->tupleToIndex;
	vector<int> &queryCount = *((Counts*)data)->queryCount;
	int &readNumber = *((Counts*)data)->readIndex;
	samfile_t *in = ((Counts*)data)->in;
	bam1_t *b = bam_init1();
	map<BitNucVector, int>::iterator searchResult;

	while (true) {
		int retval;
		retval = sem_wait(semreader);		
		retval = bam_read1(in->x.bam, b);
		if (retval <= 0) {
			sem_post(semreader);
			return;
		}
		sem_post(semreader);

		uint8_t *seq = bam1_seq(b);
		int seqLength = b->core.l_qseq;
		string str = GetSeq(b);
		int seqPos = 0;
		int res;
		BitNucVector tuple;
		res = InitializeTupleFromBam(seq, seqPos, seqLength, tuple);
		if (res >= 0) {
			seqPos = res;
			do {
				//
				// Store result.
				//
				searchResult = tupleToIndex.find(tuple);
				if (searchResult != tupleToIndex.end()) {
					sem_wait(semcount);
					queryCount[(*searchResult).second] += 1;				
					//					cerr << (*searchResult).first.ToString() << " " << readNumber << " " << (*searchResult).second << "\t" << queryCount[(*searchResult).second] << endl;
					sem_post(semcount);

				}
				
				// 
				// Move to the next nucleotide;
				//
				int res = AdvanceNucleotide(seq, seqPos, seqLength, tuple);
				if (res == 0) {
					seqPos++;
				}
				else {
					seqPos = InitializeTupleFromBam(seq, seqPos, seqLength, tuple);
				}
			} while (seqPos < seqLength - tuple.k );
		}
		sem_wait(semcount);
		if (readNumber % 1000000 == 0) {
			cerr << readNumber << endl;
		}
		++readNumber;
		sem_post(semcount);

	}



}

int main(int argc, char* argv[]) {
	string queryFileName, bamFileName, outputFileName;
	
	if (argc < 4) {
		cout << "Usage: queryBam queryTable file.bam outputFile"  << endl;
		cout << "  queryTable is in the format name seq. " << endl;
		cout << "  Each query sequence must be less than 32 nt long." << endl;
		exit(0);
	}
	queryFileName = argv[1];
	bamFileName = argv[2];
	outputFileName = argv[3];
	ofstream outFile(outputFileName.c_str());
	ifstream queryFile(queryFileName.c_str());
	//
	// Make sure there are no duplicate keys.
	//
	if (not queryFile) {
		cout << "Could not open " << queryFileName << endl; exit(1);
	}
	map<BitNucVector, int> queryList;
	set<string> uniqueQueries;
	set<string> skipQueries;
	while (queryFile) {
		string seq, tuple;
		if (! (queryFile >> seq >> tuple)  ) { break; }
		if (uniqueQueries.find(seq) != uniqueQueries.end()) {
			skipQueries.insert(seq);
		}
	}

	queryFile.close();
	queryFile.open(queryFileName.c_str());
	
	vector<int> queryCount;
	vector<string> queryNames;
	vector<BitNucVector> queryTuples;
	map<BitNucVector, int> tupleToIndex;

	int index = 0;
	while (queryFile) {
		string seq, tupleStr;
		if (! (queryFile >> seq >> tupleStr)  ) { break; }
		if (skipQueries.find(seq) != skipQueries.end()) {
			continue;
		}
		else {
			BitNucVector tuple, tuplerc;
			tuple.k = tupleStr.size();
			tuple.InitializeFromString((unsigned char*) tupleStr.c_str(), tupleStr.size());
			queryCount.push_back(0);
			queryNames.push_back(seq);
			queryTuples.push_back(tuple);
			tupleToIndex[tuple] = index;
			
			tuple.SetReverseComplement(tuplerc);
			tupleToIndex[tuplerc] = index;
			++index;
		}
		
	}
	samfile_t *in;  
	in = samopen(bamFileName.c_str(), "rb", 0);

	//  idx = bam_index_load(bamFileName.c_str());
	bam1_t *b = bam_init1();
	map<BitNucVector, int>::iterator searchResult;
	int readNumber = 0;
	
	Counts counts;
	counts.in = in;
	counts.tupleToIndex = &tupleToIndex;
	counts.queryCount = &queryCount;
	counts.readIndex = &readNumber;

	semreader     = sem_open("/semreader",     O_CREAT, 0644, 1);
	semcount      = sem_open("/semcount",     O_CREAT, 0644, 1);
	sem_init(semreader, 1, 1);
	sem_init(semcount, 1, 1);

	int nproc = 8;
	pthread_attr_t *threadAttr = new pthread_attr_t[nproc];
	int t;	

	for (t = 0; t < nproc; t++ ){
		pthread_attr_init(&threadAttr[t]);
	}
	pthread_t *threads = new pthread_t[nproc];

	for (t = 0; t < nproc; t++) {
		pthread_create(&threads[t], &threadAttr[t], (void* (*)(void*))CountWords, &counts);
	}

	for (t = 0; t < nproc; t++) {
		pthread_join(threads[t], NULL);
	}

	int i;																			 
	for (i = 0; i < queryNames.size(); i++) {
		outFile << queryNames[i] << "\t" << queryTuples[i].ToString() << "\t" << queryCount[i] << endl;
	}
	outFile.close();
}


