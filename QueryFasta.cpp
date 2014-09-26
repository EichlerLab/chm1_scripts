#include "BitNucVector.h"

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
#include <errno.h>
#include <time.h>
                         //0, 1, 2, 3, 4, 5, 6, 7, 8
const int BamToTwoBit[] = {0, 0, 1, 0, 2, 0, 0, 0, 3};
const char BamToAscii[] = {0,'A','C',0,'G',0,0,0,'T',0,0,0,0,0,0,0,0,0,0,0,'N'};

sem_t *semreader;
sem_t *semcount;
sem_t *semfastqwriter;

void rand_str(char *dest, size_t length) {
    char charset[] = "0123456789"
                     "abcdefghijklmnopqrstuvwxyz"
                     "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
		
    while (length-- > 0) {
        size_t index = (double) rand() / RAND_MAX * (sizeof charset - 1);
        *dest++ = charset[index];
    }
    *dest = '\0';
}

int AdvanceNucleotide(string &seq, int pos, int length, BitNucVector &tuple) {

	if (pos +tuple.k >= length) {
		return pos + tuple.k;
	}
	if (seq[pos] != 'N') {
		tuple.ShiftOneNucleotide(TwoBit[seq[pos+tuple.k]]);
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
	ifstream *in, *pairIn;
	ofstream *readsOut;
	bool isFastq;
	bool noDoubleCount;
	map<BitNucVector, vector<int> > *tupleToIndex;
	vector<int> *queryCount;
	int *readIndex ;
	int maxNReads;
	Counts() {
		pairIn = NULL;
		readsOut = NULL;
	}
};

int readNumber = 0;
int nHits = 0;

void CountWords(void *data) {
	map<BitNucVector, vector<int> > &tupleToIndex =  *((Counts*)data)->tupleToIndex;
	vector<int> &queryCount = *((Counts*)data)->queryCount;
	int &readNumber = *((Counts*)data)->readIndex;
	ifstream *in = ((Counts*)data)->in;
	bool isFastq = ((Counts*)data)->isFastq;
	int numSeq = 1;
	int maxNReads = ((Counts*)data)->maxNReads;
	
	ifstream *pairedIn = ((Counts*)data)->pairIn;
	ofstream *readsOut = ((Counts*)data)->readsOut;
	if (pairedIn != NULL) {
		numSeq = 2;
	}
	bool noDoubleCount = ((Counts*)data)->noDoubleCount;

	map<BitNucVector, vector<int> >::iterator searchResult;
	
	while (true) {
		int retval;
		//
		// Make sure that only one thread at a time enters the reading section
		//
		retval = sem_wait(semreader);		
		string title, seq, qualtitle, qual;
		string mptitle, mpseq, mpqualtitle, mpqual;
		vector<string> seqs(2);
		vector<string> quals(2);
		vector<string> titles(2);
		if (isFastq) {
			//			cout << "reading fastq!!" << endl;
			getline(*in, titles[0]);
			(*in) >> seqs[0];
			(*in).get();
			getline(*in, qualtitle);
			(*in) >> quals[0];
			(*in).get();
			if (numSeq == 2) {
				getline(*pairedIn, titles[1]);
				(*pairedIn) >> seqs[1];
				(*pairedIn).get();
				getline(*pairedIn, mpqualtitle);
				(*pairedIn) >> quals[1];
				(*pairedIn).get();
			}
		}
		else {
			if (numSeq == 2) {
				cout << "ERROR, paired reads only work with fastq for now!" << endl;
				exit(1);
			}
			getline(*in, title);
			while (in->peek() != '>' and in) {
				char p = in->peek();
				//				cout << "got peek: "<< p << endl;
				string line;
				(*in) >> line;
				if (line == "") {
					break;
				}
				assert(line[0] != '>');
				seq += line;
				//				cout << "seq tmp: " << seq << endl;
				if (in->peek() == '\n') {
					in->get();
				}
			}

		}
				
		if (readNumber % 1000000 == 0 or (isFastq == false and readNumber % 10000 == 0)) {
			cerr << readNumber << "\t" << nHits << endl;
		}
		++readNumber;

		if ((isFastq == true and quals[0] == "") or (isFastq == false and (titles[0] == "" or seqs[0] == "")) or (maxNReads > 0 and readNumber >= maxNReads)) {
			cout << "qual: " << quals[0] << endl << " title " << titles[0] << " qualtitle: " << qualtitle << " seq " << seqs[0] << endl;
			cout << "read index " << readNumber << endl;
			sem_post(semreader);
			return;
		}

		//
		// All done reading, release the thread lock.
		//
		sem_post(semreader);

		
		int seqi;
		bool stop = false;
		bool foundMatch = false;
		for (seqi = 0; seqi < numSeq and stop == false; ++seqi) {
			int seqPos = 0;
			int res;
			BitNucVector tuple;
			bool initializeTuple = true;
			seq = seqs[seqi];
			if (seq.size() < tuple.k) {
				continue;
			}
			do {
				//
				// Store result.
				//
				if (initializeTuple) {
					while (seqPos < seq.size() - tuple.k) {
						if (tuple.InitializeFromString((unsigned char*) &seq[seqPos], tuple.k)) {
							break;
						}
						else {
							seqPos++;
						}
					}
					initializeTuple = false;
				}

				if (seqPos >= seq.size() - tuple.k) {
					break;
				}

				searchResult = tupleToIndex.find(tuple);
				if (searchResult != tupleToIndex.end()) {
					sem_wait(semcount);
					int i;
					for (i = 0; i < (*searchResult).second.size(); i++) {
						queryCount[(*searchResult).second[i]] += 1;				
						++nHits;
						foundMatch = true;
						if (noDoubleCount) {
							seqPos = seq.size();
							stop = true;
						}
					}
					sem_post(semcount);
				}
				
				// 
				// Move to the next nucleotide;
				//
				int res = AdvanceNucleotide(seq, seqPos, seq.size(), tuple);
				if (res == 0) {
					seqPos++;
				}
				else {
					initializeTuple = true;
				}

			
			} while (seqPos < seq.size() - tuple.k);
		}
		if (foundMatch and readsOut != NULL) {
			sem_wait(semfastqwriter);
			for (seqi = 0; seqi < numSeq; seqi++) {
				*readsOut << titles[seqi] << endl;
				*readsOut << seqs[seqi] << endl;
				*readsOut << qualtitle << endl;
				*readsOut << quals[seqi] << endl;
			}
			sem_post(semfastqwriter);			
		}
	}
}

int main(int argc, char* argv[]) {
	string queryFileName, fastqFileName, outputFileName;
	
	if (argc < 4) {
		cout << "Usage: queryFasta queryTable file.fastq outputFile [-single] [-nproc n] [-fasta f] [-pair paired reads file] [-matches output.fastq]"  << endl; 
		cout << "  queryTable is in the format name seq. " << endl;
		cout << "  Each query sequence must be less than 32 nt long." << endl;
		exit(0);
	}
	queryFileName  = argv[1];
	fastqFileName  = argv[2];
	outputFileName = argv[3];
	int nproc = 8;
	bool forceFasta = false;
	bool noDoubleCount = false;
	string queryMPFileName = "";
	string matchedReadsFileName = "";
	int maxNReads = 0;
	if (argc > 4) {
		int argi = 4;
		while (argi < argc) {
			if (strcmp(argv[argi], "-nproc") == 0) {
				nproc = atoi(argv[++argi]);
			}
			else if (strcmp(argv[argi], "-fasta") == 0) {
				forceFasta = true;
			}
			else if (strcmp(argv[argi], "-single") == 0) {
				noDoubleCount = true;
			}
			else if (strcmp(argv[argi], "-pair") == 0) {
				queryMPFileName = argv[++argi];
			}
			else if (strcmp(argv[argi], "-matches") == 0) {
				matchedReadsFileName = argv[++argi];
			}
			else if (strcmp(argv[argi], "-nreads") == 0) {
				maxNReads = atoi(argv[++argi]);
			}
			++argi;
		}
	}

	ofstream outFile(outputFileName.c_str());
	ifstream queryFile(queryFileName.c_str());
	ifstream fastqFile(fastqFileName.c_str());
	
	ifstream queryMPFile;
	if (queryMPFileName != "") {
		queryMPFile.open(queryMPFileName.c_str());
	}
	
	ofstream readsOut;
	if (matchedReadsFileName != "") {
		readsOut.open(matchedReadsFileName.c_str());
	}

	bool isFastq;
	if (forceFasta) {
		isFastq = false;
	}
	else {
		if (fastqFileName.find(".fasta") != fastqFileName.npos) {
			isFastq = false;
		}
		else {
			isFastq = true;
		}
	}
	//
	// Make sure there are no duplicate keys.
	//
	if (not queryFile) {
		cout << "Could not open " << queryFileName << endl; exit(1);
	}
	map<BitNucVector, int> queryList;

	
	vector<int> queryCount;
	vector<string> queryNames;
	vector<int> queryPositions;
	vector<BitNucVector> queryTuples;
	map<BitNucVector, vector<int> > tupleToIndex;
	bool kIsInitialized = false;
	int index = 0;
	while (queryFile) {
		string seq, tupleStr;
		int pos;
		if (! (queryFile >> seq >> tupleStr >> pos)  ) { break; }


		BitNucVector tuple, tuplerc;
		if (kIsInitialized == false) {
			tuple.k = tupleStr.size();
			kIsInitialized = true;
		}
		tuple.InitializeFromString((unsigned char*) tupleStr.c_str(), tuple.k);
		queryCount.push_back(0);
		queryNames.push_back(seq);
		queryTuples.push_back(tuple);
		queryPositions.push_back(pos);
		tupleToIndex[tuple].push_back(index);
		tuple.SetReverseComplement(tuplerc);
		tupleToIndex[tuplerc].push_back(index);
		++index;
		
	}

	map<BitNucVector, int>::iterator searchResult;
	int readNumber = 0;
	
	Counts counts;
	counts.in = &fastqFile;
	counts.pairIn = &queryMPFile;
	counts.isFastq = isFastq;
	counts.tupleToIndex = &tupleToIndex;
	counts.queryCount = &queryCount;
	counts.readIndex = &readNumber;
	counts.noDoubleCount = noDoubleCount;
	counts.readsOut = &readsOut;
	counts.maxNReads = maxNReads;
	const int idLen=10;
	char id[idLen+1];
	id[idLen] = '\0';
	srand (time(NULL));
	rand_str(id, idLen);

	string readerName = string("/semreader_") + string(id);
	string countName  = string("/semcount_") + string(id);
	string semfastqwriterName = string("/semfastqwriter_") + string(id);

	semreader     = sem_open(readerName.c_str(), O_CREAT, 0644, 1);
	if (semreader == NULL) {
		cout << "ERROR opening semaphore. ERRNO " << errno << " " << readerName << endl;
		exit(1);
	}
	semcount      = sem_open(countName.c_str(), O_CREAT, 0644, 1);
	if (semreader == NULL) {
		cout << "ERROR opening semaphore. ERRNO " << errno << " " << countName << endl;
		exit(1);
	}
	semfastqwriter = sem_open(semfastqwriterName.c_str(), O_CREAT, 0644, 1);
	if (semfastqwriter == NULL) {
		cout << "ERROR opening semaphore. ERRNO " << errno << " " << semfastqwriterName << endl;
		exit(1);
	}

	sem_init(semreader, 1, 1);
	sem_init(semcount, 1, 1);
	sem_init(semcount, 1, 1);

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
		string tupleString = queryTuples[i].ToString();
		outFile << queryNames[i] << "\t" << tupleString << "\t" << queryCount[i] << "\t" << queryPositions[i] << endl;
	}

	if (matchedReadsFileName != "") {
		readsOut.close();
	}

	outFile.close();
}


