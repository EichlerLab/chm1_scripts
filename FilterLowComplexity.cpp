#include <string>
#include <map>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <iomanip>      // std::setw

#include "BitNucVector.h"
#include "FASTAReader.h"
#include "FASTASequence.h"
#include <math.h>

using namespace std;
class WordCount {
public:
	int count;
	Word w;
	int operator<(const WordCount &rhs) const {
		return count > rhs.count;
	}
};

void ResetWordQueries(vector<map<Word, int> > &queries) {
	int i;
	for (i = 0; i < queries.size(); i++) {
		map<Word, int>::iterator mapIt, mapEnd;
		mapEnd = queries[i].end();
		mapIt = queries[i].begin();
		while (mapIt != mapEnd) {
			(*mapIt).second = 0;
			++mapIt;
		}
	}
}


int main(int argc, char* argv[]) {
	if (argc < 4) {
		cout << "Usage: filterlc fasta k tableOut [options]" << endl;
		cout << "Options:" << endl;
		cout << "  -c \tPrint word counts" << endl;
		cout << "  -f FLOAT\tminimal enrichment [0,1]" << endl;
		cout << "  -n INT\tminimum number of enriched [5]" << endl;
		cout << "  -w fasta \tPrint sequences enriched for words in this fasta file." << endl;
		cout << "  -o fasta \tPrint matching sequences to this fasta file." << endl;
		exit(1);
	}

	string fastaFileName = argv[1];
	int k = atoi(argv[2]);	
	string outTableFileName = argv[3];

	int argi = 4;
	float fraction = 0.05;
	int minEnriched = 5;
	bool printWordCounts = false;
	bool printWordRates = false;
	bool printMeanRate = false;
	string wordQueryFile = "";
	string outFastaFileName = "";
	while (argi < argc) {
		if (strcmp(argv[argi], "-f") == 0) {
			++argi;
			fraction = atof(argv[argi]);
		}
		if (strcmp(argv[argi], "-n") == 0) {
			++argi;
			minEnriched = atoi(argv[argi]);
		}
		if (strcmp(argv[argi], "-c") == 0) {
			printWordCounts = true;
		}
		if (strcmp(argv[argi], "-r") == 0) {
			printWordRates = true;
		}
		if (strcmp(argv[argi], "-m") == 0) {
			printMeanRate = true;
		}
		if (strcmp(argv[argi], "-w") == 0) {
			++argi;
			wordQueryFile = argv[argi];
		}
		if (strcmp(argv[argi], "-o") == 0) {
			++argi;
			outFastaFileName = argv[argi];
		}
		++argi;
	}
	
	ifstream in;
	ofstream out, outFasta;
	
	in.open(fastaFileName.c_str());
	out.open(outTableFileName.c_str());
	BitNucVector::k = k;
	
	int numWords = 1 << (2*k);

	if (outFastaFileName != "") {
		outFasta.open(outFastaFileName.c_str());
	}

	vector <map<Word, int > > wordQueries;
	vector <string> wordQueriesTitles;
	if (wordQueryFile != "") {
		FASTAReader reader;
		reader.Initialize(wordQueryFile);
		FASTASequence seq;
		int seqIndex = 0;
		while (reader.GetNext(seq)) {
			int i;
			wordQueries.push_back(map<Word, int>());
			for (i = 0; i < seq.length - k + 1; i++) {
				BitNucVector v, key;
				if (v.InitializeFromString((unsigned char*) &seq.seq[i], k)) {
					v.SetKey(key);
					wordQueries[seqIndex][key.data] = 0;
				}
			}				
			wordQueriesTitles.push_back(seq.GetName());
			++seqIndex;
		}
		seq.Free();
	}


	FASTAReader seqReader;
	int readIndex = 0;
 	while (in) {
		++readIndex;
		if (readIndex % 10000 == 0) {
			cerr << "filerlc " << readIndex << endl;
		}
		string title;
		getline(in, title);
		string seq;
		while (in.peek() != '>' and in) {
			char p = in.peek();
			string line;
			(in) >> line;
			if (line == "") {
				break;
			}
			assert(line[0] != '>');
			seq += line;
			if (in.peek() == '\n') {
				in.get();
			}
		}
		BitNucVector tuple;
		tuple.k = k;

		int i;
		map<Word, int> wordCount;
		// 
		// Count all words in this file.
		//
		if (seq.size() < k) {
			continue;
		}
		for (i = 0; i < seq.size() - k + 1; i++) {
			if (tuple.InitializeFromString((unsigned char*) &seq.c_str()[i], tuple.k)) {				
				BitNucVector key;
				tuple.SetKey(key);
				if (wordCount.find(key.data) == wordCount.end()) {
					wordCount[key.data] = 1;
				}
				else {
					wordCount[key.data]++;
				}
			}
			else {
				i++;
			}
		}

		float pWord = 1.0/numWords;
		float logPWord = -log(pWord);
		
		int maxNumEnriched = 0;
		map<Word, int>::iterator mapIt, mapEnd;
		if (wordQueryFile != "") {
			//
			// First compute the expected number of matches.
			//
			int s;
			float maxEnrichment = 0;
			int maxEnrichmentI = 0;
			int maxNQueryWords = 0;
			int maxCountQueryWords = 0;
			int numEnriched = 0;

			for (s = 0; s < wordQueries.size(); s++ ) {
				int nQueryWords = wordQueries[s].size();
				int expNQueryWords = nQueryWords * pWord * seq.size();
				int countQueryWords = 0;
				float totalEnrichment = 0;
				mapIt = wordQueries[s].begin();
				mapEnd = wordQueries[s].end();
				
				while (mapIt != mapEnd) {
					if (  wordCount.find((*mapIt).first) != wordCount.end()) {
						countQueryWords += wordCount[mapIt->first];
						totalEnrichment += (wordCount[mapIt->first]) / (pWord * seq.size());
					}
					++mapIt;
				}
				//				totalEnrichment = countQueryWords / float(wordQueries[s].size());
				//				totalEnrichment /= seq.size();
				totalEnrichment = countQueryWords / float(seq.size());
				if (totalEnrichment  > maxEnrichment) {
					maxEnrichmentI = s;
					maxEnrichment = totalEnrichment;
					maxNQueryWords = expNQueryWords;
					maxCountQueryWords = countQueryWords;
				}
			}
			if (maxEnrichment >= fraction) {
				out << title 
						<< "\t" 
						<< wordQueriesTitles[maxEnrichmentI] 
						<< "\t" << seq.size() 
						<< "\t" << maxNQueryWords 
						<< "\t" << maxCountQueryWords 
						<< "\t" << maxEnrichment 
						<< "\t" << numEnriched << endl;
				if (outFastaFileName != "") {
					outFasta << title << endl;
					outFasta << seq << endl;
				}
			
			}

			continue;
		}
		
		
		vector<WordCount> counts;
		
		mapEnd = wordCount.end();
		counts.resize(wordCount.size());
		i = 0;
		for (mapIt = wordCount.begin(); mapIt != wordCount.end(); mapIt++) {
			counts[i].w = (*mapIt).first;
			counts[i].count = (*mapIt).second;
			i++;
		}

		sort(counts.begin(), counts.end());

		int nEnriched = 0;
		for (i = 0; i < counts.size(); i++) {
			if (counts[i].count > fraction * seq.size()) {
				nEnriched ++;
			}
		}

		if (printWordCounts) {
			int i;
			out << title << " ";
			vector<float> freqs(counts.size());
			for (i = 0; i < counts.size(); i++) {
				freqs[i] = counts[i].count / ((float)seq.size());
			}
			std::sort(freqs.begin(), freqs.end());
			std::reverse(freqs.begin(), freqs.end());
			float total = 0;
			for (i = 0; i < counts.size() and i < 10; i++) {
				//out << std::setprecision(2) << freqs[i] << " " ;
				total += freqs[i];
			}
			out << std::setprecision(2) << total << endl;
			//			out << endl;
		}
		else if (printWordRates) {
			int i;
			out << title << " ";
			for (i = 0; i < counts.size(); i++) {
				out << counts[i].count / float(counts.size()) << " ";
			}
			out << endl;

		}
		else if (printMeanRate){ 
			float sum = 0;
			for (i = 0; i < minEnriched; i++) {
				sum += counts[i].count;
			}
			out << title << "\t" << sum / (counts.size()*minEnriched) << endl;
		}
		else {
			out << nEnriched << "\t" << counts.size() << "\t" << float(nEnriched)/counts.size() << endl;
		}
		
		/*			std::setw(4);
			out << title << "\t" << nEnriched << "\t" << counts.size() << "\t" << float(nEnriched)/counts.size() << "\t";
			for (i = 0; i < 5 and i  < counts.size(); i++) {
				tuple.data = counts[i].w;
				//				std::setwidth(4);
				out << "\t" << tuple.ToString() << "\t" << counts[i].count / float(seq.size());
			}
			out << "\t" << seq.size() << "\t" << seq;
			out << endl;
			//		}
			*/
	}
	if (outFastaFileName != "") {
		outFasta.close();
	}
}
