#include <iostream>
#include <fstream>
#include <map>
#include <string>
#include <stdlib.h>
#include "FASTAReader.h"
#include "FASTASequence.h"
#include <queue>

using namespace std;
void StoreBins(vector<int> &counts, int seqLength, int kmer, int binSize, vector<int> &bins) {
	int L = seqLength - kmer + 1;
	bins.resize(seqLength / binSize+ (seqLength % binSize != 0 ? 1 : 0), 0);
	fill(bins.begin(), bins.end(), 0);
	int b = 0;
	int i;
	for (i = 0; i < L; i++) {
		assert(b < bins.size());
		bins[b] += counts[i];
		if (i > 0 and i % binSize == 0) {
			b++;
		}
	}
}

float StoreGC(DNASequence &seq, int start, int end) {
	int ngc = 0;
	int i;
	for (i = start; i < end; i++) {
		if (seq.seq[i] == 'G' or seq.seq[i] == 'C') { ngc++; }
	}
	return ((float)ngc) / (end - start);
}

float StoreIdentity(vector<int> &bins, int start, int end) {
	float sum =0;
	int i; 
	for (i = start; i < end; i++) { sum += bins[i]; }
	float m = sum/(end - start);
	float sumIdentity = 0;
	for (i = start; i < end; i++) { float b=bins[i]; sumIdentity +=  min(b,m)/max(b,m);}
	return sumIdentity/(end-start);
}

int main(int argc, char* argv[]) {
	if (argc < 4) {
		cout << "Usage: detectPeriod file.fasta k binfile [-v]" << endl;
		exit(0);
	}
	string inFile = argv[1];
	int kmer = atoi(argv[2]);
	
	ofstream binFile(argv[3]);
	int argi = 4;
	bool verbose = false;
	bool printHeader = false;
	while (argi < argc) {
		if (strcmp(argv[argi], "-v") == 0) {
			verbose = true;
		}
		if (strcmp(argv[argi], "-header") == 0) {
			printHeader = true;
		}
		++argi;
	}
	FASTAReader reader;
	reader.Initialize(inFile);
	FASTASequence seq;
	if (printHeader) {
		cout << "title\tlength\tmaxRunLength\tmaxRunStart\tmaxRunEnd\tbinSize" << endl;	
	}
	int seqIndex = 0;

	while (reader.GetNext(seq)) {
		seq.ToUpper();

		string kSeq;
		int i;
		map<string, int> wordCount;
		int L = seq.length - kmer + 1;
		for (i = 0; i < L; i++) {
			kSeq = string((char*)&seq.seq[i], kmer);
			if (wordCount.find(kSeq) == wordCount.end()) {
				wordCount[kSeq] = 1;
			}
			else {
				wordCount[kSeq]++;
			}
		}

		vector<int> counts(L, 0);
		for (i = 0; i < L; i++) {
			kSeq = string((char*)&seq.seq[i], kmer);
			assert(wordCount.find(kSeq) != wordCount.end());
			counts[i] = wordCount[kSeq];
		}
	
		int binSize; 
		int maxBinSize = min((int)seq.length, 30*kmer);
		int globalMaxRunLength= 0;
		int globalMaxRunStart = 0;
		int globalMaxRunEnd   = 1;
		int globalMaxRunBinSize = binSize;
		float globalMaxPctIdentity = 0;
		vector<int> bins;
		for (binSize = kmer; binSize < maxBinSize; binSize++) {

			StoreBins(counts, seq.length, kmer, binSize, bins);
			int maxRunLength = 1;
			int j;
			int runStart = 0;
			int runEnd = 0;
			float pctIdentity = 0;
			for (i = 0; i < bins.size() - 1; i++) {
				bool cont = true;
				if (verbose) {
					cout << binSize << " " << i * binSize << " " << bins[i] << endl;
				}
				for (j = i+1; j < bins.size() and cont; j++) {
					if (!(bins[i]*0.95 < bins[j] and bins[i] *1.053 > bins[j])) {
						
						cont = false;
					}
				}
				int runLength = j - i;
				if (runLength > maxRunLength) {
					runStart = i;
					runEnd = j;
					maxRunLength = runLength;
					pctIdentity = StoreIdentity(bins, i, j);
				}
			}

			if (globalMaxRunLength < maxRunLength) {
				globalMaxRunLength = maxRunLength;
				globalMaxRunStart  = runStart;
				globalMaxRunEnd    = runEnd;
				globalMaxRunBinSize= binSize;
				globalMaxPctIdentity = pctIdentity;
			}
		}
		float gc = StoreGC(seq, globalMaxRunStart*globalMaxRunBinSize , globalMaxRunEnd*globalMaxRunBinSize);
		string repeat((const char*) seq.seq, globalMaxRunStart*globalMaxRunBinSize, (globalMaxRunStart+1)*globalMaxRunBinSize);
		cout << seq.title << "\t" << seq.length <<"\t" << globalMaxRunLength  << "\t" << globalMaxRunStart*globalMaxRunBinSize << "\t" << globalMaxRunEnd*globalMaxRunBinSize << "\t" << globalMaxRunBinSize << "\t" << globalMaxPctIdentity << "\t" << gc << "\t" << repeat << endl;

		StoreBins(counts, seq.length, kmer, globalMaxRunBinSize, bins);
		for (i = 0; i < bins.size(); i++) {
			int label= 0;
			if (i >= globalMaxRunStart and i <= globalMaxRunEnd) {
				label = 1;
			}
			binFile << seq.title << "\t" << label << "\t" << bins[i] << endl;
		}
		seq.Free();
		++seqIndex;
	}
	binFile.close();
}
