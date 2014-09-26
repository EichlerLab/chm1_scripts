#include "FASTAReader.h"
#include "FASTASequence.h"
#include "utils.h"

#include <map>

int main(int argc, char* argv[]) {
	if (argc != 3) {
		cout << "Usage: pmme in.fasta out.fasta" << endl;
		exit(0);
	}
	string inFasta = argv[1];
	string outFasta = argv[2];
	map<string, int> nMaskedMap;
	FASTAReader reader;
	ofstream outFile;
	CrucialOpen(outFasta, outFile, std::ios::out);

	vector<FASTASequence> sequences;
	reader.Initialize(inFasta);
	reader.ReadAllSequences(sequences);
	cerr << "done reading sequences." << endl;
	int s;
	vector<int> nMaskedVect(sequences.size(), 0);
	for (s = 0; s < sequences.size(); s++) {
		int nMasked = 0;
		int i;
		for (i = 0; i < sequences[s].length; i++) {
			if (sequences[s].seq[i] == 'N') {
				nMasked++;
			}
		}
		nMaskedVect[s] = nMasked;
		string title = sequences[s].title;
		if (nMaskedMap.find(title) == nMaskedMap.end()) {
			nMaskedMap[title] = nMasked;
		}
		else {
			if (nMasked > nMaskedMap[title]) {
				nMaskedMap[title] = nMasked;
			}
		}
	}

	for (s = 0; s < sequences.size(); s++) {
		string title = sequences[s].title;
		if (nMaskedMap[title] == nMaskedVect[s]) {
			sequences[s].PrintSeq(outFile);
			nMaskedMap[title] = 1000000;
		}
	}
	outFile.close();


}
