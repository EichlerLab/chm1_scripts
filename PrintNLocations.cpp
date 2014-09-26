#include "FASTAReader.h"
#include "FASTASequence.h"

using namespace std;

int main(int argc, char* argv[]) {

	FASTAReader reader;
	if (argc < 2) {
		cout << "Usage nloc file.fasta" << endl;
		exit(1);
	}
	string name = argv[1];
	reader.Initialize(name);

	FASTASequence seq;
	while (reader.GetNext(seq)) {
		int i, j;
		i = 0;
		while (i < seq.length) {
			j = i;
			while (seq.seq[j] == 'N') {
				j++;
			}
			if (j > i) {
				cout << seq.title << "\t" << i << "\t" << j << "\t" << j - i << endl;
			}
			i = j+1;
		}
		seq.Free();
	}
}
