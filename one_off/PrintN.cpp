#include <FASTAReader.h>
#include <string>

using namespace std;

int main(int argc, char* argv[]) {
	FASTAReader reader;
	string fn = argv[1];
	reader.Initialize(fn);

	FASTASequence seq;
	while (reader.GetNext(seq)) {
		int i, j;
		i = 0;
		while (i < seq.length) {
			j = i;
			while (j < seq.length and seq.seq[j] == 'N') {
				j++;
			}
			if (j - i > 100000) {
				cout << seq.title << "\t" << i << "\t" << j << "\t" << j -i  << endl;
			}
			i=j+1;
		}
	}
}
