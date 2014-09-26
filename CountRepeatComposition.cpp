#include "FASTASequence.h"
#include "FASTAReader.h"
#include "algorithms/anchoring/MapBySuffixArray.h"
#include "datastructures/suffixarray/SuffixArrayTypes.h"
#include <set>

typedef DNASuffixArray T_SuffixArray;
int minCount = 0;
int minLength = 0;
int maxLength = 99999999;
int debug = 0;
int FindMaxPattern(FASTASequence &genome, DNASuffixArray &sarray,
									 int start, int end,
									 vector<string> &patterns, vector<int> &counts, int lcp=0) {
	int i = start;
	int curStart = start;
	if (debug) {
		cerr << start << " " << end << " " << lcp << endl;
	}

	while (i < end) {
		while (i < end and 
					 sarray[i] + lcp < genome.length and 
					 genome.seq[sarray[i]+lcp] == genome.seq[sarray[curStart] + lcp]) {
			i++;
		}
		if (i > curStart) {
			//
			// Found a matching pattern.
			//
			string pat((const char*) &genome.seq[sarray[curStart]], lcp + 1);
			int count = i - curStart;
			if (count >= minCount and pat.size() >= minLength and pat.size() <= maxLength) {
				patterns.push_back(pat);
				counts.push_back(i - curStart);
			}

			//
			// Recurse in this range to find longer patterns
			//
			if (i - curStart > minCount and 
					lcp + 1 <= maxLength) {
				FindMaxPattern(genome, sarray, curStart, i, 
											 patterns, counts, lcp+1);
			}
		}
		else {
			i++;
		}
		curStart = i;
	}
}

int main(int argc, char* argv[1]) {
	if (argc < 2) {
		cout << "Usage: countRepeatComposition genome.fasta [--minCount c] [--minLength l] [--maxLength l]" << endl;
		cout << "  Counts the number of times substrings appear in a sequence. For example, given" << endl
				 << "  the sequence TTATTATTGTTATTATTG , --minLength 2 --minCount 2 --maxLength 5 " << endl
				 << "4		2	AT\n"
"4		3	ATT\n"
"2		4	ATTA\n"
"2		4	ATTG\n"
"4		2	TA\n"
"4		3	TAT\n"
"4		4	TATT\n"
"2		2	TG\n"
"6		2	TT\n"
"4		3	TTA\n"
"4		4	TTAT\n"
"2		3	TTG" << endl;

		exit(0);
	}

	DNASuffixArray sarray;
	
	string genomeFileName = argv[1];

	int argi = 2;
	while (argi < argc) {
		if (strcmp(argv[argi], "--minCount") == 0) {
			minCount = atoi(argv[++argi]);
		}
		else if (strcmp(argv[argi], "--minLength") == 0) {
			minLength = atoi(argv[++argi]);
		}
		else if (strcmp(argv[argi], "--maxLength") == 0) {
			maxLength = atoi(argv[++argi]);
		}
		else if (strcmp(argv[argi], "--debug") == 0) {
			debug = 1;
		}

		++argi;
	}
	
	FASTAReader reader;
	FASTASequence genome;

	reader.Initialize(genomeFileName);
	reader.GetNext(genome);
	
	//
	// Build the suffix array for fast matching.
	//
	FASTASequence genomeCopy;
	genomeCopy.Copy(genome);
	genomeCopy.ToThreeBit();		
	vector<int> alphabet;
	sarray.InitThreeBitDNAAlphabet(alphabet);
	sarray.LarssonBuildSuffixArray(genomeCopy.seq, genome.length, alphabet);

	vector<string> patterns;
	vector<int> counts;
	genome.ToUpper();
	FindMaxPattern(genome, sarray, 0, genome.length, patterns, counts);

	int i;
	for (i = 0; i < patterns.size(); i++) {
		cout << counts[i] << "\t" << "\t" << patterns[i].size() << "\t" << patterns[i] << endl;
	}

	return 0;

}


