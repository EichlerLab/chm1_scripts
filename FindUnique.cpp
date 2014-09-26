#include "FASTASequence.h"
#include "FASTAReader.h"
#include "algorithms/anchoring/MapBySuffixArray.h"
#include "datastructures/suffixarray/SuffixArrayTypes.h"
#include <set>

typedef DNASuffixArray T_SuffixArray;
int main(int argc, char* argv[1]) {
	if (argc < 3) {
		cout << "Usage: findUnique genome.fasta query.fasta effective_k [options]" << endl;
		cout << "  genome.fasta.sa must exist." << endl;
		cout << "  Finds sequences at least effective_k in length that are unique." << endl;
		cout << "  -max m       Allow up to m matches" << endl;
		cout << "  -minLength l Ensure the length of the match is at least this." << endl;
		cout << "  -prefix p n  Allow up to n matches across a prefix of length p" << endl;
		cout << "  -suffix s n  Allow up to n matches across a suffix of length s" << endl;
		cout << "               Prefix and suffix options override max." << endl;
		cout << "  -out file    Print queries to this output file (query.fasta.queries)" << endl;
		exit(0);
	}

	DNASuffixArray sarray;
	
	string genomeFileName = argv[1];
	string suffixArrayFileName = genomeFileName + ".sa";
	
	FASTAReader reader;
	FASTASequence genome;

	int maxN = 0;

	int prefix = 0;
	int suffix = 0;
	int prefixN = 0;
	int suffixN = 0;
	int argi = 4;
	string outputFileName = "";
	int minLength = 0;
	while (argi < argc) {
		if (strcmp(argv[argi], "-max") == 0) {
			++argi;
			maxN = atoi(argv[argi]);
		}
		else if (strcmp(argv[argi], "-prefix") == 0) {
			++argi;
			prefix = atoi(argv[argi]);
			++argi;
			prefixN = atoi(argv[argi]);
		}
		else if (strcmp(argv[argi], "-suffix") == 0) {
			++argi;
			suffix = atoi(argv[argi]);
			++argi;
			suffixN = atoi(argv[argi]);
		}
		else if (strcmp(argv[argi], "-out") == 0) {
			++argi;
			outputFileName = argv[argi];
		}
		else if (strcmp(argv[argi], "-minLength") == 0) {
			++argi;
			minLength = atoi(argv[argi]);
		}
		++argi;
	}

	reader.Initialize(genomeFileName);
	reader.ReadAllSequencesIntoOne(genome);
	sarray.Read(suffixArrayFileName);

	FASTAReader queryReader;
	FASTASequence querySequence;
	string queryFileName = argv[2];
	int maxLength = atoi(argv[3]);
	string summaryTableFileName = queryFileName + ".summary";
	if (outputFileName == "") {
		outputFileName = queryFileName + ".queries";
	}
		
	
	ofstream summaryTable(summaryTableFileName.c_str());
	ofstream outputFile(outputFileName.c_str());

	queryReader.Initialize(queryFileName);

	while (queryReader.GetNext(querySequence)) {
		int i;
		cerr << "searching " << querySequence.title << endl;
		if (querySequence.length < maxLength) {
			continue;
		}

		int nMatches = 0;
		querySequence.ToUpper();
		int localMax;
		for (i = 0; i < querySequence.length - maxLength + 1; i++) {
			if ((i + 1) % 100000 == 0) {
				cerr << "processed: " << i + 1 << endl;
			}

			int lcpLength;
			vector<SAIndex> lcpLeftBounds, lcpRightBounds;
			vector<SAIndex> rclcpLeftBounds, rclcpRightBounds;
			localMax = maxN;
			if (i < prefix) {
				localMax = prefixN;
			}
			if (i >= querySequence.length - suffix) {
				localMax = suffixN;
			}
			if (querySequence.length - i <= maxLength) {
				continue;
			}
			if (querySequence.seq[i] == 'N') {
				continue;
			}
			lcpLength = sarray.StoreLCPBounds(genome.seq, genome.length, // The string which the suffix array is built on.
																				&querySequence.seq[i], querySequence.length-i,
																				true,
																				maxLength,
																				lcpLeftBounds, lcpRightBounds,
																				false);
			if (lcpLength < minLength) {
				continue;
			}
			if (lcpLength < maxLength or 
					lcpRightBounds.size() == 0 or 
					(lcpRightBounds.size() > 0 and 
					 lcpLeftBounds.size() > 0 and  
					 lcpRightBounds[lcpRightBounds.size() - 1] - lcpLeftBounds[lcpLeftBounds.size()-1] <= localMax)) {

				FASTASequence rc;
				DNASequence subseq;
				subseq.ReferenceSubstring(querySequence, i, maxLength);
				subseq.MakeRC(rc);
				int rclcpLength;
				int numForwardMatches;
				if (lcpLength == 0) {
					numForwardMatches = 0;
				}
				else {
					numForwardMatches = lcpRightBounds[lcpRightBounds.size() - 1] - lcpLeftBounds[lcpLeftBounds.size()-1];
				}
				rclcpLength = sarray.StoreLCPBounds(genome.seq, genome.length, // The string which the suffix array is built on.
																						rc.seq, maxLength,
																						true,
																						rclcpLength,
																						rclcpLeftBounds, rclcpRightBounds,
																						false);

				string rcstr((const char*)rc.seq, rc.length);

				if (rclcpLength < maxLength or 
						rclcpRightBounds.size() == 0 or
						(numForwardMatches + 
						 rclcpRightBounds[rclcpRightBounds.size() - 1] -
						 rclcpLeftBounds[rclcpLeftBounds.size()-1] <= localMax)) 
					{
						char* substr = new char[maxLength+1];
						substr[maxLength] = '\0';
						memcpy(substr, &querySequence.seq[i], maxLength);

						//						string substr = string((const char*) querySequence.seq, i, maxLength);
						
						outputFile << querySequence.title << "\t" << substr << "\t" << i << endl;

						++nMatches;
						delete[] substr;
						//					}
					}
				rc.Free();
			}

		}
		summaryTable << querySequence.title << "\t" << nMatches << endl;
		querySequence.Free();
	}
	outputFile.close();
	genome.Free();
}


