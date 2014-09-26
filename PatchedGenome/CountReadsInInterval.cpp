#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <map>
#include <string>
#include <time.h>

using namespace std;

int main(int argc, char* argv[]) {

	ifstream countFile(argv[1]);



	int total = 0;
	int line = 0;
	vector<string> countChrom;
	vector<int>    countPos;
	vector<int>    countValues;
	string chrom;
	int pos;
	int count;
	const int bufSize = 32768;
	char fileBuf[bufSize];
	countFile.rdbuf()->pubsetbuf(fileBuf, bufSize);
		
	while (countFile) {
		line +=1;
		if ((countFile >> chrom >> pos >> pos >> count) == 0) {
			break;
		}
		else {
			countChrom.push_back(chrom);
			countPos.push_back(pos);
			countValues.push_back(count);
		}
	}
	cerr << "done reading queries." << endl;
	int i;
	for (i = 2; i < argc; i++) {
		ifstream intvFile(argv[i]);
		map<string, vector<int> > starts, ends, counts;

		while (intvFile) {
			string chrom;
			int start;
			int end;
			if ( (intvFile >> chrom >> start >> end) == 0) {
				break;
			}
			else {
				string line;
				getline(intvFile, line);
			}
			if (starts.find(chrom) == starts.end()) {
				starts[chrom] = vector<int>();
				ends[chrom] = vector<int>();
				counts[chrom] = vector<int>();
			}
			starts[chrom].push_back(start);
			ends[chrom].push_back(end);
			counts[chrom].push_back(0);
		}

		int q;
		total = 0;
		time_t startTime;
		time(&startTime);
		for (q = 0; q < countChrom.size(); q++) {
			chrom = countChrom[q];
			pos   = countPos[q];
			count = countValues[q];
			if (starts.find(chrom) == starts.end()) {
				continue;
			}
			vector<int>::iterator startIndex;
			startIndex = lower_bound(starts[chrom].begin(), starts[chrom].end(), pos);
			if (startIndex != starts[chrom].end()) {
				int index = startIndex - starts[chrom].begin();
				
				if (pos == starts[chrom][index]) {
					total += count;
				}
				else {
					index -= 1;
					if ( index >= 0 and starts[chrom][index] <= pos and ends[chrom][index] > pos) {
						total += count;
					}
				}
			}
		}
		time_t endTime;
		time(&endTime);
		cerr << difftime(endTime, startTime) << " s" << endl;
		cout << argv[i] << "\t" << total << endl;
	}
}
	
		

