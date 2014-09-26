#include <iostream>
#include <fstream>
#include <map>
#include <vector>


int main(int argc, char* argv[]) {

	ifstream regionIn(argv[1]);
	ifstream faiFile(argv[2]);
	int read = atoi(argv[3]);
	int bin  = atoi(argv[4]);

	map<string, int> chromLen;

	while (faiFile) {
		string chrom;
		long long len, total, lineLen, lineBytes;
		if (!(fai >> chrom >> len >> total >> lineLen >> lineBytes)){ 
			break;
		}
		else {
			chromLen[chrom] = len;
		}
	}


}
