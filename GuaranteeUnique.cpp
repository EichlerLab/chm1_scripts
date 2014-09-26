#include "FASTAReader.h"
#include "FASTASequence.h"
#include <semaphore.h>
#include <errno.h>

sem_t *semreader;

class Data {
public:
	FASTASequence *genome;
	ifstream *queryFile;
};
	

void SearchGenome(void *data) {

	ifstream &queryFile = *((Data*)data)->queryFile;
	FASTASequence &genome = *((Data*)data)->genome;

	while (true) {

		int retval;

		retval = sem_wait(semreader);		

		string name, query;

		int pos;
		if (queryFile) {

			if ((queryFile >> name >> query >> pos) == 0) {

				sem_post(semreader);
				return ;
			}

		}
		else {

			sem_post(semreader);
			return ;
		}
		sem_post(semreader);

		long i;
		for (i =  0; i < query.size(); i++) { query[i] = toupper(query[i]);}

		string rc = query;

		for (i =  0; i < query.size(); i++) { 
			switch(query[i]) {
			case 'A': 
				rc[query.size() - i - 1] = 'T';
				break;
			case 'T':
				rc[query.size() - i - 1] = 'A';
				break;
			case 'C':
				rc[query.size() - i - 1] = 'G';
				break;
			case 'G':
				rc[query.size() - i - 1] = 'C';
				break;
			}
		}
		
		const char *queryPtr = query.c_str();
		const char *rcPtr = rc.c_str();
		int queryLen = query.size();
		int middle = query.size() / 2;
		int count = 0;
		
		for (i = 0; i < genome.length- queryLen; i++) {
			if (genome.seq[i+0] == queryPtr[0] &&
					genome.seq[i+middle] == queryPtr[middle] &&
					genome.seq[i+queryLen-1] == queryPtr[queryLen-1]) {
				if (strncmp((const char*) &genome.seq[i], queryPtr, queryLen) == 0) {
					count++;
				}
			}
			if (genome.seq[i+0] == rcPtr[0] &&
					genome.seq[i+middle] == rcPtr[middle] && 
					genome.seq[i+queryLen-1] == rcPtr[queryLen-1]) {
				if (strncmp((const char*) &genome.seq[i], rcPtr, queryLen) == 0) {
					count++;
				}
			}				
		}
		sem_wait(semreader);
		cout << count << "\t" << name << "\t" << pos << "\t" << query << endl;
		sem_post(semreader);
	}
	

}
int main(int argc, char* argv[]) {
	semreader     = sem_open("/semreader",     O_CREAT, 0644, 1);
	sem_init(semreader, 1, 1);
	string genomeFileName = argv[1];
	string queryFileName = argv[2];
	FASTAReader reader;
	FASTASequence genome;
	
	reader.Initialize(genomeFileName);
	reader.ReadAllSequencesIntoOne(genome);
	genome.ToUpper();

	ifstream queryFile(queryFileName.c_str());
	int nproc =12;
	pthread_attr_t *threadAttr = new pthread_attr_t[nproc];
	Data data;
	data.genome = &genome;
	data.queryFile = &queryFile;
	int t;	
	for (t = 0; t < nproc; t++ ){
		pthread_attr_init(&threadAttr[t]);
	}
	pthread_t *threads = new pthread_t[nproc];


	for (t = 0; t < nproc; t++) {
		pthread_create(&threads[t], &threadAttr[t], (void* (*)(void*))SearchGenome, &data);
	}

	for (t = 0; t < nproc; t++) {
		pthread_join(threads[t], NULL);
	}


	
}


		
