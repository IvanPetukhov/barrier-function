#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <ctime>
#include <cassert>
#include <math.h>
#include <windows.h>
#include <mpi.h>

using namespace::std;

int myRank;
int numProcs;
char host[256];

int degree2(int size)
{
	int degree = 0;
	int twoPowered = 1;
	while(twoPowered <= size){
		twoPowered *= 2;
		degree++;
	}
	degree--;
	return degree;
}

int power2(int n){
	int res = 1;
	for(int i = 0; i < n; i++){
		res *= 2;
	}
	return res;
}

void barrier(int rank, int size)
{
	int buf[2] = {0, 1};
	MPI_Status status;
	int readRank = 1;
	int maxN = degree2(size) + 1;
	if (rank == 0) {
		for (int i = 0; i < maxN; ++i){
			MPI_Recv(buf, 2, MPI_INT, readRank, 0, MPI_COMM_WORLD, &status);
			readRank *= 2;
		}
		readRank /= 2;
	}
	for (int i = 1; i <= maxN; ++i){
		if(rank % power2(i) == power2(i - 1)){
			for (int j = 0; j < i - 1; ++j){
				if (rank + power2(j) < size){
    				MPI_Recv(buf, 2, MPI_INT, rank + power2(j), 0, MPI_COMM_WORLD, &status);
				}
			}
			if (rank % power2(degree2(rank)) == 0){
				MPI_Send(buf, 2, MPI_INT, 0, 0, MPI_COMM_WORLD);
			} else {
				MPI_Send(buf, 2, MPI_INT, rank - power2(i - 1), 0, MPI_COMM_WORLD);
			}
			break;
		}
	}
	readRank = 1;
	if (rank == 0) {
		for (int i = 0; i < maxN; ++i){
			MPI_Send(buf, 2, MPI_INT, readRank, 0, MPI_COMM_WORLD);
			readRank *= 2;
		}
		readRank /= 2;
	}
	for (int i = 1; i <= maxN; ++i){
		if(rank % power2(i) == power2(i - 1)){
			if (rank % power2(degree2(rank)) == 0){
				MPI_Recv(buf, 2, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
			} else {
				MPI_Recv(buf, 2, MPI_INT, rank - power2(i - 1), 0, MPI_COMM_WORLD, &status);
			}
			for (int j = 0; j < i - 1; ++j){
				if (rank + power2(j) < size){
					MPI_Send(buf, 2, MPI_INT, rank + power2(j), 0, MPI_COMM_WORLD);
				}
			}
			break;
		}
	}
}

int main(int argc, char* argv[])
{
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
  gethostname(host, sizeof(host)/sizeof(host[0]));
  cout << "Process " << myRank << " of " << numProcs << " is running on '" << host << "'." << endl;
  barrier(myRank, numProcs);
  MPI_Finalize();
  return 0;
}