#include <mpi.h>
#include <string>

using namespace std;

const MPI_Comm comm = MPI_COMM_WORLD;
const int root = 0;

int main (int argc, char **argv) {
    int rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(comm, &rank);

    if (rank == root) {

    } else {

    }

    MPI_Finalize();
    return 0;
}
