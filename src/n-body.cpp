#include <mpi.h>
#include <string>

using namespace std;

const MPI_Comm comm = MPI_COMM_WORLD;
const int root = 0;

typedef struct s {
    double xPos;
    double yPos;
    double mass;
    double force;
} Body;

int main (int argc, char **argv) {
    int rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(comm, &rank);

    if (rank == root) {
        int n;
        cin >> n;

        Body bodies[n];

        for (int i = 0; i < n; i++) {
            Body body;

            cin >> body.xPos;
            cin >> body.yPos;
            cin >> body.mass;
            cin >> body.force;

            bodies[i] = body;
        }
    } else {

    }

    MPI_Finalize();
    return 0;
}
