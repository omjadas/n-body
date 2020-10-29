#include <math.h>
#include <mpi.h>
#include <string>

using namespace std;

const MPI_Comm comm = MPI_COMM_WORLD;
const int root = 0;
const double G = 6.67408 * pow(10, -11);

typedef struct s {
    double xPos;
    double yPos;
    double mass;
    double force;
} Body;

double distance(Body b1, Body b2);
double force(Body b1, Body b2);

int main(int argc, char **argv) {
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

double distance(Body b1, Body b2) {
    return sqrt(pow(b2.xPos - b1.xPos, 2) + pow(b2.yPos - b1.yPos, 2));
}

double force(Body b1, Body b2) {
    return (G * b1.mass * b2.mass) / pow(distance(b1, b2), 2);
}
