#include <math.h>
#include <mpi.h>
#include <iostream>
#include <string>

using namespace std;

const MPI_Comm comm = MPI_COMM_WORLD;
const int root = 0;
const double G = 6.67408 * pow(10, -11);
const double dt = 1;
const double softening = 0.00000001;

typedef struct Vec2_s {
    double x;
    double y;
} Vec2;

typedef struct Body_s {
    Vec2 pos;
    Vec2 velocity;
    double mass;
} Body;

double acceleration(double force, Body b);
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

            cin >> body.pos.x;
            cin >> body.pos.y;
            cin >> body.velocity.x;
            cin >> body.velocity.y;
            cin >> body.mass;

            bodies[i] = body;
        }
    } else {
    }

    MPI_Finalize();
    return 0;
}

double acceleration(double force, Body b) {
    return force / b.mass;
}

double distance(Body b1, Body b2) {
    const double dx = b2.pos.x - b1.pos.x;
    const double dy = b2.pos.y - b1.pos.y;
    return sqrt(pow(dx, 2) + pow(dy, 2));
}

double force(Body b1, Body b2) {
    const double dist = distance(b1, b2);
    return (G * b1.mass * b2.mass) / (pow(dist, 2) + pow(softening, 2));
}

Vec2 acceleration(Body b1, Body b2) {
    const double F = force(b1, b2);
    const Vec2 acc = {
        b2.pos.x - b1.pos.x * F,
        b2.pos.y - b1.pos.y * F,
    };

    return acc;
}
