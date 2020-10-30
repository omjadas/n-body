#include <math.h>
#include <mpi.h>
#include <iostream>
#include <string>

using namespace std;

const MPI_Comm comm = MPI_COMM_WORLD;
const int root = 0;
const double G = 6.67408 * pow(10, -11);
const double dt = 1;
const double softening = 1 * pow(10, -5);

struct Vec2 {
    double x;
    double y;
};

struct Body {
    int id;
    Vec2 pos;
    Vec2 velocity;
    double mass;
};

struct QuadTree {
    QuadTree *topLeft = nullptr;
    QuadTree *topRight = nullptr;
    QuadTree *bottomLeft = nullptr;
    QuadTree *bottomRight = nullptr;
    Body *body = nullptr;
    Vec2 anchor;
    double width;

    QuadTree(Vec2 a, double w) {
        anchor = a;
        width = w;
    }

    static QuadTree root() {
        return QuadTree((Vec2){
            0,
            0
        }, 100);
    }

    Vec2 getPos() {
        if (isLeaf()) {
            return body->pos;
        }

        return {};
    }

    void insert(Body *body) {
        if (isLeaf() && this->body == nullptr) {
            this->body = body;
            return;
        } else if (this->body != nullptr) {
            *topLeft = QuadTree((Vec2){
                anchor.x,
                anchor.y + width / 2,
            }, width / 2);
            *topRight = QuadTree((Vec2){
                anchor.x + width / 2,
                anchor.y + width / 2,
            }, width / 2);
            *bottomLeft = QuadTree((Vec2){
                anchor.x,
                anchor.y,
            }, width / 2);
            *bottomRight = QuadTree((Vec2){
                anchor.x + width / 2,
                anchor.y,
            }, width / 2);

            if (topLeft->contains(*this->body)) {
                topLeft->insert(this->body);
            } else if (topRight->contains(*this->body)) {
                topRight->insert(this->body);
            } else if (bottomLeft->contains(*this->body)) {
                bottomLeft->insert(this->body);
            } else if (bottomRight->contains(*this->body)) {
                bottomRight->insert(this->body);
            }

            if (topLeft->contains(*body)) {
                topLeft->insert(body);
            } else if (topRight->contains(*body)) {
                topRight->insert(body);
            } else if (bottomLeft->contains(*body)) {
                bottomLeft->insert(body);
            } else if (bottomRight->contains(*body)) {
                bottomRight->insert(body);
            }
        }
    }

    bool contains(Body body) {
        if (
            anchor.x < body.pos.x
            && (anchor.x + width) > body.pos.x
            && anchor.y < body.pos.y
            && (anchor.y + width) > body.pos.y
        ) {
            return true;
        }

        return false;
    }

    double getMass() {
        if (body != nullptr) {
            return body->mass;
        }

        double mass = 0;

        if (topLeft != nullptr) {
            mass += topLeft->getMass();
        }

        if (topRight != nullptr) {
            mass += topRight->getMass();
        }

        if (bottomLeft != nullptr) {
            mass += bottomLeft->getMass();
        }

        if (bottomRight != nullptr) {
            mass += bottomRight->getMass();
        }

        return mass;
    }

    bool isLeaf() {
        return topLeft == nullptr
            && topRight == nullptr
            && bottomLeft == nullptr
            && bottomRight == nullptr;
    }
};

Vec2 acceleration(double force, Body b);
double distance(Body b1, Body b2);
double force(Body b1, Body b2);
void doMPITask(int rank, int size);

int main(int argc, char **argv) {
    int rank;
    int size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    if (rank == root) {
        QuadTree qTree = QuadTree::root();

        int iterations;
        int n;

        cin >> iterations;
        cin >> n;

        MPI_Bcast(&iterations, 1, MPI_INT, root, comm);
        MPI_Bcast(&n, 1, MPI_INT, root, comm);

        Body bodies[n];

        for (int i = 0; i < n; i++) {
            Body body;

            body.id = i;
            cin >> body.pos.x;
            cin >> body.pos.y;
            cin >> body.velocity.x;
            cin >> body.velocity.y;
            cin >> body.mass;

            MPI_Bcast(&body, sizeof(Body), MPI_BYTE, root, comm);
        }

    } else {
        doMPITask(rank ,size);
    }

    MPI_Finalize();
    return 0;
}

Vec2 acceleration(Body b1, Body b2) {
    const double F = force(b1, b2);
    return {
        b2.pos.x - b1.pos.x * F,
        b2.pos.y - b1.pos.y * F,
    };
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

void doMPITask(int rank, int size) {
    int iterations;
    int n;

    MPI_Bcast(&iterations, 1, MPI_INT, root, comm);
    MPI_Bcast(&n, 1, MPI_INT, root, comm);

    QuadTree qTree = QuadTree((Vec2){ 0, 0 }, 100);

    for (int i = 0; i < n; i++) {
        Body body;
        MPI_Bcast(&body, sizeof(Body), MPI_BYTE, root, comm);
        qTree.insert(&body);
    }
}
