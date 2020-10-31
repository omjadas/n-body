#include <iostream>
#include <limits>
#include <math.h>
#include <mpi.h>
#include <sstream>
#include <string>
#include <sys/time.h>

using namespace std;

// Constants

const MPI_Comm comm = MPI_COMM_WORLD;
const int root = 0;
const double G = 6.67408 * pow(10, -11);
const double dt = 1 * pow(10, -2);
const double softening = 1 * pow(10, -5);
const double theta = 0.5;

/**
 * Struct representing a vector with two numeric elements.
 */
struct Vec2 {
    double x;
    double y;

    Vec2() {
        x = 0;
        y = 0;
    }

    Vec2(double x1, double y1) {
        x = x1;
        y = y1;
    }

    bool operator==(Vec2 v) {
        return x == v.x && y == v.y;
    }

    Vec2 operator+(Vec2 v) {
        return {
            x + v.x,
            y + v.y,
        };
    }

    Vec2 operator-(Vec2 v) {
        return {
            x - v.x,
            y - v.y,
        };
    }

    Vec2 operator*(double d) {
        return {
            x * d,
            y * d,
        };
    }

    Vec2 operator/(double d) {
        return {
            x / d,
            y / d,
        };
    }

    void operator+=(Vec2 v) {
        x += v.x;
        y += v.y;
    }

    void operator-=(Vec2 v) {
        x -= v.x;
        y -= v.y;
    }

    void operator*=(double d) {
        x *= d;
        y *= d;
    }

    void operator/=(double d) {
        x /= d;
        y /= d;
    }
};

/**
 * Struct representing a body in space.
 */
struct Body {
    int id;
    Vec2 pos;
    Vec2 velocity;
    double mass;

    Body() { }

    Body(int i, Vec2 p, Vec2 v, double m) {
        id = i;
        pos = p;
        velocity = v;
        mass = m;
    }

    Body(Vec2 p, Vec2 v, double m) {
        pos = p;
        velocity = v;
        mass = m;
    }

    /**
     * @return string representation of the body
     */
    string str() {
        ostringstream s;
        s << pos.x << " "
            << pos.y << " "
            << velocity.x << " "
            << velocity.y << " "
            << mass;
        return s.str();
    }
};

/**
 * Struct representing a node in a QuadTree.
 */
struct QuadTree {
    QuadTree *topLeft = nullptr;
    QuadTree *topRight = nullptr;
    QuadTree *bottomLeft = nullptr;
    QuadTree *bottomRight = nullptr;
    Body *body = nullptr;
    Vec2 anchor;
    double width;

    QuadTree(Vec2 *a, double w) {
        anchor = *a;
        width = w;
    }

    QuadTree() {
        anchor = Vec2(0, 0);
        width = 100;
    }

    QuadTree(double w) {
        anchor = Vec2(0, 0);
        width = w;
    }

    /**
     * @return the centre of mass for the current node in the QuadTree
     */
    Vec2 getPos() {
        if (isEmpty()) {
            return {
                anchor.x + width / 2,
                anchor.y + width / 2,
            };
        } else if (isLeaf() && body != nullptr) {
            return body->pos;
        }

        Vec2 pos = { 0, 0 };

        pos += topLeft->getPos() * topLeft->getMass();
        pos += topRight->getPos() * topRight->getMass();
        pos += bottomLeft->getPos() * bottomLeft->getMass();
        pos += bottomRight->getPos() * bottomRight->getMass();

        return pos / getMass();
    }

    /**
     * @return the centre of mass velocity for the current node in the QuadTree
     */
    Vec2 getVelocity() {
        if (isLeaf() && body != nullptr) {
            return body->velocity;
        } else if (isLeaf()) {
            return { 0, 0 };
        }

        Vec2 vel = { 0, 0 };
        vel += topLeft->getVelocity() * topLeft->getMass();
        vel += topRight->getVelocity() * topRight->getMass();
        vel += bottomLeft->getVelocity() * bottomLeft->getMass();
        vel += bottomRight->getVelocity() * bottomRight->getMass();
        return vel / getMass();
    }

    /**
     * @return the total mass of the current node of the QuadTree
     */
    double getMass() {
        if (isEmpty()) {
            return 0;
        } else if (isLeaf() && body != nullptr) {
            return body->mass;
        }

        double mass = 0;

        mass += topLeft->getMass();
        mass += topRight->getMass();
        mass += bottomLeft->getMass();
        mass += bottomRight->getMass();

        return mass;
    }

    /**
     * Insert a body into the QuadTree
     *
     * @param body pointer to the body to insert
     */
    void insert(Body *body) {
        if (isLeaf() && this->body == nullptr) {
            this->body = body;
            return;
        } else if (isLeaf() && this->body != nullptr) {
            topLeft = new QuadTree(new Vec2(
                anchor.x,
                anchor.y + width / 2
            ), width / 2);
            topRight = new QuadTree(new Vec2(
                anchor.x + width / 2,
                anchor.y + width / 2
            ), width / 2);
            bottomLeft = new QuadTree(new Vec2(
                anchor.x,
                anchor.y
            ), width / 2);
            bottomRight = new QuadTree(new Vec2(
                anchor.x + width / 2,
                anchor.y
            ), width / 2);

            if (topLeft->contains(*this->body)) {
                topLeft->insert(this->body);
            } else if (topRight->contains(*this->body)) {
                topRight->insert(this->body);
            } else if (bottomLeft->contains(*this->body)) {
                bottomLeft->insert(this->body);
            } else if (bottomRight->contains(*this->body)) {
                bottomRight->insert(this->body);
            }

            this->body = nullptr;
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

    /**
     * @return the current node of the QuadTree represented as a single body
     */
    Body asBody() {
        return Body(getPos(), getVelocity(), getMass());
    }

    /**
     * @param body body to determine if it is in the current node
     * @return whether the current node of the QuadTree geometrically contains a
     *         given body.
     */
    bool contains(Body body) {
        if (
            body.pos.x >= anchor.x
            && body.pos.x < (anchor.x + width)
            && body.pos.y >= anchor.y
            && body.pos.y < (anchor.y + width)
        ) {
            return true;
        }

        return false;
    }

    /**
     * @return whether the current node of the QuadTree is empty
     */
    bool isEmpty() {
        if (!isLeaf()) {
            return topLeft->isEmpty()
                && topRight->isEmpty()
                && bottomLeft->isEmpty()
                && bottomRight->isEmpty();
        }

        return body == nullptr;
    }

    /**
     * @return whether the current node of the QuadTree is a leaf
     */
    bool isLeaf() {
        return topLeft == nullptr
            && topRight == nullptr
            && bottomLeft == nullptr
            && bottomRight == nullptr;
    }
};

// Function declarations

Vec2 acceleration(Body b1, Body b2);
Vec2 acceleration(Body body, QuadTree qTree, Vec2 acc = { 0, 0 });
double distance(Body b1, Body b2);
double distance(Vec2 v1, Vec2 v2);
Vec2 force(Body b1, Body b2);
void updateBodies(int rank, int size, int iterations, int n, Body bodies[]);
Body updateBody(Body body, QuadTree qTree);
void recvBodies(int n, Body bodies[]);
void initialiseTree(int n, Body bodies[], QuadTree *qTree);
uint64_t GetTimeStamp();

/**
 * Main function.
 *
 * @param argc number of command line arguments
 * @param argv command line arguments
 * @return exit code
 */
int main(int argc, char **argv) {
    int rank;
    int size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    int iterations;
    int n;

    if (rank == root) {
        cin >> iterations;
        cin >> n;

        MPI_Bcast(&iterations, 1, MPI_INT, root, comm);
        MPI_Bcast(&n, 1, MPI_INT, root, comm);

        Body bodies[n];

        // Read bodies from stdin and broadcast them to other processes
        for (int i = 0; i < n; i++) {
            Body body;

            body.id = i;
            cin >> body.pos.x;
            cin >> body.pos.y;
            cin >> body.velocity.x;
            cin >> body.velocity.y;
            cin >> body.mass;

            bodies[i] = body;

            MPI_Bcast(&body, sizeof(Body), MPI_BYTE, root, comm);
        }

        uint64_t start = GetTimeStamp();

        updateBodies(rank, size, iterations, n, bodies);

        printf("Time: %ld us\n", (uint64_t)(GetTimeStamp() - start));

        for (int i = 0; i < n; i++) {
            cout << bodies[i].str() << endl;
        }
    } else {
        MPI_Bcast(&iterations, 1, MPI_INT, root, comm);
        MPI_Bcast(&n, 1, MPI_INT, root, comm);

        Body bodies[n];

        recvBodies(n, bodies);
        updateBodies(rank, size, iterations, n, bodies);
    }

    MPI_Finalize();
    return 0;
}

/**
 * Receive bodies from the root node.
 *
 * @param number of bodies to receive
 * @param bodies array to insert received bodies into
 */
void recvBodies(int n, Body bodies[]) {
    for (int i = 0; i < n; i++) {
        Body body;
        MPI_Bcast(&body, sizeof(Body), MPI_BYTE, root, comm);
        bodies[i] = body;
    }
}

/**
 * Initialise a QuadTree with the given bodies.
 *
 * @param n      number of bodies
 * @param bodies bodies to insert into the tree
 * @param qTree  pointer to the QuadTree to insert the bodies into
 */
void initialiseTree(int n, Body bodies[], QuadTree *qTree) {
    double minx = numeric_limits<double>::max();
    double miny = numeric_limits<double>::max();
    double maxx = numeric_limits<double>::min();
    double maxy = numeric_limits<double>::min();

    for (int j = 0; j < n; j++) {
        minx = min(minx, bodies[j].pos.x);
        miny = min(miny, bodies[j].pos.y);
        maxx = max(maxx, bodies[j].pos.x);
        maxy = max(maxy, bodies[j].pos.y);
    }

    *qTree = QuadTree(new Vec2(minx, miny),
                    (double) max(maxx - minx + 1, maxy - miny + 1));

    for (int j = 0; j < n; j++) {
        qTree->insert(&bodies[j]);
    }
}

/**
 * Update all bodies for the given number of iterations.
 *
 * @param rank       rank of the running process
 * @param size       number of processes running in the cluster
 * @param iterations number of iterations to calculate
 * @param n          number of bodies
 * @param bodies     bodies to simulate
 */
void updateBodies(int rank, int size, int iterations, int n, Body bodies[]) {
    for (int i = 0; i < iterations; i++) {
        Body updatedBodies[n];
        QuadTree qTree;
        initialiseTree(n, bodies, &qTree);

        for (int j = rank; j < n; j += size) {
            updatedBodies[j] = updateBody(bodies[j], qTree);
        }

        for (int j = rank; j < n; j += size) {
            bodies[j] = updatedBodies[j];
        }

        if (size > 1) {
            for (int j = 0; j < size; j++) {
                for (int k = j; k < n; k += size) {
                    MPI_Bcast(&bodies[k], sizeof(Body), MPI_BYTE, j, comm);
                }
            }
        }
    }
}

/**
 * Update a body given the current state of the QuadTree.
 *
 * @param body  body to update
 * @param qTree current state of the simulation
 * @return updated body
 */
Body updateBody(Body body, QuadTree qTree) {
    Vec2 acc = acceleration(body, qTree);
    Body newBody = Body(body.id, body.pos, body.velocity, body.mass);
    newBody.velocity += acc * dt;
    newBody.pos += newBody.velocity * dt;

    return newBody;
}

/**
 * Calculate the acceleration of a body due to another body.
 *
 * @param b1 body to calculate acceleration for
 * @param b2 other body
 * @return accelerator vector for b1
 */
Vec2 acceleration(Body b1, Body b2) {
    return force(b1, b2) / b1.mass;
}

/**
 * Recursively calculate the total acceleration of a body due to all other
 * bodies.
 *
 * @param body  to calculate acceleration for
 * @param qTree current state of the simulation
 * @param acc   acceleration of the body
 * @return acceleration vector
 */
Vec2 acceleration(Body body, QuadTree qTree, Vec2 acc) {
    if (qTree.isEmpty()) {
        return acc;
    } else if (qTree.isLeaf()) {
        if (body.id == qTree.body->id) {
            return acc;
        }
        acc += acceleration(body, *qTree.body);
    } else if ((qTree.width / distance(body.pos, qTree.getPos())) < theta) {
        acc += acceleration(body, qTree.asBody());
    } else {
        acc += acceleration(body, *qTree.topLeft, acc);
        acc += acceleration(body, *qTree.topRight, acc);
        acc += acceleration(body, *qTree.bottomLeft, acc);
        acc += acceleration(body, *qTree.bottomRight, acc);
    }

    return acc;
}

/**
 * Calculate the distance between two bodies.
 *
 * @param b1 first body
 * @param b2 second body
 * @return distance between the two bodies
 */
double distance(Body b1, Body b2) {
    return distance(b1.pos, b2.pos);
}

/**
 * Calculate the distance between two Vec2s.
 *
 * @param v1 first Vec2
 * @param v2 second Vec2
 * @return distance between the two Vec2s
 */
double distance(Vec2 v1, Vec2 v2) {
    const double dx = v1.x - v2.x;
    const double dy = v1.y - v2.y;
    return sqrt(pow(dx, 2) + pow(dy, 2));
}

/**
 * Calculate the force between two bodies.
 *
 * @param b1 first body
 * @param b2 second body
 * @return the force vector for the two bodies
 */
Vec2 force(Body b1, Body b2) {
    const double dist = distance(b1, b2);
    return (b2.pos - b1.pos) * (G * b1.mass * b2.mass) / (pow(dist, 2) + pow(softening, 2));
}

/**
 * @return current time
 */
uint64_t GetTimeStamp() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec * (uint64_t)1000000 + tv.tv_usec;
}
