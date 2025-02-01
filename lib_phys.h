#ifndef PHYS_H
#define PHYS_H
#include "lib_core.h"

enum figure {BOX, BALL};

typedef struct {
	float distance;
	Vector2 point;
} Intersection;

typedef struct {
    float minX;
    float maxX;
    float minY;
    float maxY;
} Shadow;

typedef struct {
    bool isCollision;
    Vector2 cp;
    Vector2 normal;
} CollisionPoint;

typedef struct shape_t {
    enum figure typ;
    Vector2 location;
    Vector2 velocity;
    float angVelocity;
    Vector2 accel;
    float angAccel;
    float mass;
    float inertia;
    bool marked;
    float radius;
    Vector2 orientation;
    Vector2 vertices[5];
    void (*funcDraw)(struct shape_t*, float, Color);
    void (*funcUpdate)(struct shape_t*);
    void (*funcResetPos)(struct shape_t*, Vector2);
} Shape;

Intersection intersectLine(Vector2 start_a, Vector2 end_a, Vector2 start_b, Vector2 end_b);

Shape Box(float x, float y, float w, float h);
Shape Ball(float x, float y, float r);
Shadow createShadow(Shape* shape);

void _rotateBox(Shape *box, float angle);

void shapeDraw(Shape* shape, float thick, Color c);
void shapeUpdate(Shape* shape);
void applyForce(Shape* shape, Vector2 force, float angForce);
void applyFriction(Shape* shape);
void shapeResetPos(Shape* shape, Vector2 v);
Vector2 checkKick(Shape* shape);

CollisionPoint detectCollBox(Shape* boxA, Shape* boxB);
CollisionPoint detectCollBall(Shape* ballA, Shape* ballB);
CollisionPoint detectCollBallBox(Shape* ball, Shape* box);
void resolveCollBox(Shape* boxA, Shape* boxB, Vector2 cp, Vector2 normal);
void resolveCollBall(Shape* ballA, Shape* ballB, Vector2 normal);
void resolveCollBallBox(Shape* ball, Shape* box, Vector2 cp, Vector2 normal);
void checkColl(Shape* shapeA, Shape* shapeB);


#endif
