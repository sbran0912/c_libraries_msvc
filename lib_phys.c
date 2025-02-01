#include "lib_phys.h"

Intersection intersectLine(Vector2 start_a, Vector2 end_a, Vector2 start_b, Vector2 end_b) {
	Vector2 a = vec2_sub(end_a, start_a);
	Vector2 b = vec2_sub(end_b, start_b);
	float cross1 = vec2_cross(a, b);
	float cross2 = vec2_cross(b, a);
	if (fabs(cross1 - 0.0f) > 0.01) { //Float kann man nicht direkt auf 0 testen!!!
		float s = vec2_cross(vec2_sub(start_b, start_a), b) / cross1;
		float u = vec2_cross(vec2_sub(start_a, start_b), a) / cross2;
		if (s > 0.0001 && s < 1 && u > 0.0001 && u < 1) {
			return (Intersection){s, vec2_add(start_a, vec2_scale(a, s))};
		}
	}
	return (Intersection){0.0f, (Vector2){0.0f, 0.0f}};
}

void _drawBox(Shape *box, float thick, Color c) {
    for (int i = 0; i < 4; i++)
    {
        DrawLineEx(box->vertices[i], box->vertices[i + 1], thick, c);
    }
    DrawCircleV(box->location, 5, c);
}

void _rotateBox(Shape *box, float angle) {
    for (int i = 0; i < 5; i++)
    {
        box->vertices[i] = vec2_rotate(box->vertices[i], box->location, angle);
    }
}

void _updateBox(Shape *box) {
    box->velocity = vec2_add(box->velocity, box->accel);
    box->velocity =  vec2_limit(box->velocity, 10.0f);
    box->accel = (Vector2){0.0f, 0.0f};

    box->angVelocity += box->angAccel;
    box->angVelocity = limit_num(box->angVelocity, 0.05f);
    box->angAccel = 0.0f;

    box->location = vec2_add(box->location, box->velocity);
    box->vertices[0] = vec2_add(box->vertices[0], box->velocity);
    box->vertices[1] = vec2_add(box->vertices[1], box->velocity);
    box->vertices[2] = vec2_add(box->vertices[2], box->velocity);
    box->vertices[3] = vec2_add(box->vertices[3], box->velocity);
    box->vertices[4] = vec2_add(box->vertices[4], box->velocity);

    _rotateBox(box, box->angVelocity);
}

void _resetPosBox(Shape *box, Vector2 v) {
    if (box->mass != INFINITY)
    {
        box->location = vec2_add(box->location, v);
        box->vertices[0] = vec2_add(box->vertices[0], v);
        box->vertices[1] = vec2_add(box->vertices[1], v);
        box->vertices[2] = vec2_add(box->vertices[2], v);
        box->vertices[3] = vec2_add(box->vertices[3], v);
        box->vertices[4] = vec2_add(box->vertices[4], v);
    }
}

void _drawBall(Shape *ball, float thick, Color c) {
    DrawRing(ball->location, ball->radius - 3, ball->radius, 0, 360, 1, c);
    DrawCircleV(ball->location,3 , c);
    DrawLineEx(ball->location, ball->orientation, thick, c);
}

void _rotateBall(Shape* ball, float angle) {
    ball->orientation = vec2_rotate(ball->orientation, ball->location, angle);
}

void _updateBall(Shape* ball) {
    ball->velocity = vec2_add(ball->velocity, ball->accel);
    ball->velocity = vec2_limit(ball->velocity, 10.0f);
    ball->accel = (Vector2){0.0f, 0.0f};

    ball->angVelocity += ball->angAccel;
    ball->angVelocity = limit_num(ball->angVelocity, 0.05f);
    ball->angAccel = 0.0f;

    ball->location = vec2_add(ball->location, ball->velocity);
    ball->orientation = vec2_add(ball->orientation, ball->velocity);

    _rotateBall(ball, ball->angVelocity);
}

void _resetPosBall(Shape *ball, Vector2 v) {
    if (ball->mass != INFINITY) {
        ball->location = vec2_add(ball->location, v);
        ball->orientation = vec2_add(ball->orientation, v);
    }
}

Shadow createShadow(Shape* shape) {
    Shadow shadow;
    if (shape->typ == BALL) {
        shadow = (Shadow){
            .minX = shape->location.x - shape->radius, 
            .maxX = shape->location.x + shape->radius, 
            .minY = shape->location.y - shape->radius, 
            .maxY = shape->location.y + shape->radius
            };
    } else {
        shadow = (Shadow){
            .minX = INFINITY, 
            .maxX =-INFINITY, 
            .minY = INFINITY, 
            .maxY =-INFINITY
            };
        for (int i = 0; i < 4; i++) {
            if (shape->vertices[i].x < shadow.minX) {
                shadow.minX = shape->vertices[i].x;
            } 
            if (shape->vertices[i].y < shadow.minY) {
                shadow.minY = shape->vertices[i].y;
            } 
            if (shape->vertices[i].x > shadow.maxX) {
                shadow.maxX = shape->vertices[i].x;
            } 
            if (shape->vertices[i].y > shadow.maxY) {
                shadow.maxY = shape->vertices[i].y;
            } 
        }    
    }
    
    return shadow;
}

Shape Box(float x, float y, float w, float h) {
    Shape result = {
        .typ = BOX,
        .marked = false,
        .location = {x + w / 2, y + h / 2},
        .mass = (w + h) * 2,
        .inertia = w * h * w,
        .velocity = {0, 0},
        .angVelocity = 0,
        .accel = {0, 0},
        .angAccel = 0,
        .vertices[0] = {x, y},
        .vertices[1] = {x + w, y},
        .vertices[2] = {x + w, y + h},
        .vertices[3] = {x, y + h},
        .vertices[4] = {x, y},
        .funcDraw = &_drawBox,
        .funcUpdate = &_updateBox,
        .funcResetPos = &_resetPosBox
        };
    return result;
}

Shape Ball(float x, float y, float r) {
    Shape result = {
        .typ = BALL,
        .marked = false,
        .location = {x, y},
        .radius = r,
        .mass = r * 2,
        .inertia = r * r * r / 2,
        .velocity = {0, 0},
        .angVelocity = 0,
        .accel = {0, 0},
        .angAccel = 0,
        .orientation = {r + x, y},
        .funcDraw = &_drawBall,
        .funcUpdate = &_updateBall,
        .funcResetPos = &_resetPosBall
        };
    return result;
}

void shapeDraw(Shape *shape, float thick, Color c) {
    shape->funcDraw(shape, thick, c);
}

void shapeUpdate(Shape *shape) {
    shape->funcUpdate(shape);
}

void applyForce(Shape *shape, Vector2 force, float angForce) {
    shape->accel = vec2_add(shape->accel, vec2_div(force, shape->mass));
    shape->angAccel += angForce / shape->mass;
}

void applyFriction(Shape* shape) {
    float coefficient = 0.5f;
    Vector2 frictForce = vec2_norm(shape->velocity);
    frictForce = vec2_scale(frictForce, coefficient * -1); // in Gegenrichtung
    frictForce = vec2_limit(frictForce, vec2_mag(shape->velocity));

    float frictAngDirection = shape->angVelocity < 0 ? 1 : -1; // in Gegenrichtung
    float frictAngForce = limit_num(coefficient * 0.05 * frictAngDirection, abs(shape->angVelocity));

    applyForce(shape, frictForce, frictAngForce);
}

void shapeResetPos(Shape *shape, Vector2 v) {
    shape->funcResetPos(shape, v);
}

Vector2 checkKick(Shape* shape) {

    Vector2 mousePos = { (float)GetMouseX(), (float)GetMouseY() };

    if (IsMouseButtonDown(MOUSE_BUTTON_LEFT)) {

        if (vec2_dist(shape->location, mousePos) < 10) {
            shape->marked = true;
        }

        if (shape->marked) {
            DrawLineEx(shape->location, mousePos, 3, RED); //draw arrow
        }
    }

    if (IsMouseButtonReleased(MOUSE_BUTTON_LEFT) && shape->marked) {
        shape->marked = false;
        Vector2 force = vec2_sub(mousePos, shape->location);
        force = vec2_scale(force, 5);
        return force;
    }

    return (Vector2){0, 0};
}

CollisionPoint detectCollBox(Shape* boxA, Shape* boxB) {
    // Geprüft wird, ob eine Ecke von boxA in die Kante von boxB schneidet
    // Zusätzlich muss die Linie von Mittelpunkt boxA und Mittelpunkt boxB durch Kante von boxB gehen
    // i ist Index von Ecke und j ist Index von Kante
    // d = Diagonale von A.Mittelpunkt zu A.vertices(i)
    // e = Kante von B(j) zu B(j+1)
    // z = Linie von A.Mittelpunkt zu B.Mittelpunkt
    // _perp = Perpendicularvektor
    // scalar_d Faktor von d für den Schnittpunkt d/e
    // scalar_z Faktor von z für den Schnittpunkt z/e
    // mtv = minimal translation vector (überlappender Teil von d zur Kante e)

    for (int i = 0; i < 4; i++) {            
        for (int j = 0; j < 4; j++) {
            // Prüfung auf intersection von Diagonale d zu Kante e
            Intersection isd = intersectLine(boxA->location, boxA->vertices[i], boxB->vertices[j], boxB->vertices[j + 1]);   
            if (isd.distance > 0.0f) {
                // Prüfung auf intersection Linie z zu Kante e
                Intersection isz = intersectLine(boxA->location, boxB->location, boxB->vertices[j], boxB->vertices[j + 1]);
                if (isz.distance > 0.0f) {
                    // Collision findet statt
                    // Objekte zurücksetzen und normal_e berechnen. Kollisionspunkt ist Ecke i von BoxA
                    Vector2 e = vec2_sub(boxB->vertices[j + 1], boxB->vertices[j]);
                    Vector2 e_perp = {-(e.y), e.x};
                    Vector2 d = vec2_sub(boxA->vertices[i], boxA->location);  
                    d = vec2_scale(d, 1-isd.distance);
                    e_perp = vec2_norm(e_perp);
                    float distance = vec2_dot(e_perp, d);
                    e_perp = vec2_scale(e_perp, -distance); //mtv
                    shapeResetPos(boxA, vec2_scale(e_perp, 0.5f));
                    shapeResetPos(boxB, vec2_scale(e_perp, -0.5f));
                    e_perp = vec2_norm(e_perp); // normal_e
                    return (CollisionPoint){true, boxA->vertices[i], e_perp};                }
            }
        }
    }
    return (CollisionPoint) {false, {0.0f, 0.0f}, {0.0f, 0.0f}};
};

CollisionPoint detectCollBall(Shape* ballA, Shape* ballB) {
    //Distanz ermitteln
    float radiusTotal = ballA->radius + ballB->radius;
    float dist = vec2_dist(ballA->location, ballB->location);
    if (dist < radiusTotal) {
        //Treffer
        float space = (radiusTotal - dist);
        Vector2 collisionLine = vec2_sub(ballA->location, ballB->location);
        collisionLine = vec2_setMag(collisionLine, space);
        shapeResetPos(ballA, vec2_scale( collisionLine, 0.5));
        shapeResetPos(ballB, vec2_scale( collisionLine, -0.5));
        collisionLine = vec2_norm(collisionLine);
        return (CollisionPoint) {true, {0.0f, 0.0f}, collisionLine};
    }
    return (CollisionPoint) {false, {0.0f, 0.0f}, {0.0f, 0.0f}};
}

CollisionPoint detectCollBallBox(Shape* ball, Shape* box) {
    for (int i = 0; i < 4; i++) {
        //Kante der Box
        Vector2 e = vec2_sub(box->vertices[i+1], box->vertices[i]);
        //Vektor von Ecke der Box zum Ball
        Vector2 VerticeToBall = vec2_sub(ball->location, box->vertices[i]);
        //Kollision mit Ecken abfangen
        if (vec2_mag(VerticeToBall) < ball->radius) {
            return (CollisionPoint){true, box->vertices[i], VerticeToBall};
        }
        float mag_e = vec2_mag(e);
        e = vec2_norm(e);
        //Scalarprojektion von Vektor VerticeToBall auf Kante e
        float scalar_e = vec2_dot(VerticeToBall, e);
        if (scalar_e > 0 && scalar_e <= mag_e) {
            //Senkrechte von Ball trifft auf Kante e der Box
            //e2 = Kante e mit der Länge von scalar_e
            Vector2 e2 = vec2_scale(e, scalar_e);
            //Senkrechte von e zum Ball = VerticeToBall - e2
            Vector2 e_perp = vec2_sub(VerticeToBall, e2);

            if (vec2_mag(e_perp) < ball->radius) {
                //Ball berührt Box
                //Abstand wieder herstellen mit mtv (minimal translation vector)
                Vector2 mtv = e_perp;
                Vector2 p = vec2_add(box->vertices[i], e2);
                mtv = vec2_setMag(mtv, ball->radius - vec2_mag(e_perp));
                //e_perp und damit mtv zeigt von Kante zu Ball
                shapeResetPos(ball, mtv);
                //vor Berechnung muss e_perp normalisiert werden
                e_perp = vec2_norm(e_perp);
                //resolveCollisionBallBox(ball, box, p, e_perp)
                return (CollisionPoint){true, p, e_perp};
            }
        }
    }
    return (CollisionPoint) {false, {0.0f, 0.0f}, {0.0f, 0.0f}};
}

void resolveCollBox(Shape* boxA, Shape* boxB, Vector2 cp, Vector2 normal) {
    // rAP = Linie von A.location zu Kollisionspunkt (Ecke i von BoxA)
    Vector2 rAP = vec2_sub(cp, boxA->location);
    // rBP = Linie von B.location zu Kollisionspunkt (ebenfalls Ecke i von BoxA)
    Vector2 rBP = vec2_sub(cp, boxB->location);
    Vector2 rAP_perp = {-rAP.y, rAP.x};
    Vector2 rBP_perp = {-rBP.y, rBP.x};
    Vector2 VtanA = vec2_scale(rAP_perp, boxA->angVelocity);
    Vector2 VtanB = vec2_scale(rBP_perp, boxB->angVelocity);
    Vector2 VgesamtA = vec2_add(boxA->velocity, VtanA);
    Vector2 VgesamtB = vec2_add(boxB->velocity, VtanB);
    Vector2 velocity_AB = vec2_sub(VgesamtA, VgesamtB);
    if (vec2_dot(velocity_AB, normal) < 0) { // wenn negativ, dann auf Kollisionskurs
        float e = 1.0f; //inelastischer Stoß
        float j_denominator = vec2_dot(vec2_scale(velocity_AB, -(1+e)), normal);
        float j_divLinear = vec2_dot(normal, vec2_scale(normal, (1/boxA->mass + 1/boxB->mass)));
        float j_divAngular = (float)pow(vec2_dot(rAP_perp, normal), 2) / boxA->inertia + (float)pow(vec2_dot(rBP_perp, normal), 2) / boxB->inertia;
        float j = j_denominator / (j_divLinear + j_divAngular);
        // Grundlage für Friction berechnen (t)
        Vector2 t = {-normal.y, normal.x};
        float t_scalarprodukt = vec2_dot(velocity_AB, t);
        t = vec2_norm(vec2_scale(t, t_scalarprodukt));
        

        //apply Force        
        Vector2 force = vec2_add(vec2_scale(normal, (j/boxA->mass)), vec2_scale(t, (0.2*-j/boxA->mass)));
        float force_ang = vec2_dot(rAP_perp, vec2_add(vec2_scale(normal, j/boxA->inertia), vec2_scale(t, 0.2*-j/boxA->inertia)));
        boxA->accel = vec2_add(boxA->accel, force);
        boxA->angAccel += force_ang;

        force = vec2_add(vec2_scale(normal, (-j/boxB->mass)), vec2_scale(t, (0.2*j/boxB->mass)));
        force_ang = vec2_dot(rAP_perp, vec2_add(vec2_scale(normal, -j/boxB->inertia), vec2_scale(t, 0.2*j/boxB->inertia)));
        boxB->accel = vec2_add(boxB->accel, force);
        boxB->angAccel += force_ang;

    }
}

void resolveCollBall(Shape* ballA, Shape* ballB, Vector2 normal) {
    Vector2 rA = vec2_scale(normal, -ballA->radius);
    Vector2 rA_perp = {-rA.y, rA.x};
    Vector2 rB = vec2_scale(normal, ballB->radius);
    Vector2 rB_perp = {-rB.y, rB.x};
    Vector2 VtanA = vec2_scale(rA_perp, ballA->angVelocity);
    Vector2 VtanB = vec2_scale(rB_perp, ballB->angVelocity);
    Vector2 VgesamtA = vec2_add(ballA->velocity, VtanA);
    Vector2 VgesamtB = vec2_add(ballB->velocity, VtanB);
    Vector2 velocity_AB = vec2_sub(VgesamtA, VgesamtB);   
    
    if (vec2_dot(velocity_AB, normal) < 0) { // wenn negativ, dann auf Kollisionskurs
        float e = 1; //inelastischer Stoß
        float j_denominator = vec2_dot(vec2_scale(velocity_AB, -(1+e)), normal);
        float j_divLinear = vec2_dot(normal, vec2_scale(normal, (1/ballA->mass + 1/ballB->mass)));
        float j = j_denominator / j_divLinear;
        // Grundlage für Friction berechnen
        Vector2 t = {-normal.y, normal.x};
        float t_scalarprodukt = vec2_dot(velocity_AB, t);
        t = vec2_norm(vec2_scale(t, t_scalarprodukt));

        //apply Force
        Vector2 force = vec2_add(vec2_scale(normal, (0.8*j/ballA->mass)), vec2_scale(t, (0.2*-j/ballA->mass)));
        float force_ang = vec2_dot(rA_perp, vec2_scale(t, 0.1*-j/ballA->inertia));
        ballA->accel = vec2_add(ballA->accel, force);
        ballA->angAccel += force_ang;

        force = vec2_add(vec2_scale(normal, (0.8*-j/ballB->mass)), vec2_scale(t, (0.2*j/ballB->mass)));
        force_ang = vec2_dot(rB_perp, vec2_scale(t, 0.1*j/ballB->inertia));
        ballB->accel = vec2_add(ballB->accel, force);
        ballB->angAccel += force_ang;
    }
}

void resolveCollBallBox(Shape* ball, Shape* box, Vector2 cp, Vector2 normal) {
    Vector2 rA = vec2_scale(normal, -ball->radius);
    Vector2 rA_perp = {-rA.y, rA.x};
    Vector2 rBP = vec2_sub(cp, box->location);
    Vector2 rBP_perp = {-rBP.y, rBP.x};
    Vector2 VtanA = vec2_scale(rA_perp, ball->angVelocity);
    Vector2 VgesamtA = vec2_add(ball->velocity, VtanA);
    Vector2 VtanB = vec2_scale(rBP_perp, box->angVelocity);
    Vector2 VgesamtB = vec2_add(box->velocity, VtanB);
    Vector2 velocity_AB = vec2_sub(VgesamtA, VgesamtB);

    if (vec2_dot(velocity_AB, normal) < 0) { // wenn negativ, dann auf Kollisionskurs

        float e = 1; //inelastischer Stoß
        float j_denominator = vec2_dot(vec2_scale(velocity_AB, -(1+e)), normal);
        float j_divLinear = vec2_dot(normal, vec2_scale(normal, (1/ball->mass + 1/box->mass)));
        float j_divAngular = (float)pow(vec2_dot(rBP_perp, normal), 2) / box->inertia; //nur für Box zu rechnen
        float j = j_denominator / (j_divLinear + j_divAngular);
        // Grundlage für Friction berechnen
        Vector2 t = {-normal.y, normal.x};
        float t_scalarprodukt = vec2_dot(velocity_AB, t);
        t = vec2_norm(vec2_scale(t, t_scalarprodukt));
        
        // Apply Force
        Vector2 force = vec2_add(vec2_scale(normal, (0.8*j/ball->mass)), vec2_scale(t, (0.05*-j/ball->mass)));
        float force_ang = vec2_dot(rA_perp, vec2_scale(t, 0.05*-j/ball->inertia));
        ball->accel = vec2_add(ball->accel, force);
        ball->angAccel += force_ang;
           
        force = vec2_add(vec2_scale(normal, (-j/box->mass)), vec2_scale(t, (0.05*j/box->mass)));
        force_ang = vec2_dot(rBP_perp, vec2_add(vec2_scale(normal, -j/box->inertia), vec2_scale(t, 0.05*j/box->inertia)));
        box->accel = vec2_add(box->accel, force);
        box->angAccel += force_ang;
    }
}

void checkColl(Shape* shapeA, Shape* shapeB) {
    //Shadow berechnen von Element i und Element j 
    Shadow shadowA = createShadow(shapeA);
    Shadow shadowB = createShadow(shapeB);
    //Überschneidung prüfen
    if (shadowA.maxX >= shadowB.minX && shadowA.minX <= shadowB.maxX && shadowA.maxY >= shadowB.minY && shadowA.minY <= shadowB.maxY) {  
        //dann Überschneidung
        // Testcode
        DrawLineV(shapeA->location, shapeB->location, GREEN);
        // Ende Testcodew

        if (shapeA->typ == BALL) {
            if (shapeB->typ == BALL) {
                CollisionPoint cp = detectCollBall(shapeA, shapeB);
                if (cp.isCollision) {
                    resolveCollBall(shapeA, shapeB, cp.normal);
                }
            } else {
                CollisionPoint cp = detectCollBallBox(shapeA, shapeB);
                if (cp.isCollision) {
                    resolveCollBallBox(shapeA, shapeB, cp.cp, cp.normal);
                }
            }
        }

        if (shapeA->typ == BOX) {
            if (shapeB->typ == BOX) {
                // beide Boxen müssen geprüft werden, ob sie auf
                // die jeweils andere trefen könnte
                CollisionPoint cp = detectCollBox(shapeA, shapeB);
                if (cp.isCollision) {
                    resolveCollBox(shapeA, shapeB, cp.cp, cp.normal);  
                } else {
                    CollisionPoint cp = detectCollBox(shapeA, shapeB);    
                    if (cp.isCollision) {
                        resolveCollBox(shapeA, shapeB, cp.cp, cp.normal);
                    }
                }
            } else {
                CollisionPoint cp = detectCollBallBox(shapeB, shapeA);
                if (cp.isCollision) {
                    resolveCollBallBox(shapeB, shapeA, cp.cp, cp.normal);
                }
            }            
        }
    }
}