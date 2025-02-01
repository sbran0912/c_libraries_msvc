#include "lib_core.h"

typedef struct vehicle_t {
    Vector2 location;
    Vector2 velocity;
    Vector2 accel;
    Vector2 heading;
    float radius;
    float max_steer;
    float max_desired;
    float health;
    float dna[4];
} Vehicle;

typedef struct target_t {
    Vector2 location;
    float zone;
} Target;

Vehicle createVehicle(float x, float y, float r) {
    Vehicle result = {
        .location = {x, y},
        .radius = r,
        .velocity = {0, 0},
        .accel = {0, 0},
        .max_steer = 5.0f,
        .max_desired = 8.0f,
        .health = 1,
        .dna[0] = 1.0, //Force to Poison
        .dna[1] = 0.5, //Force to good Food
        .dna[2] = 200,  //Radius to Poison
        .dna[3] = 300   //Radius to good Food
        };
    return result;
} 

void applyForce(Vehicle *vehic, Vector2 force) {
    vehic->accel = vec2_add(vehic->accel, vec2_div(force, vehic->radius*2));
}

void seek(Vehicle* vehic, Vector2 target, bool is_badFood) {

    Vector2 desired = vec2_setMag(vec2_sub(target, vehic->location), vehic->max_desired);
    if (is_badFood) {
        desired = vec2_scale(desired, vehic->dna[0]);    
    } 
    else {
        desired = vec2_scale(desired, vehic->dna[1]);    
    }

    Vector2 steer = vec2_limit(vec2_sub(desired, vehic->velocity), vehic->max_steer); 
    applyForce(vehic, steer);
}

void updateVehicle(Vehicle* vehic) {
    vehic->velocity = vec2_limit(vec2_add(vehic->velocity, vehic->accel), 13.0f);
    vehic->accel = (Vector2){0.0f, 0.0f};
    vehic->location = vec2_add(vehic->location, vehic->velocity);
}

void drawVehicle(Vehicle* vehic) {

    Vector2 points[3];
    Vector2 back = vec2_sub(vehic->location, vec2_setMag(vehic->velocity, vehic->radius));
    Vector2 perp = vec2_setMag((Vector2){-(vehic->velocity.y), vehic->velocity.x}, vehic->radius/2);
    
    points[0] = vehic->location;
    points[1] = vec2_add(back, perp);
    points[2] = vec2_sub(back, perp);

    Color c = LIGHTGRAY;
    if (vehic->health < 0.3) {
        c = RED;
    }
    DrawTriangle(points[0], points[2],points[1], c);
    DrawRing(points[0], vehic->dna[2] - 3, vehic->dna[2], 0, 360, 1, RED);
    DrawRing(points[0], vehic->dna[3] - 3, vehic->dna[3], 0, 360, 1, GREEN);    
}

Vector2* createFood(size_t count) {
    Vector2* food = arr_create(count, sizeof(Vector2));
    for (size_t i = 0; i < count; i++) {
        Vector2 pos = {random_num(1, GetScreenWidth()), random_num(1, GetScreenHeight()) };
        food = arr_push(food, &pos);
        
    }
    //printf("cap: %lu len: %lu \n", arr_cap(food), arr_len(food));
    return food;
}

Vector2* respawnFood(Vector2* food, size_t min, size_t count) {
    if (arr_len(food) < min) {
        printf("Neues Futter ... \n");
        for (size_t i = 0; i < count; i++) {
            Vector2 pos = {random_num(1, GetScreenWidth()), random_num(1, GetScreenHeight()) };
            food = arr_push(food, &pos);
        } 
    }
    return food;
}

void drawFood(Vector2* food, Color c) {
    for (size_t i = 0; i < arr_len(food); i++) {   
        DrawCircleV(food[i], 4, c);
    }
    
}

void eatFood(Vehicle* vehic, Vector2* food, bool is_badFood) {
    float mindist = INFINITY;
    int idx = -1;
    
    for (size_t i = 0; i < arr_len(food); i++) {
        float filter = is_badFood ? vehic->dna[2] : vehic->dna[3];
        float distance = vec2_dist(vehic->location, food[i]);
        if ( distance < filter && distance < mindist) {
            mindist = distance;
            idx = i;
        }     
    }

    if (idx > -1) {
        if (vec2_dist(vehic->location, food[idx]) < 3) {
            arr_pop(food, idx);  //eat food
            if (is_badFood) {
                vehic->health -= 0.1;
            }
            else {
                vehic->health += 0.1;
            }
        } else {
            seek(vehic, food[idx], is_badFood);
        }
    }
}

int main() {
    const int screenWidth = 1600;
    const int screenHeight = 1200; 
    InitWindow(screenWidth, screenHeight, "seeking and arrive");
    SetTargetFPS(60);               // Set our game to run at 60 frames-per-second

    Vehicle mover = createVehicle(100, 400, 20);
    mover.velocity = (Vector2){5, -1};    
    Vector2* food = createFood(50);
    Vector2* poison = createFood(50);
    
    while (!WindowShouldClose()) {

        eatFood(&mover, food, false);
        eatFood(&mover, poison, true);
        updateVehicle(&mover);
        food = respawnFood(food,40, 10);
        poison = respawnFood(poison, 40, 80);
        
        BeginDrawing();                
        ClearBackground(DARKBLUE);
        drawVehicle(&mover);
        drawFood(food, GREEN);
        drawFood(poison, RED);
        EndDrawing();         
    }

    CloseWindow();   

    return EXIT_SUCCESS;
}
