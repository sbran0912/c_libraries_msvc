#include "lib_core.h"

int random_num(int start, int end) {

    unsigned int random_value = 0;
    if (rand_s(&random_value) != 0) {
        // Fehler beim Generieren einer Zufallszahl
        return -1;
    }

    // Berechnung der Zufallszahl im gewünschten Bereich
    return (random_value % (end - start + 1)) + start;
}


float limit_num(float number, float limit) {
	int vorzeichen = (number < 0) ? -1 : 1;
	float numberMag = number * vorzeichen;
	if (numberMag > limit) numberMag = limit;

	return numberMag * vorzeichen;
}

float constrain_num(float value, float min, float max) {
  return fminf(max, fmaxf(min, value));
}

float map_num(float n, float range_old, float range_new) {
  return (n / range_old * range_new);
}

ListStr lists_create(size_t capacity) {
    return (ListStr){.array = malloc(capacity * sizeof(char*)), .capacity = capacity, .count = 0};
}

void lists_append(ListStr* l, char val[]) {
    l->count++;
    if ((l->count) > (l->capacity)) {
        l->capacity = l->capacity + 10;
        l->array = realloc(l->array, l->capacity * sizeof(char*));
    }
    l->array[l->count-1] = malloc(strlen(val)+1);
    strcpy(l->array[l->count-1], val);
}

void lists_free(ListStr* l) {
    for (size_t i = 0; i < l->count; i++) {
        free(l->array[i]);
    }
    free(l->array);
    l->array = NULL;
}

void lists_reduce(ListStr* l) {    
    if (l->count < l->capacity) {
        l->capacity = l->count;   //Kapazität auf Länge setzen
        l->array = realloc(l->array, l->capacity * sizeof(char*));
    }
}

ListStr split(char string[], const char delimiter) {
    ListStr result = lists_create(10);
    char d;
    d = delimiter;

    char *token = strtok(string, d);
    while (token != NULL) {
        lists_append(&result, token);
        token = strtok(NULL, d);
    }
    return result;
}

ListFloat2d listf2d_create(size_t rows, size_t cols, float val) {
    float** arr = malloc(rows * sizeof(float*));
    for (size_t i = 0; i < rows; i++) {
        arr[i] = malloc(cols * sizeof(float));
        for (size_t j = 0; j < cols;j++) {
            arr[i][j] = val;
        }
    }
    return (ListFloat2d){.array = arr, .rows = rows, .cols = cols};
}

void listf2d_free(ListFloat2d* l) {
    for (size_t i = 0; i < l->rows; i++) {
        free(l->array[i]);
    }
    free(l->array);
    l->array = NULL;
}

void* arr_create(size_t capacity, size_t elemsize) {
    size_t *raw = malloc(3 * sizeof(size_t) + (capacity * elemsize));
    raw[0] = 0;                 // Number of elements currently in the array
    raw[1] = capacity;         // capacity - how many elements can the array currently store
    raw[2] = elemsize;           // Size of elements
    return (void*) &raw[3];     // Pointer to the beginning of the array
}

size_t arr_len(void* arr) {
    return *((size_t*)arr - 3);
}

size_t arr_cap(void* arr) {
    return *((size_t*)arr - 2);
}

void* arr_push(void* arr, void* val) {
    size_t* raw = (size_t*)(void*)arr - 3;
    size_t len = raw[0];
    size_t cap = raw[1];
    size_t elemsize = raw[2];

    len++;                  //increase counter
    raw[0] = len;           //save counter in array

    if (len > cap) {
        size_t new_cap = cap + 10;
        raw[1] = new_cap;      //calculate new capacity
        raw = realloc(raw, 3 * sizeof(size_t) + (new_cap * elemsize));
        arr = (void*)&raw[3]; //Array-Pointer wird wieder auf richtigen Startpunkt gestellt
    } 
    
    char* new_item = (char*)&raw[3] + ((len-1) * elemsize);   //count 1 = index 0, count 2 = index 1 etc.
    memcpy(new_item, val, elemsize);
    return arr;
}

void arr_pop(void* arr, size_t idx) {
    size_t* raw = (size_t*)(void*)arr - 3;
    size_t len = raw[0];
    //size_t cap = raw[1];
    size_t elemsize = raw[2];
    while (idx < len -1) {
        char* val = (char*)&raw[3] + ((idx + 1) * elemsize); //Pointer auf nächstes Element
        char* upd_item = (char*)&raw[3] + (idx * elemsize); //Pointer aktuelles Element
        memcpy(upd_item, val, elemsize); //Kopiere Wert vom nächsten Element auf aktuelles Element
        idx++;
    }
    len--;
    raw[0] = len;
}

void* arr_reduce(void* arr) {
    size_t* raw = (size_t*)(void*)arr - 3;
    size_t len = raw[0];
    size_t cap = raw[1];
    size_t elemsize = raw[2];
    
    if (len < cap) {
        raw[1] = len; //Kapazität auf Länge setzen
        raw = realloc(raw, 3 * sizeof(size_t) + (len * elemsize)); //Speicherbereich kürzen
        arr = (void*)&raw[3]; //Array-Pointer wird wieder auf richtigen Startpunkt gestellt
    }
    return arr;
}

void arr_free(void* arr) {
    size_t* raw = (size_t*)(void*)arr - 3;
    free(raw);
}

void llist_push(LinkedList* list, void* data, size_t datasize) {
    Node* new = malloc(sizeof(Node));
    new->data = malloc(sizeof(datasize));
    memcpy(new->data, data, datasize);
    new->next = NULL;
    
    if (list->size == 0) {
        list->head = new;
    } else {
        list->tail->next = new;
    }
    list->tail = new;
    list->size++;
    
}

void llist_pop(LinkedList* list) {
    if (list->size > 0) {
        Node* tmp = list->head->next;
        free(list->head->data);
        free(list->head);
        list->head = tmp;
        list->size--;

        if (list->size == 0) {
            list->tail = NULL;
        }
    }
}

Vector2 vec2_scale(Vector2 v, float n) {
	return (Vector2) {v.x * n, v.y * n};
}

Vector2 vec2_div(Vector2 v, float n) {
	Vector2 result = {v.x/n, v.y/n};
	return result;

}

Vector2 vec2_add(Vector2 v1, Vector2 v2) {
	Vector2 result = {v1.x + v2.x, v1.y + v2.y};
	return result;
}

Vector2 vec2_sub(Vector2 v1, Vector2 v2) {
	Vector2 result = {v1.x - v2.x, v1.y - v2.y};
	return result;
}

Vector2 vec2_norm(Vector2 v) {
    Vector2 result = {0};
    float len = sqrtf((v.x*v.x) + (v.y*v.y));
    if (len > 0.001) {
        result.x = v.x / len;
        result.y = v.y / len;
    }
    return result;	
}

Vector2 vec2_limit(Vector2 v, float max) {
	Vector2 result = v;
	float len = sqrtf(result.x * result.x + result.y * result.y);
	if (len > max) {
		result.x = result.x / len * max;
		result.y = result.y / len * max;
	}
	return result;
}

Vector2 vec2_setMag(Vector2 v, float magnitude) {
    Vector2 result = {0};
    float len = sqrtf((v.x*v.x) + (v.y*v.y));
    if (len > 0.001) {
        result.x = v.x / len * magnitude;
        result.y = v.y / len * magnitude;
    }
    return result;	

}

float vec2_dot(Vector2 v1, Vector2 v2) {
	return v1.x * v2.x + v1.y * v2.y;
}

float vec2_cross(Vector2 v1, Vector2 v2) {
	return v1.x * v2.y - v1.y * v2.x;
}

Vector2 vec2_perp(Vector2 v) {
	return (Vector2){-v.y, v.x};
}

float vec2_mag(Vector2 v) {
	return sqrtf(v.x * v.x + v.y * v.y);
}

float vec2_magsq(Vector2 v) {
	return v.x * v.x + v.y * v.y;
}

float vec2_dist(Vector2 v1, Vector2 v2) {
	Vector2 vdist = vec2_sub(v1, v2);
	return vec2_mag(vdist);
}

Vector2 vec2_rotate(Vector2 v, Vector2 base, float n) {
	Vector2 direction = vec2_sub(v, base);
	float x = direction.x * cosf(n) - direction.y * sinf(n);
	float y = direction.x * sinf(n) + direction.y * cosf(n);
	return (Vector2){x = x + base.x, y = y + base.y};
}

float vec2_heading(Vector2 v) {
	return atan2f(v.y, v.x);
}

float vec2_angle(Vector2 v1, Vector2 v2) {
    float result = 0.0f;
    float dot = v1.x*v2.x + v1.y*v2.y;
    float det = v1.x*v2.y - v1.y*v2.x;
    result = atan2f(det, dot);

    return result;
}

float vec2_minDist(Vector2 p, Vector2 start_a, Vector2 end_a) {
	float dist = -1.0f;

	//Vektor start_a to end_a (line_a)
	Vector2 line_a = vec2_sub(end_a, start_a);
	//Vektor imaginary line start_a to p
	Vector2 start_a_to_p = vec2_sub(p, start_a);
	//Magnitude of line_a
	float magnitude = vec2_mag(line_a);

	//Scalarprojecton from line (start_a to p) on line_a
	line_a = vec2_norm(line_a);
	float sp = vec2_dot(line_a, start_a_to_p);

	//Scalarprojection in magnitude of line_a?
	if (sp > 0.0001 && sp <= magnitude) {
		line_a= vec2_scale(line_a, sp);
		dist = vec2_mag(vec2_sub(start_a_to_p, line_a));
	}
	return dist;
}

Vector2 vec2_random() {
    float angle = (float) random_num(0, 100) / 100 * PI * 2;
    float x = cosf(angle);
    float y = sinf(angle);

    return (Vector2){x, y};
}