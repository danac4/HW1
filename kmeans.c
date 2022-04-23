#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <math.h>

typedef struct{
    int n;
    int dim;
} Info;

typedef struct{
    int count;
    double* newPoints;
    double* centroid;
} Cluster;

typedef struct{
    double* vector;
} Point;


Cluster* createClusters(int dim, Point* points, int k) {
    int i, j;
    Cluster *clusters;
    Cluster clus;
    clusters = (Cluster*) calloc(k, sizeof(Cluster));
    if (clusters == NULL){ 
        printf("An Error Has Occurred\n");
        exit(1);
    }
    for (i = 0; i < k; i++) {
        int count = 0;
        double *newPoints = (double*) calloc(dim, sizeof(double));
        double *centroid = (double*) calloc(dim, sizeof(double));
        if (newPoints == NULL|| centroid == NULL){
            printf("An Error Has Occurred\n");
            exit(1);
        }
        for (j = 0; j < dim; j++) {
            centroid[j] = points[i].vector[j];
        }
        clus.count = count;
        clus.newPoints = newPoints;
        clus.centroid = centroid;
        *(clusters + i) = clus;
    }
    return clusters;
}

Point* createPoints(int dim, int n){
    int i;
    double *v;
    Point* points = (Point*)calloc(n, sizeof(Point));
    if (points == NULL){ 
        printf("An Error Has Occurred\n");
        exit(1);
    }
    for(i = 0; i < n; i++){
        v = (double*) calloc(dim,sizeof(double));
        if (v == NULL){ 
            printf("An Error Has Occurred\n");
            exit(1);
        }
        points[i].vector = v;
    }
    return points;
}

void processFile(Info info, Point* points,const char* file_path){ 
    int i,j;
    FILE* file = NULL;
    file = fopen(file_path, "r");
    if (file == NULL){ 
        printf("Invalid Input!\n");
        exit(1);
    }
    for(i = 0; i < info.n; i++){
        for(j = 0; j < info.dim; j++){
            double coordinate;
            fscanf(file, "%lf", &coordinate);
            points[i].vector[j] = coordinate;
            getc(file);
        }
    }
    fclose(file);
}

Info extractInfo(const char* file_path){
    char c;
    int dim, n;
    FILE* file = NULL;
    Info inf;
    dim = 1;
    n = 1;
    file = fopen(file_path, "r");
    if (file == NULL){ 
        printf("Invalid Input!\n");
        exit(1);
    }
    for (c = getc(file); c != '\n'; c = getc(file)){
        if (c == ','){ 
            dim++;
        }
    }
    for (c = getc(file); c != EOF; c = getc(file)){
        if (c == '\n'){ 
            n++;
        }
    }
    inf.dim = dim;
    inf.n = n;
    fclose(file);
    return inf;
}

double distance(Point* point1, const double* centroid, int dim){
    int i;
    double sum, p;
    sum = 0;
    for (i = 0; i < dim; i++){
        p = point1->vector[i] - *(centroid + i);
        sum += p*p;
    }
    return sqrt(sum);
}

int min_distance(Point* point, Cluster* clusters, int dim, int k){
    int min_index, i;
    double min_val, curr;
    min_index = 0;
    min_val = distance(point, clusters[0].centroid, dim);
    for (i = 1; i < k; i++){
        curr = distance(point, clusters[i].centroid, dim);
        if (curr < min_val){
            min_val = curr;
            min_index = i;
        }
    }
    return min_index;
}

double euclidean_norm(const double* vector, int dim){
    int i;
    double result,p;
    result = 0;
    for (i = 0; i < dim; i++){
        p = vector[i];
        result += p*p;
    }
    return sqrt(result);
}


void add_point(Point* point, Cluster* cluster, int dim){
    int i;
    for (i = 0; i < dim; i++) {
        cluster->newPoints[i] += point->vector[i];
    }
    cluster->count += 1;
}

int centroid_update(Cluster* cluster, int dim, double *tmp_vector){
    int has_changed, i, l;
    double norm_check;
    has_changed = 1;
    if (cluster->count == 0){
        return 1;
    }
    for (i = 0; i < dim; i++) {
        tmp_vector[i] = cluster->newPoints[i]/cluster->count;
    }
    norm_check = euclidean_norm(cluster->centroid, dim) - euclidean_norm(tmp_vector, dim);
    if (norm_check >= 0.001 || norm_check <= -0.001){
        has_changed = 0;
    }
    for (l = 0; l < dim; l++) {
        cluster->centroid[l] = tmp_vector[l];
        cluster->newPoints[l] = 0;
    }
    cluster-> count = 0;
    return has_changed;
}

int clusters_update(Cluster* clusters, int k, int dim) {
    int changed, i, epsilon_indicator;
    double *tmp_vector;
    changed = 1;
    tmp_vector = (double *) calloc(dim, sizeof(double));
    if (tmp_vector == NULL){ 
        printf("An Error Has Occurred\n");
        exit(1);
    }
    for (i = 0; i < k; i++) {
        epsilon_indicator = centroid_update(&clusters[i], dim, tmp_vector);
        changed = ((changed) && (epsilon_indicator));
    }
    free(tmp_vector);
    return changed;
}

int is_num(char* arg){
    while (*arg != '\0'){
        if (*arg < '0' || *arg > '9')
        {
            return 0;
        }
        arg++;
    }
    return 1;
}

void write_to_output(Cluster* clusters, const char* file_path, int k, int dim){ 
    int i, j;
    FILE* file;
    file = fopen(file_path, "w");
    if (file == NULL){
        printf("Invalid Input!\n");
        exit(1);
    }
    for(i = 0; i < k; i++){
        for(j = 0; j < dim-1; j++){
            fprintf(file, "%0.4f,", clusters[i].centroid[j]);
        }
        fprintf(file, "%0.4f", clusters[i].centroid[dim-1]);
        fprintf(file, "%s", "\n");
    }
    fclose(file);
}

void kmeans(int max_iter, int n, Cluster* clusters, Point* points, int dim, int k){
    int epsilon_check, iter,i, index;
    epsilon_check = 0;
    iter = 0;
    while ((iter < max_iter) && (1 - epsilon_check)) {
        for (i = 0; i < n; i++) {
            index = min_distance(&points[i], clusters, dim, k);
            add_point(&points[i], &clusters[index], dim);
        }
        epsilon_check = clusters_update(clusters, k, dim);
        iter++;
    }

}

void free_memory(int k, int n, Point* points, Cluster* clusters){
    int i, j;
    for (i = 0; i < n; i++){
        free(points[i].vector);
    }
    free(points);
    for (j = 0; j < k; j++){
        free(clusters[j].newPoints);
        free(clusters[j].centroid);
    }
    free(clusters);

}
int main(int argc, char **argv) {
    char* verifier;
    char* file_path_input;
    char* file_path_output;
    int max_iter, dim, k, n;
    Info info;
    Point* points;
    Cluster* clusters;
    max_iter = 200;
    if (argc < 4 || argc > 5) {
        printf("Invalid Input!\n");
        return 1;
    }
    if (argc > 4) {
        verifier = argv[2];
        if (!is_num(verifier)){
            printf("Invalid Input!\n");
            return 1;
        }
        max_iter = atoi(argv[2]);
        if (!(max_iter > 0)){
            printf("Invalid Input!\n");
            return 1;
        }
        file_path_input = argv[3];
        file_path_output = argv[4];
    }
    else{
        file_path_input = argv[2];
        file_path_output = argv[3];
    }
    
    info = extractInfo(file_path_input);
    dim = info.dim;
    n = info.n;
    verifier = argv[1];
    if (!is_num(verifier)){
        printf("Invalid Input!\n");
        return 1;
    }

    k = atoi(argv[1]);
    if ((k < 1)||(k > n)){ 
        printf("Invalid Input!\n");
        return 1;
    }

    points = createPoints(dim, n);
    processFile(info, points, file_path_input);

    clusters = createClusters(dim, points, k);
    kmeans(max_iter, n, clusters, points, dim, k);

    write_to_output(clusters, file_path_output, k, dim);
    free_memory(k, n, points, clusters);
    return 0; 
}
