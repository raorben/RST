//
//  residenceRadii.c
//
//  Created by Irina Tolkova on 7/28/15.
//
//  Function for computing residence time and distance given x-position, y-position, and time arrays, and values for the radius and threshold.
//
//  When using, call residenceTD, which loops through residenceMetric for every pair of radius and threshold values.


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void residenceMetric(int *n, double *x, double *y, double *time, double *radius, double *threshold, double *resid_time, double *resid_dist);

void residenceTD(int *n, double *x, double *y, double *time, int *r, double *radius_vec, double *threshold_vec, int *i1, int *i2, double *all_resid_time, double *all_resid_dist);

double dist(double x1, double y1, double x2, double y2);

void printVec(int n, double *arr);

void printVecBreaks(int n_each, int num, double *arr);


 
int main() {
    int n = 10;
    
    double x[10] = {-784.136, -783.195, -781.982, -776.238, -769.72, -763.236, -756.572, -750.713, -745.29, -738.64};
    double y[10] = {4939.199, 4939.678, 4940.191, 4941.898, 4943.828, 4945.085, 4946.441, 4948.267, 4949.437, 4951.715};
    double time[10] = {0.000000, 6.066667, 12.116667, 18.266667, 25.083333, 32.566667, 38.533333, 44.550000, 50.550000, 56.533333};
    
    int r = 5;
    double radius[5] = {0, 1, 5, 10, 20};
    double threshold[5] = {0, 0, 0, 0, 0};
    
    double* all_resid_time = (double*) calloc(n * r, sizeof(double));
    double* all_resid_dist = (double*) calloc(n * r, sizeof(double));
    
    
    printf("Before...\n");
    printf("Residence distance:\n");
    printVec(n * r, &all_resid_dist[0]);
    printf("Residence time:\n");
    printVec(n * r, &all_resid_time[0]);
    
    int i1 = -1;
    int i2 = -1;
    
    residenceTD(&n, &x[0], &y[0], &time[0], &r, &radius[0], &threshold[0], &i1, &i2, all_resid_time, all_resid_dist);
    
    printf("-------------------------------------");
    printf("\nResulting arrays: \n");
    printf("x positions:\n");
    printVec(n, &x[0]);
    printf("y positions:\n");
    printVec(n, &y[0]);
    printf("time:\n");
    printVec(n, &time[0]);
    printf("Residence distance:\n");
    printVecBreaks(n, r, &all_resid_dist[0]);
    printf("Residence time:\n");
    printVecBreaks(n, r, &all_resid_time[0]);
    return 0;
}




void residenceTD(int *n, double *x, double *y, double *time, int *r, double *radius_vec, double *threshold_vec, int *i1, int *i2, double *all_resid_time, double *all_resid_dist) {
    int radius_ind;
    for (radius_ind = 0; radius_ind < *r; radius_ind++) {
        printf("### Calculating residence values for a radius of %.3f. ###\n", radius_vec[radius_ind]);
        residenceMetric(n, x, y, time, &radius_vec[radius_ind], &threshold_vec[radius_ind], &all_resid_time[*n * radius_ind], &all_resid_dist[*n * radius_ind]);
    }
}


void residenceMetric(int *n, double *x, double *y, double *time, double *radius, double *threshold, double *resid_time, double *resid_dist) {
    
    int i;
    int j;
    
    int i1;
    int i2;
    
    // set start index
    i = 0;
    double dist_diff = 0;
    while (dist_diff < *radius) {
        i = i + 1;
        dist_diff = dist(x[0], y[0], x[i], y[i]);
    }
    i1 = i;
    for (j = 0; j < i1; j++) {
        resid_time[j] = -1;
        resid_dist[j] = -1;
    }

    // set end index
    i = *n - 1;
    dist_diff = 0;
    while (dist_diff < *radius) {
        i = i - 1;
        dist_diff = dist(x[*n - 1], y[*n - 1], x[i], y[i]);
    }
    i2 = i;
    for (j = i2 + 1; j < *n; j++) {
        resid_time[j] = -1;
        resid_dist[j] = -1;
    }
    
    
    // calculating residence values
    
    for (i = i1; i <= i2; i++) {
        int *track = (int*) calloc(*n, sizeof(int));
        track[i] = 1;
        
        double x1 = x[i];
        double y1 = y[i];
        
        double distance, step, dist_outside, start_time, start_dist_index;
        
        /* points up from the trajectory */
        
        distance = 0;
        dist_outside = 0;
        
        start_time = time[i];
        start_dist_index = i;
        
        double resid_time_up = 0;
        double resid_dist_up = 0;
        
        j = i + 1;
        
        // printf(" Current iteration: %d\n\n", i);
        
        // printf(" Going up! \n \n");
        
        while (dist_outside <= *threshold && j < *n) {
            double x2 = x[j];
            double y2 = y[j];
            
            distance = dist(x1, y1, x2, y2);
            step = dist(x2, y2, x[j-1], y[j-1]);
            
            // printf("j is %d, distance is %f, step is %f \n", j, distance, step);
            
            if (distance < *radius) {
                if (dist_outside > 0) {
                    start_time = time[j];
                } else {
                    resid_dist_up = resid_dist_up + step;
                }
                dist_outside = 0;
                track[j] = 1;
            } else {
                if (dist_outside == 0) {
                    // printf("--- resid_time_up is %f, time[j-1] is %f, start_time is %f \n", resid_time_up, time[j-1], start_time);
                    resid_time_up = resid_time_up + time[j-1] - start_time ;
                }
                track[j] = 2;
                dist_outside = dist_outside + step;
            }
            
            // printf("--- resid_time_up is %f, resid_dist_up is %f\n", resid_time_up, resid_dist_up);
            
            j = j + 1;
        }
        
        
        /* points down from the trajectory */
        
        distance = 0;
        dist_outside = 0;
        
        start_time = time[i];
        start_dist_index = i;
        
        double resid_time_down = 0;
        double resid_dist_down = 0;
        
        j = i - 1;
        
        // printf("\n Going down! \n \n");
        
        while (dist_outside <= *threshold && j >= 0) {
            double x2 = x[j];
            double y2 = y[j];
            
            distance = dist(x1, y1, x2, y2);
            step = dist(x2, y2, x[j+1], y[j+1]);
            
            // printf("j is %d, distance is %f, step is %f \n", j, distance, step);
            
            if (distance < *radius) {
                if (dist_outside > 0) {
                    start_time = time[j];
                } else {
                    resid_dist_down = resid_dist_down + step;
                }
                dist_outside = 0;
                track[j] = 1;
            } else {
                if (dist_outside == 0) {
                    resid_time_down = resid_time_down + start_time - time[j+1];
                }
                track[j] = 2;
                dist_outside = dist_outside + step;
            }
            j = j - 1;
            
            // printf("--- resid_time_down is %f, resid_dist_down is %f\n", resid_time_down, resid_dist_down);
        }

        resid_time[i] = resid_time_up + resid_time_down;
        resid_dist[i] = resid_dist_up + resid_dist_down;

        free(track);
    }
}


double dist(double x1, double y1, double x2, double y2) {
    return(sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2)));
}

void printVec(int n, double *arr) {
    int i;
    for (i = 0; i < n; i++) {
        printf("%f, ", arr[i]);
    }
    printf("\n\n");
}


void printVecBreaks(int n_each, int num, double *arr) {
    int i;
    int j;
    for (i = 0; i < num; i++) {
        for (j = 0; j < n_each; j++) {
            printf("%f, ", arr[n_each * i + j]);
        }
        printf("\n --- \n");
    }
    printf("\n\n");
}