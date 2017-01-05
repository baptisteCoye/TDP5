#ifndef DATATYPE_H
#define DATATYPE_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct {
  double m;
  double px;
  double py;
  double vx;
  double vy;
  double ax;
  double ay;
} particule;

void print_particule(particule p){
  printf("%lf %lf %lf %lf %lf %lf %lf\n", p.m, p.px, p.py, p.vx, p.vy, p.ax, p.ay);
}

typedef struct {
  double x;
  double y;
} vecteur;

typedef struct {
  particule * p;
  int h;
  int * begin;
  int * tailles_blocs;
  int max_part;
} quadtree;

void allocate_quadtree(quadtree* q, int h, int max_part){
  q->h = h;
  q->begin = malloc(sizeof(int) * (h+2));

  q->begin[0] = 0;
  for (int i = 1; i <= h+1; i++){
    q->begin[i] = q->begin[i-1] + pow(4, i-1);
  }
  q->max_part = max_part;

  q->p = malloc(sizeof(particule) * (q->begin[h+1] + max_part*pow(4,h)));
  q->tailles_blocs = malloc(sizeof(int) * pow(4,h));
}

void destroy_quadtree(quadtree* q){
  free(q->p); free(q->begin); free(q->tailles_blocs);
  q->h = -1; 
  q->max_part = -1;
}

particule BAD_PARTICLE = {-1,-1,-1,-1,-1,-1,-1};
vecteur BAD_VECTEUR = {-1,-1};

#endif
