#ifndef DATATYPE_H
#define DATATYPE_H

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
  q->begin = malloc(sizeof(int) * (h+1));

  q->begin[0] = 0;
  for (int i = 1; i <= h; i++){
    q->begin[i] = q->begin[i-1] + pow(4, i-1);
  }
  q->max_part = max_part;

  q->p = malloc(sizeof(particule) * (q->begin[h] + max_part*pow(4,h)));
  q->tailles_blocs = malloc(sizeof(int) * pow(4,h));
}

void destroy_quadtree(quadtree* q){
  free(q->p); free(q->begin); free(q->tailles_blocs);
  q->h = -1; 
  q->max_part = -1;
}

particule get_noeud(quadtree q, int hauteur, int i, int j){
  if ((hauteur >= q.h) || (hauteur < 0)){
    fprintf(stderr,"mauvaise profondeur dans le quadtree : demandé : %d, profondeur des feuilles : %d\n", hauteur, q.h);
    return p;
  }
  if ((i >= pow(2, hauteur)) || (j >= pow(2,h))){
    fprintf(stderr,"mauvais indices dans le quadtree, i = %d, j = %d, max = %d", i, j, pow(2, hauteur));
    return p;
  }

  return q.p[q.begin[h]+i*pow(2,h)+j];
}

particule * get_bloc(quadtree q, int i, int j){
  return &(q.p[q.begin[q.h] + 2000*(i*pow(2,h) + j)]);
}

#endif
