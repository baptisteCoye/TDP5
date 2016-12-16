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
  int nb_feuilles;
  int nb_noeuds;
  int * tailles_blocs;
  int max_part;
} quadtree;

void allocate_quadtree(quadtree* q, int h, int max_part){
  int nb_feuilles = pow(4,h);
  q->nb_feuilles = nb_feuilles;
  q->h = h;
  q->max_part = max_part;

  int nb_noeuds = 0;
  int tmp = nb_feuilles;
  do {
    tmp /= 4;
    nb_noeuds += tmp;
  } while (tmp > 1);

  q->p = malloc(sizeof(particule) * (nb_noeuds + max_part*nb_feuilles));
  
  q->tailles_blocs = malloc(sizeof(int) * nb_feuilles));
}

void destroy_quadtree(quadtree* q){
  free(q->p); free(q->tailles_blocs);
  q->h = -1; 
  q->nb_feuilles = -1; q->nb_noeuds = -1;
  q->max_part = -1;
}

particule get_noeud(int hauteur, ){}


#endif
