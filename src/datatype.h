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

#endif
