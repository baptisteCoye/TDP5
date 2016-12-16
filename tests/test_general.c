#include <stdio.h>
#include <stdlib.h>
#include <perf.h>
#include <math.h>

#include "datatype.h"
#include "util.h"
#include "TDPConfig.h"

#define FILENAME "../tests/data/file_test2.txt"
#define NB_ITER 100
#define NBPART_PER_BLOCK 1000
#define MAX_SIZE_BLOC 2000
#define h 3
#define teta 0.707

int main(int argc, char ** argv){
  int nb_blocks = pow(4,h);
  int n = NB_PART_PER_BLOCK;
  int N = n*nb_blocks;

  particule * total = malloc(sizeof(particule) * N);  
  vecteur * force = malloc(sizeof(vecteur) * N);
  double * distMin = malloc(sizeof(double) * N);
  double dt, dt1, dt2;

  fill_random_bloc(total, N, 0, 1600, 1600, 1600);

  printf("bloc 1 : \n");
  for (int i = 0; i < n; i++){
    print_particule(total[i]);
  }
  printf("bloc 2 : \n");
  for (int i = n; i < N; i++){
    print_particule(total[i]);
  }
 
  // Initialisation des forces a 0 et des distances à -1
  for (int i = 0; i < N; i++){
    force[i].x = 0.0;
    force[i].y = 0.0;
    distMin[i] = -1.0;
  }

  // P2P du bloc 1
  P2P(&force[0], &total[0], n, &distMin[0]);
  
  // P2P du bloc 2
  P2P(&force[n], &total[n], n, &distMin[n]);

  // P2M du bloc 1
  particule mp1 = P2M(&total[0], n);

  // P2M du bloc 2
  particule mp2 = P2M(&total[n], n);

  // M2P du bloc 1 vers le bloc 2
  M2P(&force[n], mp1, &total[n], n);

  // M2P du bloc 2 vers le bloc 1
  M2P(&force[0], mp2, &total[0], n);

  // accelere les particules du bloc 1
  accelerate(&total[0],&force[0],n);

  // accelere les particules du bloc 2
  accelerate(&total[n],&force[n],n);

  // determine le dt dans le bloc 1
  dt1 = determine_dt_forall(&total[0], &force[0], n, &distMin[0], 1);

  // determine le dt dans le bloc 2
  dt2 = determine_dt_forall(&total[n], &force[n], n, &distMin[n], 1);

  // calcule le bon dt
  dt = min(dt1, dt2);

  //bouge les particules dans le bloc 1
  move_particules(&total[0],&force[0], n, dt);

  //bouge les particules dans le bloc 2
  move_particules(&total[n],&force[n], n, dt);


  printf("bloc 1 : \n");
  for (int i = 0; i < n; i++){
    print_particule(total[i]);
  }
  printf("bloc 2 : \n");
  for (int i = n; i < N; i++){
    print_particule(total[i]);
  }
}
