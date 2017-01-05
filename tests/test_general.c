#include <stdio.h>
#include <stdlib.h>
#include <perf.h>
#include <math.h>

#include "datatype.h"
#include "util.h"
#include "TDPConfig.h"

#define FILENAME "../tests/data/file_test2.txt"
#define NB_ITER 100
#define NBPART_PER_BLOCK 1
#define MAX_SIZE_BLOCK 2000
#define h 2


int main(int argc, char ** argv){
  int nb_blocks = pow(4,h);
  int n = NBPART_PER_BLOCK;
  int N = n*nb_blocks;
  vecteur * centres = malloc(sizeof(vecteur) * nb_blocks);
  
  quadtree q;
  allocate_quadtree(&q, h, MAX_SIZE_BLOCK);
  
  fill(&(q.p[q.begin[h+1]]), nb_blocks*MAX_SIZE_BLOCK, N, 0, 0, 100, 100, h, centres, nb_blocks);

  for (int j = 0; j < nb_blocks; j++){
    q.p[q.begin[h] + j] = P2M(&(q.p[q.begin[h+1] + MAX_SIZE_BLOCK*j]), n);
    q.tailles_blocs[j] = n;
  }
  
  for (int i = h-1; i >= 0; i--){
    for (int j = 0; j < pow(4,i); j++){
      q.p[q.begin[i] + j] = M2M(&(q.p[q.begin[i+1] + j]), 4);
    }
  }
  
  vecteur * force = malloc(sizeof(vecteur) * N);
  double * distMin = malloc(sizeof(double) * N);
  double dt = DT_MAX, dtmp;
   
  for(int i = 0; i < nb_blocks; i++){
    print_particule(q.p[q.begin[h+1]+2000*i]);
    }
  printf("\n");
  for(int i = 0; i < 4; i++){
    print_particule(q.p[q.begin[h-1]+i]);
  }
  print_particule(q.p[q.begin[h-2]]);
  // Initialisation des forces à 0 et des distances à -1
  for (int i = 0; i < N; i++){
    force[i].x = 0.0;
    force[i].y = 0.0;
    distMin[i] = -1.0;
  }

  int tmpf = 0;
  for (int i = 0; i < nb_blocks; i++){
    printf("bloc numéro: %d\n", i);
    rec_calc(&(q.p[q.begin[h+1] + MAX_SIZE_BLOCK*i]), &(force[tmpf]), q.tailles_blocs[i], &(distMin[tmpf]), centres[i], q, 0, 0, 100);
    
    accelerate(&(q.p[q.begin[h+1] + MAX_SIZE_BLOCK*i]),&(force[tmpf]),n);
    dtmp = determine_dt_forall(&(q.p[q.begin[h+1] + MAX_SIZE_BLOCK*i]), &(force[tmpf]), n, &(distMin[tmpf]), 1);   
    if (dtmp < dt){
      dt = dtmp;
    }

    move_particules(&(q.p[q.begin[h+1] + MAX_SIZE_BLOCK*i]), &(force[tmpf]), n, dt);
    
    tmpf += q.tailles_blocs[i];
  }
}
