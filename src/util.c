#include"util.h"

#define G 6.67
#define teta 0.707
#define MAX_RAND 100

double distance(particule A, particule B){
  return sqrt((A.px - B.px)*(A.px - B.px) + (A.py - B.py)*(A.py - B.py));
}

double generate_random(double a, double b){
  //srand((unsigned int) time(NULL));
  return (rand()/(double)RAND_MAX)*(b-a)+a;
}

vecteur force_interaction(particule A, particule B, double* distanceMinTmp){
  vecteur u, res;
  double dAB = distance(A, B);
  *distanceMinTmp = dAB;
  u.x = (B.px-A.px)/dAB;
  u.y = (B.py-A.py)/dAB;

  res.x = (G * A.m * B.m / (dAB * dAB)) * u.x;
  res.y = (G * A.m * B.m / (dAB * dAB)) * u.y;

  return res;
}

int readData(char* filename, int nbProc, int myRank, particule ** data, int * nbPart, MPI_Datatype PARTICULE){
  int i, j, k; 
  MPI_Status status;
  FILE * file;
  int nbPartPerProc;
  int err;

  if (myRank == 0){
    file = fopen(filename, "r");
    if (file == NULL){
      fprintf(stderr, "erreur lors de l'ouverture du fichier.\n");
      return -1;
    }
    
    err = fscanf(file, "%d\n", nbPart);
    if (err == EOF)
      return -1;
    
    if (*nbPart < nbProc){
      fprintf(stderr, "Pas assez de particules pour occuper tous les processeurs. Diminuez le nombre de processeurs.\n");
      return -1;
    }
    if ((*nbPart % nbProc) != 0){
      fprintf(stderr, "Le nombre de particules n'est pas un multiple du nombre de processeurs, le cas n'est pas géré.\n");
      return -1;
    }
  }

  MPI_Bcast(nbPart, 1, MPI_INT, 0, MPI_COMM_WORLD);
  nbPartPerProc = *nbPart / nbProc;
  *(data) = malloc(sizeof(particule) * nbPartPerProc);

  if (myRank == 0){
    for (k = 0; k < nbProc-1; k++){

      for (i = 0; i < nbPartPerProc; i++){
	err = fscanf(file, "%lf %lf %lf %lf %lf\n", &((*data)[i].m), &((*data)[i].px), &((*data)[i].py), &((*data)[i].vx), &((*data)[i].vy));
	if (err == EOF)
	  return -1;
      }
      MPI_Send(*data, nbPartPerProc, PARTICULE, k+1, 100, MPI_COMM_WORLD);
    }

    for (i = 0; i < nbPartPerProc; i++){
      err = fscanf(file, "%lf %lf %lf %lf %lf\n", &((*data)[i].m), &((*data)[i].px), &((*data)[i].py), &((*data)[i].vx), &((*data)[i].vy));
      if (err == EOF)
	return -1;
    }
    fclose(file);    
  } else {
    MPI_Recv(*data, nbPartPerProc, PARTICULE, 0, 100, MPI_COMM_WORLD, &status);
  }


  return 0;
}

void P2P(vecteur* force, particule* data, int N, double* distMin){

  vecteur tmp;
  double distTmp;

  for (int i = 0; i < N; i++){
    for (int j = i+1; j < N; j++){
      tmp = force_interaction(data[i], data[j], &distTmp);

      if ((distTmp < distMin[i]) || (distMin[i] < 0)){
	distMin[i] = distTmp;
      }
      if ((distTmp < distMin[j]) || (distMin[j] < 0)){
	distMin[j] = distTmp;
      }

      
      force[i].x += tmp.x;
      force[i].y += tmp.y;
      force[j].x -= tmp.x;
      force[j].y -= tmp.y;      
      
    }
  }
}

void P2P_ext(vecteur* force, particule* buffer, particule* data, int N, double* distMin){
  vecteur tmp;
  double distTmp;

  for (int i = 0; i < N; i++){
    for (int j = 0; j < N; j++){
      tmp = force_interaction(data[i], buffer[j], &distTmp);
      force[i].x += tmp.x;
      force[i].y += tmp.y;
      if (distTmp < distMin[i] || distMin[i] == -1){
        distMin[i] = distTmp;
      }
    }
  }
}

void copier(particule* buffer, particule* data, int N){
  for (int i = 0; i < N; i++){
    buffer[i] = data[i];
  }
}

void accelerate(particule * data, vecteur * force, int N){
  for (int i = 0; i < N; i++){
    data[i].ax = force[i].x / data[i].m;
    data[i].ay = force[i].y / data[i].m;
  }
}

void move_particules(particule * data, vecteur * force, int N, double dt){
  for (int i = 0; i < N; i++){
    data[i].px += data[i].vx*dt + data[i].ax*dt*dt/2;
    data[i].py += data[i].vy*dt + data[i].ay*dt*dt/2;

    data[i].vx += data[i].ax*dt;
    data[i].vy += data[i].ay*dt;
  }
}

int save_results(particule * data, const int N, char * filename, int nbProc, int myRank, MPI_Datatype PARTICULE){
  int i, j, k; 
  MPI_Status status;
  FILE * file;
  particule * buffer;

  if (myRank == 0){
    buffer = malloc(sizeof(particule) * N);
    file = fopen(filename, "w+");
    if (file == NULL){
      fprintf(stderr, "erreur lors de l'ouverture du fichier.\n");
      return -1;
    }
    
    //    fprintf(file, "%d\n", N*nbProc);
    fprintf(file, "X,Y\n");

    for (i = 0; i < N; i++){
      //      fprintf(file, "%lf %lf %lf %lf %lf\n", data[i].m, data[i].px, data[i].py, data[i].vx, data[i].vy);
      fprintf(file, "%lf,%lf\n", data[i].px, data[i].py);
    }

    for (k = 0; k < nbProc - 1; k++){

      MPI_Recv(buffer, N, PARTICULE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
      for (i = 0; i < N; i++){
	//	fprintf(file, "%lf %lf %lf %lf %lf\n", buffer[i].m, buffer[i].px, buffer[i].py, buffer[i].vx, buffer[i].vy);
	fprintf(file, "%lf %lf\n", buffer[i].px, buffer[i].py);
      }

    }

    fclose(file);    
  } else {
    MPI_Send(data, N, PARTICULE, 0, 100, MPI_COMM_WORLD);
  }

  return 0;
}

double determine_rac_max(double a, double b, double c){
    double det = b*b - 4*a*c;
    if(det < 0){
        return 0;
    }
    else return ((-b + sqrt(det))/(2*a));
}

double determine_dt(particule data, vecteur force, double distMin){
    double dtx;
    double dty;

    dtx = max(determine_rac_max((data.ax/2),data.vx,-(0.1*distMin)), determine_rac_max(-(data.ax/2),-data.vx,-(0.1*distMin)));
    dty = max(determine_rac_max((data.ay/2),data.vy,-(0.1*distMin)), determine_rac_max(-(data.ay/2),-data.vy,-(0.1*distMin)));

    if(dtx < dty)
      return dtx;
    else 
      return dty;
}

double determine_dt_forall(particule* data, vecteur* force, int N, double* distMin, int nbProc){
    int i;
    double dt = DT_MAX;
    double dtTmp;
    double dtTot;

    for (i = 0; i < N; i++){
      dtTmp = determine_dt(data[i], force[i], distMin[i]/sqrt(2));
      if (dtTmp < dt){
	    dt = dtTmp;
      }
    }
    /* if (dt < DT_MIN) */
      /* dt = DT_MIN; */

    if (nbProc > 1)
      MPI_Allreduce((void *) &dt, (void *) &dtTot, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    else 
      dtTot = dt;

    return dtTot;
}


void fill_random_bloc(particule* bloc, int N, double x, double y, double dx, double dy){
  int i;  
  
  for(i = 0; i < N; i++){
    bloc[i].m = generate_random(1,100);
    bloc[i].px = generate_random(0,dx)+ x;
    bloc[i].py = generate_random(0,dy)+ y;
    bloc[i].vx = 0;//(double) (rand() % ((int) dx) - (dx/2));
    bloc[i].vy = 0;//(double) (rand() % ((int) dy) - (dy/2));
  }
}

void fill(particule * bloc, int size, int nb_part, double x, double y, double dx, double dy, int h, vecteur * centres, int nb){
  if (h == 0){
    fill_random_bloc(bloc, nb_part, x, y, dx, dy);
    centres->x = (x+dx/2);
    centres->y = (y+dy/2);
    return;
  }

  for (int i = 0; i < 4; i++){
    fill(&(bloc[i*size/4]), size/4, nb_part/4, x+((i%2)*dx/2), y+(1-i/2)*dy/2, dx/2, dy/2, h-1, &(centres[i*nb/4]), nb/4);
  }
}

void monopole(particule* bloc, int N, double * m){
  *m = 0;
  for(int i = 0; i < N; i++){
    *m += bloc[i].m;
  }
}

particule P2M(particule* bloc, int N){
  double m;
  particule resultat;

  monopole(bloc, N, &m);
  
  resultat.m = m;
  resultat.px = 0;
  resultat.py = 0;
  resultat.vx = 0;
  resultat.vy = 0;

  for(int i = 0; i < N; i++){
    resultat.px += (bloc[i].m)*(bloc[i].px);
    resultat.py += (bloc[i].m)*(bloc[i].py);    
  }
  resultat.px /= m;
  resultat.py /= m;
 
  return resultat;
}

double dist(vecteur A, particule B){
  return sqrt((A.x - B.px)*(A.x - B.px) + (A.y - B.py)*(A.y - B.py));
}

void rec_calc(particule * bloc, vecteur * force, int N, double * distMin,	      vecteur centre,
	      quadtree q,  int pos, int h, int size){
  printf("call to rec_calc : \n  bloc = %d, N = %d\n  centre = [%lf,%lf]\n  pos = %d\n  h = %d\n  size = %d\n", bloc, N, centre.x, centre.y, pos, h, size);

    printf("dist : %lf\nteta*size : %lf\n", dist(centre, q.p[pos]), teta*size);
    print_particule(q.p[pos]);
  if ((h != q.h) && (dist(centre, q.p[pos]) >= teta*size)){
    printf("DONE M2P\n");
    M2P(force, q.p[pos], bloc, N);
    return;
  }

  if (h == q.h){
    int r = q.begin[h+1] + q.max_part*(pos - q.begin[h]);
    
    // printf("%d %d\n",&(q.p[r]), &(bloc[0]));
    if (&(q.p[r]) == &(bloc[0])){
      printf("DONE P2P\n");
      P2P(force, bloc, N, distMin);
    }
    else{
      printf("DONE P2P_EXT\n");
      P2P_ext(force, &(q.p[r]), bloc, N, distMin);
    }
    return;
  }
  int nextIndexBegin = 4*pos + 1;
  for(int i = nextIndexBegin; i < nextIndexBegin + 4; i++){
    rec_calc(bloc, force, N, distMin, centre,
	     q, i, h+1, size/2);
  }
}

void M2P(vecteur* force, particule mp, particule* data, int N){
  vecteur tmp;
  double distTmp;
  for(int i = 0; i < N; i++){
    tmp = force_interaction(data[i], mp, &distTmp);
    force[i].x += tmp.x;
    force[i].y += tmp.y;
  }
}

particule M2M(particule * bloc, int N){
  return P2M(bloc, N);
}
