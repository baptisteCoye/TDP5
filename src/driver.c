#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <time.h>
#include "datatype.h"
#include "util.h"
#include "perf.h"
#include "TDPConfig.h"
#define TAG 100

int main(int argc, char **argv){
  clock_t begin = clock();
  float timebegin = (float) begin / CLOCKS_PER_SEC;
  float realbeg;


  /*  les arguments doivent etre argv[1] = nb_iterations & argv[2] = fileName  */
  if(argc != 3){
    fprintf(stderr, "Le programme doit contenir exactement deux arguments : le nombre desire d'iterations puis le nom du fichier ou recuperer les donnees.\n");
    return EXIT_FAILURE;
  }

  //////////////////////////////////////////////////////////////////
  ///           Initialisation des variables et du MPI           ///
  //////////////////////////////////////////////////////////////////

  int err;
  int nbProc, rank;
  particule * data;
  particule * buffer[2];
  vecteur * force;
  int nbParticules;
  int tag;
  MPI_Status status;
  MPI_Request sendRequests[2];
  MPI_Request recvRequests[2];
  double * distMin;
  double dt;

  MPI_Init(NULL,NULL);

  MPI_Comm_size(MPI_COMM_WORLD, &nbProc);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);


  MPI_Reduce(&timebegin, &realbeg, 1, MPI_FLOAT,
	     MPI_MIN, 0, MPI_COMM_WORLD);


  //////////////////////////////////////////////////////////////////
  ///             Creation du MPI_Datatype PARTICULE             ///
  //////////////////////////////////////////////////////////////////

  particule part[1];

  MPI_Datatype types[7] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
  int blocklengths[7] = {1,1,1,1,1,1,1};

  MPI_Aint disps[7];

  MPI_Aint i1,i2;

  MPI_Get_address(&part[0], &i1);

  MPI_Get_address(&part[0].m, &i2); disps[0] = i2-i1;

  MPI_Get_address(&part[0].px, &i2); disps[1] = i2-i1;
  MPI_Get_address(&part[0].py, &i2); disps[2] = i2-i1;

  MPI_Get_address(&part[0].vx, &i2); disps[3] = i2-i1;
  MPI_Get_address(&part[0].vy, &i2); disps[4] = i2-i1;

  MPI_Get_address(&part[0].ax, &i2); disps[5] = i2-i1;
  MPI_Get_address(&part[0].ay, &i2); disps[6] = i2-i1;
  
  MPI_Datatype PARTICULE;

  err =  MPI_Type_struct(7, blocklengths, disps, types, &PARTICULE);
  if (err != 0){
    fprintf(stderr, "erreur lors du MPI_Type_struct() : err = %d", err);
    return EXIT_FAILURE;
  }

  err = MPI_Type_commit(&PARTICULE);
  if (err != 0){
    fprintf(stderr, "erreur lors du MPI_Type_commit() : err = %d", err);
    return EXIT_FAILURE;
  }

  //////////////////////////////////////////////////////////////////
  ///           Creation du MPI_Datatype V_PARTICULE             ///
  //////////////////////////////////////////////////////////////////

  particule p[1];

  MPI_Datatype type[3] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
  int blocklen[3] = {1,1,1};

  MPI_Aint disp[3];

  MPI_Get_address(&p[0], &i1);
  MPI_Get_address(&p[0].m, &i2); disp[0] = i2-i1;
  MPI_Get_address(&p[0].px, &i2); disp[1] = i2-i1;
  MPI_Get_address(&p[0].py, &i2); disp[2] = i2-i1;
  
  MPI_Datatype V_PARTICULE;

  err =  MPI_Type_struct(3, blocklen, disp, type, &V_PARTICULE);
  if (err != 0){
    fprintf(stderr, "erreur lors du MPI_Type_struct() : err = %d", err);
    return EXIT_FAILURE;
  }

  err = MPI_Type_commit(&V_PARTICULE);
  if (err != 0){
    fprintf(stderr, "erreur lors du MPI_Type_commit() : err = %d", err);
    return EXIT_FAILURE;
  }

  /////////////////////////////////////////////////////////////////
  /// Lecture du fichier et allocation des tableaux et buffers  ///
  /////////////////////////////////////////////////////////////////

  err = readData(argv[2], nbProc, rank, &data, &nbParticules, PARTICULE);
  if (err != 0){
    fprintf(stderr, "erreur lors de la lecture du fichier.\n");
    return EXIT_FAILURE;
  }

  const int nbPartPerProc = nbParticules / nbProc;

  buffer[0] = malloc(sizeof(particule) * nbPartPerProc);
  buffer[1] = malloc(sizeof(particule) * nbPartPerProc);
  force = malloc(sizeof(vecteur) * nbPartPerProc);
  distMin = malloc(sizeof(double) * nbPartPerProc);

  /////////////////////////////////////////////////////////////////
  ///                  Affichage des input                      ///
  /////////////////////////////////////////////////////////////////

#if VERBOSE >= 1
  printf(":%d: nombre de particules dans le fichier : %d\n", rank, nbParticules);
  printf(":%d: nombre de processus actifs : %d\n", rank, nbProc);
  
  printf(":%d: contient %d testicules\n", rank, nbPartPerProc);
#endif

#ifdef SAVE_RESULTS

#if VERBOSE >= 1
  printf("save first result\n");
#endif

  char * initfilename = malloc(sizeof(char) * 20);
  snprintf(initfilename, 20, "save_par_0.csv");
    
  err = save_results(data, nbPartPerProc, initfilename, nbProc, rank, PARTICULE);
  if (err != 0){
    fprintf(stderr, "erreur lors de la sauvegarde du fichier.\n");
    return EXIT_FAILURE;
  }
  free(initfilename);
  
#endif /* TDP2_SAVE_RESULTS */

  /////////////////////////////////////////////////////////////////
  ///               Creation des requetes MPI                   ///
  /////////////////////////////////////////////////////////////////

  if (nbProc > 1){

    // Envoi depuis le buffer 0
    err = MPI_Send_init(buffer[0], nbPartPerProc, V_PARTICULE, (rank+1)%nbProc, TAG, MPI_COMM_WORLD, &sendRequests[0]);
    if (err != 0){
      fprintf(stderr, "Erreur detectee dans MPI_Send_init 0. err = %d\n", err);
      return EXIT_FAILURE;
    }
    // Envoi depuis le buffer 1
    err = MPI_Send_init(buffer[1], nbPartPerProc, V_PARTICULE, (rank+1)%nbProc, TAG, MPI_COMM_WORLD, &sendRequests[1]);
    if (err != 0){
      fprintf(stderr, "Erreur detectee dans MPI_Send_init 1. err = %d\n", err);
      return EXIT_FAILURE;
    }

    // Réception sur le buffer 0
    err = MPI_Recv_init(buffer[0], nbPartPerProc, V_PARTICULE, (rank-1)%nbProc, TAG, MPI_COMM_WORLD, &recvRequests[0]);
    if (err != 0){
      fprintf(stderr, "Erreur detectee dans MPI_Recv_init 0. err = %d\n", err);
      return EXIT_FAILURE;
    }

    // Réception sur le buffer 1
    err = MPI_Recv_init(buffer[1], nbPartPerProc, V_PARTICULE, (rank-1)%nbProc, TAG, MPI_COMM_WORLD, &recvRequests[1]);
    if (err != 0){
      fprintf(stderr, "Erreur detectee dans MPI_Recv_init 1. err = %d\n", err);
      return EXIT_FAILURE;
    }
  }

  ///////////////////////////////////////////////////////////////////
  ///              Debut des iterations de calcul                 ///
  ///////////////////////////////////////////////////////////////////

  for (int n = 0; n < atoi(argv[1]); n++){
#if VERBOSE >= 1
    if (rank == 0)
      printf("iteration %d\n", n);
#endif

    // Initialisation des forces a 0 et des distances à -1
    for (int i = 0; i < nbPartPerProc; i++){
      force[i].x = 0.0;
      force[i].y = 0.0;
      distMin[i] = -1.0;
    }

    if (nbProc > 1){
      // On copie les donnees possedees par le processus dans son buffer.
      copier(buffer[0], data, nbPartPerProc);

      // On demarre les requetes MPI pour l'envoi et la reception des donnees.
      err = MPI_Start(&sendRequests[0]);
      if (err != 0){
	fprintf(stderr, "Erreur detectee dans MPI_start_Send 0. err = %d\n", err);
	return EXIT_FAILURE;
      }
      err = MPI_Start(&recvRequests[1]);
      if (err != 0){
	fprintf(stderr, "Erreur detectee dans MPI_start_recv 1. err = %d\n", err);
	return EXIT_FAILURE;
      }
    }

#if VERBOSE >= 1
    printf("     :%d: calcul des forces ...\n", rank);
#endif

    // Pendant que le MPI gere l'envoi des donnees, on calcule
    // les forces entre les particules dans data et elles memes.
    calcul_local(force, data, nbPartPerProc, distMin);

    if (nbProc > 1){
      // Une fois que le calcul est fini, on attend la fin des echanges.
      err = MPI_Wait(&sendRequests[0], &status);
      if (err != 0){
	fprintf(stderr, "Erreur detectee dans MPI_wait_Send 0. err = %d\n", err);
	return EXIT_FAILURE;
      }
      err = MPI_Wait(&recvRequests[1], &status);
      if (err != 0){
	fprintf(stderr, "Erreur detectee dans MPI_wait_recv 1. err = %d\n", err);
	return EXIT_FAILURE;
      }
    }

    // Debut des iterations de calcul sur les donnees des autres processus.
    for (int i = 1; i < nbProc; i++){

      // Demarrage des requetes MPI sur pour l'envoi et la reception des donnees.
      err = MPI_Start(&sendRequests[i%2]);
      if (err != 0){
	fprintf(stderr, "Erreur detectee dans MPI_start_Send %d. err = %d\n", i%2, err);
	return EXIT_FAILURE;
      }
      MPI_Start(&recvRequests[(i+1)%2]);
      if (err != 0){
	fprintf(stderr, "Erreur detectee dans MPI_start_recv %d. err = %d\n", (i+1)%2, err);
	return EXIT_FAILURE;
      }
            
      // Pendant que le MPI gere l'envoi des donnees, on calcule
      // les forces entre les particules de data et celles reçues precedemment.
      calcul_lointain(force, buffer[i%2], data, nbPartPerProc, distMin);

      // Une fois le calcul fini, on attend la fin des echanges.
      MPI_Wait(&sendRequests[i%2], &status);
      if (err != 0){
	fprintf(stderr, "Erreur detectee dans MPI_wait_send %d. err = %d\n", i%2, err);
	return EXIT_FAILURE;
      }
      MPI_Wait(&recvRequests[(i+1)%2], &status);
      if (err != 0){
	fprintf(stderr, "Erreur detectee dans MPI_wait_recv %d. err = %d\n", (i+1)%2, err);
	return EXIT_FAILURE;
      }
    }

#if VERBOSE >= 1
    printf("     :%d: determination du dt ... \n", rank);
#endif

#if VERBOSE >= 2
    printf("           distMin :");
    for (int i = 0; i < nbPartPerProc; i++)
      printf("  %lf", distMin[i]);
    printf("\n");
#endif

    accelerate(data,force,nbPartPerProc);

    dt = determine_dt_forall(data, force, nbPartPerProc, distMin, nbProc);

#if VERBOSE >= 2
    printf("     :%d: dt = %lf", rank, dt);
#endif

#if VERBOSE >= 1
    printf("     :%d: deplacement des particules\n", rank);
#endif

    // Au bout de nbProc-1 iterations, on a l'ensemble des forces.
    // On applique ces forces aux particules.
    move_particules(data,force, nbPartPerProc, dt);
    // On enregistre les resultats si necessaire.
#ifdef SAVE_RESULTS
    char * filename = malloc(sizeof(char) * 20);
    snprintf(filename, 20, "save_par_%d.csv", n+1);

#if VERBOSE >= 1
    printf("     :%d: sauvegarde des donnees\n", rank);
#endif /* VERBOSE >= 1 */

    err = save_results(data, nbPartPerProc, filename, nbProc, rank, PARTICULE);
    if (err != 0){
      fprintf(stderr, "erreur lors de la sauvegarde du fichier.\n");
      return EXIT_FAILURE;
    }

    free(filename);
#endif /* TDP2_SAVE_RESULTS */
  }

  /////////////////////////////////////////////////////////////////////
  ///                 Liberation des ressources                     ///
  /////////////////////////////////////////////////////////////////////

  // desallocation
  free(data);
  free(buffer[0]);
  free(buffer[1]);
  free(force);


  clock_t end = clock();
  float timeend = (float) end / CLOCKS_PER_SEC;
  float realend;
  MPI_Reduce(&timeend, &realend, 1, MPI_FLOAT,
	     MPI_MAX, 0, MPI_COMM_WORLD);

  if(rank == 0){
    printf("Temps d'execution: %f secondes\n", realend - realbeg);
  }

  // Fin du MPI
  MPI_Type_free(&PARTICULE);
  MPI_Type_free(&V_PARTICULE);
  MPI_Finalize();

  return 0;
}
