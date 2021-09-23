#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <errno.h>
#include <error.h>
#include <mpi.h>
#include "arralloc.h"
#include "write_ppm.h"
#include "read_options.h"

typedef int (*in_set_fn_t)(const float, const float, const float, const float, const int);

#define EMPTY_TASK -1
#define MIN(X, Y) ((X) < (Y) ? (X) : (Y))

static inline int point_in_mandelbrot_set(const float x0,
                                          const float y0,
					  const float cre,
					  const float cim,
                                          const int max_iter)
{
  /*
   * Note that cre and cim are not needed for Mandelbrot set.
   * They are only needed for Julia set, but both functions must have same
   * prototype.
   */

    int i;
    float x;
    float y;
    float x2;
    float y2;
    x = x0;
    y = y0;
    for ( i = 0; i < max_iter; i++ ) {
        x2 = x * x;
        y2 = y * y;
        /* z = (z*z) + c */
        if ( x2 + y2 > 4.0 ) {
            return i;
        } else {
            y = y0 + (2.0 * x * y);
            x = x0 + x2 - y2;
        }
    }
    return max_iter;
}

static inline int point_in_julia_set(const float x0,
                                     const float y0,
				     const float cre,
				     const float cim,
                                     const int max_iter)
{
    int i;
    float x;
    float y;

    x = x0;
    y = y0;
    for ( i = 0; i < max_iter; i++ ) {
        if ( x * x + y * y > 4.0 ) {
            return i;
        } else {
            float tmp = y * y;
            y = 2 * x * y + cim;
            x = x * x - tmp + cre;
        }
    }
    return max_iter;
}

static inline int point_in_new_julia_set(const float x0,
                                     const float y0,
                                     const int max_iter)
{
  /*
   * Uses a different function: z = z * exp(z) + c
   */

    int i;
    float x;
    float y;
    const float cim = 0.00;
    const float cre = 0.04;
    x = x0;
    y = y0;
    for ( i = 0; i < max_iter; i++ ) {
        if ( x * x + y * y > 4.0 ) {
            return i;
        } else {
	  float expx   = exp(x);
	  float expiyim = sin(y);
	  float expiyre = cos(y);
	  float ynew;

	  ynew = expx * (x*expiyim + y*expiyre);
	  x    = expx * (x*expiyre - y*expiyim) + 0.04;

	  y = ynew;
        }
    }
    return max_iter;
}

/* Initialise data for image array.  Data is stored in "scanline
 * order", i.e. x dimension varies fastest.  You get an array with
 * shape image[grid_size_y][grid_size_x] from this function. */
void initialise_image(int ***image, const int grid_size_x, const int grid_size_y)
{
    int i;
    int j;
    *image = (int**)arralloc(sizeof(int), 2, grid_size_y, grid_size_x);

    if ( NULL == *image ) {
        error(1, errno, "Unable to allocate memory for image\n");
    }
    /* initalise results array to black */
    for ( i = 0; i < grid_size_y; i++ ) {
        for ( j = 0; j < grid_size_x; j++ ) {
            (*image)[i][j] = -1;
        }
    }
}

static inline void send_image_pixels(int **image_pixels,
                                     int loc[4],
                                     MPI_Comm comm)
{
    MPI_Send(loc, 4, MPI_INT, 0, 0, comm);
    MPI_Send(&(image_pixels[0][0]),
             (loc[2] - loc[0] + 1)*(loc[3] - loc[1] + 1),
             MPI_INT, 0, 0, comm);
}

static inline int recv_image_pixels(int **image,
                                    int max_iter,
                                    int *load,
                                    MPI_Comm comm)
{
    int loc[4];
    MPI_Status status;
    int i, j, worker;

    int xsize, ysize;
    int **imagesection;

    /* Receive information about pixel location from any sender */
    MPI_Recv(loc, 4, MPI_INT, MPI_ANY_SOURCE, 0, comm, &status);
    /* Then get the pixels from the same send (using status.MPI_SOURCE
     * to find the source of the previous message) */

    worker = status.MPI_SOURCE;
    xsize = loc[2] - loc[0] + 1;
    ysize = loc[3] - loc[1] + 1;

    initialise_image(&imagesection, xsize, ysize);

    MPI_Recv(&(imagesection[0][0]), xsize*ysize,
	       MPI_INT, worker, 0, comm, &status);

    /* Copy data into main image */
    
    for (i=0; i < ysize; i++)
      {
	for (j=0; j < xsize; j++)
	  {
	    /* Hack to encode worker number in values */

	    // map value from [-1, maxiter] to [0, max_iter+1]
	    // then add multiples of (max_iter+2)
	    // add mapped result to load
	    *load = *load+imagesection[i][j]+1;
	    image[loc[1]+i][loc[0]+j] =
	      imagesection[i][j] + 1 + worker*(max_iter+2);
	  }
      }

    free(imagesection);

    /* Return where we got the pixels from, necessary to send a new task out */
    return worker;
}

int **compute_pixels(in_set_fn_t in_set_fn,
                    const int loc[4],
                    const float xmin, const float xmax,
                    const float ymin, const float ymax,
		    const float cre, const float cim,
                    const int grid_size_x, const int grid_size_y,
                    const int max_iter)
{
    int i;
    int j;
    int xstart, xstop, ystart, ystop;
    int xsize, ysize;
    int **imagesection;;
    float x0;
    float y0;

    xstart = loc[0];
    xstop  = loc[2];
    ystart = loc[1];
    ystop  = loc[3];

    xsize = xstop - xstart + 1;
    ysize = ystop - ystart + 1;

    initialise_image(&imagesection, xsize, ysize);

    for ( i = 0; i < ysize; i++ ) {
      for ( j = 0; j < xsize; j++ ) {
            /* calculate coordinates x0,y0 of current pixel */
	x0 = xmin + (xstart+j) * ((xmax - xmin) / grid_size_x);
	y0 = ymin + (ystart+i) * ((ymax - ymin) / grid_size_y);
	imagesection[i][j] = (*in_set_fn)(x0, y0, cre, cim, max_iter);
        }
    }
    return imagesection;
}

/* A task is a square block of pixels, not necessarily contiguous,
 * which can be described by an array of four ints: the start line,
 * the start pixel, the end line, and the end pixel.  The task
 * variable holds the location of the next task to be processed. */

static inline void build_task(int loc[4], int task[2], int grid_size_x, int grid_size_y, int pixels_in_task)
{
  if ( task[0] < grid_size_x && task[1] < grid_size_y )
    {
      loc[0] = task[0];
      loc[1] = task[1];
      loc[2] = loc[0] + pixels_in_task - 1;
      loc[3] = loc[1] + pixels_in_task - 1;;

      task[0] = loc[2]+1;
      task[1] = loc[1];

      if ( loc[2] >= grid_size_x-1 ) {
	loc[2] = grid_size_x - 1;
	task[0] = 0;
	task[1] += pixels_in_task;
      }
      if ( loc[3] >= grid_size_y - 1) {
	loc[3] = grid_size_y - 1;
      }
    }
  else
    {
      loc[0] = EMPTY_TASK;
      loc[1] = EMPTY_TASK;
      loc[2] = EMPTY_TASK;
      loc[3] = EMPTY_TASK;
    }
}

/* Helper to get time as a double */
static inline double getTime(){

  return MPI_Wtime();

}


/* write run configuration */
void writeConfig(const int *nProcs, const int * gridX, const int * gridY,
		 const int *iter, const float * xMin, const float * xMax,
		 const float * yMin, const float * yMax,
		 char function_name, const float *cre, const float *cim,
		 const int *pixels_in_task){

  printf("\n\n--------- CONFIGURATION OF THE TASKFARM RUN ---------\n\n");
  printf("Number of processes:\t\t\t %d\n",*nProcs);
  printf("Image size:\t\t\t %d x %d \n",*gridX, *gridY);
  printf("Task size:\t\t\t %d x %d (pixels)\n",*pixels_in_task,*pixels_in_task);
  printf("Number of iterations:\t\t %d\n",*iter);
  printf("Coordinates in X dimension:\t %f to %f\n",*xMin, *xMax);
  printf("Coordinates in Y dimension:\t %f to %f\n",*yMin, *yMax);
  if(function_name == FRACTAL_FUNCTION)
    {
      printf("Fractal function is:\t\t Mandelbrot set\n");
    }
  else
    {
      printf("Fractal function is:\t\t Julia set\n");
      printf("Value of C is:\t\t\t (%f, %f)\n",*cre, *cim);
    }
}

/* Write out workload information */
void writeInformation(const double elapsedTime, const int nProcs,const int load,const int maxLoad,const int minLoad,const  int numberOfTasks, const int * cpuLoads, const int* cpuMaxLoads,const int* cpuTasks)
{
    int count;
    int nworker = nProcs - 1;
    int aveLoad = (int) ((double) load)/((double) nworker);

    printf("\n-----Workload Summary (number of iterations)---------\n\n");
    printf("Total Number of Workers: %d\n", nworker);
    printf("Total Number of Tasks:   %d\n", numberOfTasks);
    printf("\n");
    printf("Total   Worker Load: %d\n", load);
    printf("Average Worker Load: %d\n", aveLoad);
    printf("Maximum Worker Load: %d\n", maxLoad);
    printf("Minimum Worker Load: %d\n", minLoad);
    printf("\n");
    printf("Time taken by %d workers was %f (secs)\n",nworker, elapsedTime);
    printf("Load Imbalance Factor: %f\n", ((double) maxLoad)/((double) aveLoad));
    printf("\n");
}

void master_loop(int **image,
                 const float xmin,
                 const float xmax,
                 const float ymin,
                 const float ymax,
		 char function_name,
                 const float cre,
                 const float cim,
                 const int grid_size_x,
                 const int grid_size_y,
                 const int max_iter,
                 const int pixels_in_task,
                 MPI_Comm comm)
{
    int current_task[2];
    int size;
    int loc[4];
    int i,j;
    int outstanding_tasks;
    int temp;

    MPI_Comm_size(comm, &size);
   
    /* Workload Gathering */
    int maxBlockLoad=-1, minBlockLoad=-1, currentLoad = 0, totalLoad=0;
    int maxWorkerLoad, minWorkerLoad;
    double startTime, endTime, elapsedTime;
    int cpuLoad[size];
    int cpuTasks[size];
    int cpuMaxLoad[size];
    int totalTasks =0;
    for( temp =0; temp < size; temp++)
      {
        cpuLoad[temp] = 0;
        cpuTasks[temp]=0;
        cpuMaxLoad[temp]=0;
      }
    writeConfig(&size, &grid_size_x, &grid_size_y, &max_iter,  &xmin, &xmax, &ymin, &ymax, function_name, &cre, &cim, &pixels_in_task);
    startTime = getTime();
    /* End Gathering */


    /* Broadcast constant problem data (extent of image, number of
     * pixels, iteration count) */
    {
      float data[6] = {xmin, xmax, ymin, ymax, cre, cim};
        MPI_Bcast(data, 6, MPI_FLOAT, 0, comm);
    }
    {
        int data[3] = {grid_size_x, grid_size_y, max_iter};
        MPI_Bcast(data, 3, MPI_INT, 0, comm);
    }

    outstanding_tasks = 0;
    current_task[0] = 0;  // Start in first line ...
    current_task[1] = 0;  // ... at the first pixel.
    /* Farm out first round of tasks (all workers are idle) */
    for ( i = 1; i < size; i++ ) {
      build_task(loc, current_task, grid_size_x, grid_size_y, pixels_in_task);
        MPI_Send(loc, 4, MPI_INT, i, 0, comm);
        /* Don't have to poll for a response if the task was empty */
        if ( loc[0] != EMPTY_TASK )
            outstanding_tasks++;
    }

    /* Sit waiting for data, receive it, then hand out new task */
    while ( outstanding_tasks ) {
      int worker = recv_image_pixels(image, max_iter, &currentLoad, comm);
      /* Gathering Information for workloads */
      if (maxBlockLoad < 0 || minBlockLoad < 0) {
	maxBlockLoad = currentLoad; minBlockLoad = currentLoad;
      }

      totalLoad = totalLoad + currentLoad;
      if(currentLoad < minBlockLoad) { minBlockLoad = currentLoad;}
      if(currentLoad > maxBlockLoad) { maxBlockLoad = currentLoad;}
      cpuLoad[worker] = cpuLoad[worker] + currentLoad;
      cpuTasks[worker]++;
      if(currentLoad > cpuMaxLoad[worker]) { cpuMaxLoad[worker] = currentLoad;}
      totalTasks++;
      currentLoad = 0;
      /* End Gathering Section */
        --outstanding_tasks;
        build_task(loc, current_task, grid_size_x, grid_size_y, pixels_in_task);
        MPI_Send(loc, 4, MPI_INT, worker, 0, comm);
        if ( loc[0] != EMPTY_TASK ) {
            ++outstanding_tasks;
        }
    }
    
    /* Timing and Load WriteOut */
    endTime = getTime();
    elapsedTime = endTime-startTime;

    maxWorkerLoad = cpuLoad[1];
    minWorkerLoad = cpuLoad[1];
    
    for (i=2; i < size-1; i++) {
      if (cpuLoad[i] > maxWorkerLoad) maxWorkerLoad = cpuLoad[i];
      if (cpuLoad[i] < minWorkerLoad) minWorkerLoad = cpuLoad[i];
    }

    writeInformation(elapsedTime, size, totalLoad, maxWorkerLoad, minWorkerLoad, totalTasks, cpuLoad, cpuMaxLoad, cpuTasks);
}

void worker_loop(in_set_fn_t in_set_fn, MPI_Comm comm)
{
    float xmin, xmax, ymin, ymax, cre, cim;
    int grid_size_x, grid_size_y, max_iter;
    MPI_Status status;
    int loc[4];

    /* Receive constant problem data */
    {
        float data[6];
        MPI_Bcast(data, 6, MPI_FLOAT, 0, comm);
        xmin = data[0];
        xmax = data[1];
        ymin = data[2];
        ymax = data[3];
	cre =  data[4];
	cim =  data[5];
    }
    {
        int data[3];
        MPI_Bcast(data, 3, MPI_INT, 0, comm);
        grid_size_x = data[0];
        grid_size_y = data[1];
        max_iter = data[2];
    }

    do {
        /* Receive task description */
        MPI_Recv(loc, 4, MPI_INT, 0, 0, comm, &status);
        if ( loc[0] == EMPTY_TASK) {
            /* Nothing do, so return */
            return;
        } else {
            /* Compute the image pixels */
            int **data = compute_pixels(in_set_fn, loc,
					xmin, xmax, ymin, ymax, cre, cim,
					grid_size_x, grid_size_y, max_iter);
            send_image_pixels(data, loc, comm);
            free(data);
        }
    } while ( 1 )
        ;
}

void compute_set(in_set_fn_t in_set_fn,
                 int **image,
                 const float xmin,
                 const float xmax,
                 const float ymin,
                 const float ymax,
		 char function_name,
		 const float cre,
		 const float cim,
                 const int grid_size_x,
                 const int grid_size_y,
                 const int max_iter,
                 const int pixels_in_task,
                 MPI_Comm comm)
{
    int rank;
    int size;

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);
    if (size > 1) {

      if ( rank == 0 ) {

        master_loop(image, xmin, xmax, ymin, ymax, function_name, cre, cim,
                    grid_size_x, grid_size_y, max_iter, pixels_in_task,
                    comm);
      } else {
        worker_loop(in_set_fn, comm);
      }
    } else {

      if (rank == 0) printf("\nERROR: need at least two processes for the task farm!\n\n");
    }

}

int main(int argc, char** argv)
{
    int grid_size_x;
    int grid_size_y;
    int max_iter;
    float xmin;
    float xmax;
    float ymin;
    float ymax;
    float cre;
    float cim;
    int pixels_in_task;
    int shading;
    int size, numworker;
    int **image;
    int rank;
    MPI_Comm comm;
    in_set_fn_t fp;
    char function_name;

    MPI_Init(&argc, &argv);

    comm = MPI_COMM_WORLD;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
   
    read_options(argc, argv, &grid_size_x, &grid_size_y, &max_iter,
                 &xmin, &xmax, &ymin, &ymax, &function_name, &cre, &cim,
		 &pixels_in_task, &shading);
   
    if(function_name == FRACTAL_FUNCTION)
      {
         fp = &point_in_mandelbrot_set;
      }
    else
      {
	fp = &point_in_julia_set;
      }

   
 
    if ( rank == 0 ) {
        initialise_image(&image, grid_size_x, grid_size_y);
    }

    compute_set(fp, image,
                xmin, xmax, ymin, ymax, function_name, cre, cim,
                grid_size_x, grid_size_y, max_iter, pixels_in_task,
                comm);

    if ( rank == 0 ) {
      if (shading == SHADING_ON)
	{
	  numworker = size-1;
	}
      else
	{
	  numworker = 1;
	}

        write_ppm("output.ppm", image, grid_size_x, grid_size_y, max_iter, numworker);
        free(image);
    }

    MPI_Finalize();


    return 0;
}



