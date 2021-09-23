#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "read_options.h"

void read_options(int argc, char**argv,
                  int* gridSizeX,
                  int* gridSizeY,
                  int* iter,
                  float* xMin,
                  float* xMax,
                  float* yMin,
                  float* yMax,
                  char* function,
		  float *cre,
		  float* cim,
		  int* taskSize,
                  int* shading)
{

    int count;
    int tmp;
    float ftmp;

    *gridSizeX = GRIDSIZE_X;
    *xMin = XMIN;
    *xMax = XMAX;
    *yMin = YMIN;
    *yMax = YMAX;
    *iter = ITERATIONS;
    *function = FRACTAL_FUNCTION;
    *cim = CIM;
    *cre = CRE;
    *taskSize = DEFAULT_TASK;
    *shading = SHADING_ON;
    /* loop over command line arguments */
    while ( (count = getopt(argc, argv, "s:S:i:x:X:y:Y:f:c:C:t:nh")) != -1 ) {
        switch (count) {

        case 'S':
            /* option -S: grid size in X dimension */
            tmp = atoi(optarg);
            if ( tmp >= 4 ) *gridSizeX = tmp;
            else fprintf(stderr, "Using default grid size.\n");
            break;
        case 'i':
            /* option -i: number of iterations for the calculation */
            tmp = atoi(optarg);
            if ( tmp > 0 ) *iter = tmp;
            else fprintf(stderr, "Using default number of iterations.\n");
            break;

        case 'x':
            /* option -x:
               for Mandelbrot set: xmin coordinate of the image */
            ftmp = atof(optarg);
            *xMin = ftmp;
            break;

        case 'X':
            /* option -X: xmax coordinate of the image */
            ftmp = atof(optarg);
            *xMax = ftmp;
            break;

        case 'y':
            /* option -y:
               for Mandelbrot set: ymin coordinate of the image
            */
            ftmp = atof(optarg);
            *yMin = ftmp;
            break;

        case 'Y':
            /* option -Y: ymax coordinate of the image */
            ftmp = atof(optarg);
            *yMax = ftmp;
            break;
        
        case 'f':
	  /* option -f: fractal function */
          *function = optarg[0];
	  break;

        case 'c':
            /* option -c:
               for Julia set: value of cre */
            ftmp = atof(optarg);
            *cre = ftmp;
            break;

        case 'C':
            /* option -C:
               for Julia set: value of cim */
            ftmp = atof(optarg);
            *cim = ftmp;
            break;

	case 't':
	  tmp = atoi(optarg);
          *taskSize = tmp;
          break;
        case 'n':
          *shading = SHADING_OFF;
          break;
        case '?':
        case 'h':
            printf("Usage:\n");
            printf("fractal_set [-SixXyYfth]\n");
            printf("   -S NPIXEL    Set number of pixels in X dimension of image\n");
            printf("                Y dimension is scaled to ensure pixels are square\n");
            printf("   -i ITS       Set max number of iterations for a point to be inside\n");
            printf("   -x XMIN      Set xmin coordinate\n");
            printf("   -X XMAX      Set xmax coordinate\n");
            printf("   -y YMIN      Set ymin coordinate\n");
            printf("   -Y YMAX      Set ymax coordinate\n");
	    printf("   -f FRACTAL_FUNCTION      Set fractal function (M for Mandelbrot, J for Julia)\n");
            printf("   -c CRE       Set real part of C = cre + i cim for Julia set\n");
            printf("   -C CIM       Set imag part of C = cre + i cim for Julia set\n");
            printf("   -Y YMAX      Set ymax coordinate\n");
	    printf("   -t DEFAULT_SIZE      Set task size\n");
	    printf("   -n           Do not shade pixels according to the worker process which computed them.\n");
	    printf("   -h           Show this help\n");
            exit(-1);
        }


    }


    /* check that xMin is smaller than xMax, if not, reassign the default values */
    if ( *xMin >= *xMax ) {
        fprintf(stderr, "WARNING: xMin has to be smaller than xMax.\n Reassign default values.\n");
        *xMin = XMIN;
        *xMax = XMAX;
    }

    /* check that yMin is smaller than yMax, if not, reassign the default values */
    if ( *yMin >= *yMax ) {
        fprintf(stderr, "WARNING: yMin has to be smaller than yMax.\n Reassign default values.\n");
        *yMin = YMIN;
        *yMax = YMAX;
    }

    *gridSizeY = (int)(*gridSizeX * (*yMax - *yMin)/(*xMax - *xMin));
}
