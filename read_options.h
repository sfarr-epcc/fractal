#ifndef READ_OPTIONS_H
#define READ_OPTIONS_H

void read_options(int, char**, int*, int*, int *,
                  float*, float*, float*, float*, char*, float*, float*,
		  int*, int*);


/* default coordinates for calculation */
#define XMIN (-2.0)
#define XMAX (1.0)
#define YMIN (-1.5)
#define YMAX (1.5)

/* size of entire grid in X dimension (Y dimension is scaled
 * appropriately) */
#define GRIDSIZE_X 768

/* default number of iterations */
#define ITERATIONS 5000

/* default fractal function */
#define FRACTAL_FUNCTION 'M'

/* default default Julia set parameters */
#define CIM 0.01;
#define CRE 0.285;

/* default task size */
#define DEFAULT_TASK 192

/* shading options */
#define SHADING_OFF 0
#define SHADING_ON 1

#endif
