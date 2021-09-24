# Fractal code for training purposes. 

Makes an image of the Mandelbrot set or the Julia set  and demonstrates 
parallel load balancing in a task farm.

To compile use the provided makefile

	make

Note that this code contains some Linux specific C header files so will only compile on Linux.

The program must be run using mpi with at least 2 processes, e.g.

	mpirun -np 2 fractal

or

	srun -n 2 fractal

To see all options use the help flag

	fractal -h

Which gives the following output
~~~
Usage:
fractal [-SixXyYfth]
   -S NPIXEL    Set number of pixels in X dimension of image
                Y dimension is scaled to ensure pixels are square
   -i ITS       Set max number of iterations for a point to be inside
   -x XMIN      Set xmin coordinate
   -X XMAX      Set xmax coordinate
   -y YMIN      Set ymin coordinate
   -Y YMAX      Set ymax coordinate
   -f FRACTAL_FUNCTION      Set fractal function (M for Mandelbrot, J for Julia)
   -c CRE       Set real part of C = cre + i cim for Julia set
   -C CIM       Set imag part of C = cre + i cim for Julia set
   -Y YMAX      Set ymax coordinate
   -t DEFAULT_SIZE      Set task size
   -n           Do not shade pixels according to the worker process which computed them.
   -h           Show this help
~~~


