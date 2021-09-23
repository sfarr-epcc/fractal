#ifndef WRITE_PPM_H
#define WRITE_PPM_H
#define MAX_COLOUR_VALS 255

enum { R = 0, G, B };

void write_ppm(const char *, int **, int, int, int, int);

#endif
