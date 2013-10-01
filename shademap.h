#ifndef __CLOCK_SUNPOS_H__
#define __CLOCK_SUNPOS_H__

//#include <glib.h>

#include <png.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

typedef struct {
    uint8_t gray;
    uint8_t alpha;
} pixel_t;


typedef struct {
    double lat;
    double lon;
    int size;
    char file[255];
} marker_t;

/* A picture. */
    
typedef struct  {
    pixel_t *pixels;
    int width;
    int height;
} bitmap_t;

void sun_position(time_t unix_time, double *lat, double *lon);
time_t now;

#endif
