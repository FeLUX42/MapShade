/*
 * Copyright (C) 2013 Felix Passenberg
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 2 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
 * 02110-1301, USA.
 *
 * Autor:
 * 		Felix Passenberg <fc.passenberg@gmail.com>
 *
 * Vorlage:
 *     Jonathan Blandford <jrb@redhat.com>
 *     Matthias Clasen <mclasen@redhat.com>
 
 Version 0.9: basic use
 Version 1.0:
 	(x) accept size as arguments 
 	(x)	print sunspot
 	(x)	print location markers
 	(x)	selectable maps
 	(x)	select output format
 Version 1.1:
 	(x) resize maps
 	(x) auto resize shade
 	(x) size for default marker
 	(x) size for default sun
 Version 1.2:
 	(x) bugs (sun pos in winter, wrong pos when not using -v)
 	(x) specified time
 
 ToDo:
 	autofit maps
 	use imagemagick api 
 
 */
 // examele arguments 
 //./a.out -v -l 49.449167°,8.672222° 3 -r -d 5400 -s -sf sun_clipart.gif -n dnb_land_ocean_ice.2012.3600x1800.jpg  -m world.200404.3x5400x2700.jpg 

 // compiled with  gcc -Wall shademap.c -g -O2 -lm -lpng -o daylightmap.out

#include <time.h>
#include <png.h>
#include <stdlib.h>
//#include <gtk.h>
#include <stdio.h>
#include <math.h>
#include "shademap.h"


/* Calculated with the methods and figures from "Practical Astronomy With Your
 * Calculator, version 3" by Peter Duffet-Smith.
 */
/* Table 6.  Details of the Sun's apparent matecorba at epoch 1990.0 */

#define EPOCH          2447891.5  /* days */    /* epoch 1990 */
#define UNIX_EPOCH     2440586.5  /* days */    /* epoch 1970 */
#define EPSILON_G      279.403303 /* degrees */ /* ecliptic longitude at epoch 1990.0 */
#define MU_G           282.768422 /* degrees */ /* ecliptic longitude at perigee */
#define ECCENTRICITY   0.016713                 /* eccentricity of matecorba */
#define R_0            149598500  /* km */      /* semi-major axis */
#define THETA_0        0.533128   /* degrees */ /* angular diameter at r = r_0 */
#define MEAN_OBLIQUITY 23.440592  /* degrees */ /* mean obliquity of earth's axis at epoch 1990.0 */

#define NORMALIZE(x) \
  while (x>360) x-=360; while (x<0) x+= 360;

#define DEG_TO_RADS(x) \
  (x * M_PI/180.0)

#define RADS_TO_DEG(x) \
  (x * 180.0/M_PI)

/* Calculate number of days since 4713BC.
 */
static double
unix_time_to_julian_date (int unix_time)
{
  return UNIX_EPOCH + (double) unix_time / (60 * 60 * 24);
}

/* Finds an iterative solution for [ E - e sin (E) = M ] for values of e less
   than 0.1.  Page 90  */

#define ERROR_ACCURACY 1e-6 /* radians */

static double
solve_keplers_equation (double e,
			double M)
{
  double d, E;

  /* start with an initial estimate */
  E = M;
  
  d = E - e * sin (E) - M;

  while (fabs(d) > ERROR_ACCURACY)
    {
      E = E - (d / (1 - e * cos (E)));
      d = E - e * sin (E) - M;
    }

  return E;
}

static void
getsize( char * file, int * x, int * y){
	char imagemagick[255];
	sprintf(imagemagick, "identify %s | sed -e \'s| |\\n|g\' | head -n 3 | tail -n 1", file);
	FILE * res;
	res = popen ( imagemagick,"r");
	
	if (res!=NULL){
		char arres[130];
		fgets( arres, sizeof arres, res);
		
		char * pch = strtok(arres, "x"); // grap position of x
        int dim[] = {0,0};
        int k = 0;
        while (pch != NULL && k <= 1)
		{
			dim [k] = atoi(pch);
			pch = strtok (NULL, "x");
			k++;
		}
       *x = dim[0];
       *y = dim[1];
       
	}else{
		*x = 0;
		*y = 0;
	}
	pclose (res);
}


	//convert lon lat to x,y
static void
flatten( double lon , double lat, int width, int height, int * x, int * y)
{
//	lat ^
//	lon >

	*x = (int) (180 + lon);
	*y = (int) (90 - lat);
	
	if (*x > 360)  *x -= 360;
	
	*x = *x * width / 360;
	*y = *y * height / 180;
}


  /* convert the ecliptic longitude to right ascension and declination.  Section 27.  */
static void
ecliptic_to_equatorial (double  lambda,
			double  beta,
			double *ra,
			double *dec)
{
  double cos_mo;
  double sin_mo;


  sin_mo = sin (DEG_TO_RADS (MEAN_OBLIQUITY));
  cos_mo = cos (DEG_TO_RADS (MEAN_OBLIQUITY));

  *ra = atan2 (sin (lambda) * cos_mo - tan (beta) * sin_mo, cos (lambda));
  *dec = asin (sin (beta) * cos_mo + cos (beta) * sin_mo * sin (lambda));
}

/* calculate GST.  Section 12  */
static double
greenwich_sidereal_time (double unix_time)
{
  double u, JD, T, T0, UT;

  u = fmod (unix_time, 24 * 60 * 60);
  JD = unix_time_to_julian_date (unix_time - u);
  T = (JD - 2451545) / 36525;
  T0 = 6.697374558 + (2400.051336 * T) + (0.000025862 * T * T);
  T0 = fmod (T0, 24);
  UT = u / (60 * 60);
  T0 = T0 + UT * 1.002737909;
  T0 = fmod (T0, 24);

  return T0;
}

/* Calculate the position of the sun at a given time.  pages 89-91 */
void
sun_position (time_t unix_time, double *lat, double *lon)
{
  double jd, D, N, M, E, x, v, lambda;
  double ra, dec;
  jd = unix_time_to_julian_date (unix_time);

  /* Calculate number of days since the epoch */
  D = jd - EPOCH;

  N = D*360/365.242191;

  /* normalize to 0 - 360 degrees */
  NORMALIZE (N);

  /* Step 4: */
  M = N + EPSILON_G - MU_G;
  NORMALIZE (M);

  /* Step 5: convert to radians */
  M = DEG_TO_RADS (M);

  /* Step 6: */
  E = solve_keplers_equation (ECCENTRICITY, M);

  /* Step 7: */
  x = sqrt ((1 + ECCENTRICITY)/(1 - ECCENTRICITY)) * tan (E/2);

  /* Step 8, 9 */
  v = 2 * RADS_TO_DEG (atan (x));
  NORMALIZE (v);

  /* Step 10 */
  lambda = v + MU_G;
  NORMALIZE (lambda);

  /* convert the ecliptic longitude to right ascension and declination */
  ecliptic_to_equatorial (DEG_TO_RADS (lambda), 0.0, &ra, &dec);

  ra = ra - (M_PI/12) * greenwich_sidereal_time (unix_time);
  ra = RADS_TO_DEG (ra);
  dec = RADS_TO_DEG (dec);
  NORMALIZE (ra);
  NORMALIZE (dec);

  *lat = dec;
  *lon = ra;
}

static void
clock_map_compute_vector (double lat, double lon, double *vec)
{
        double lat_rad, lon_rad;
        lat_rad = lat * (M_PI/180.0);
        lon_rad = lon * (M_PI/180.0);

        vec[0] = sin(lon_rad) * cos(lat_rad);
        vec[1] = sin(lat_rad);
        vec[2] = cos(lon_rad) * cos(lat_rad);
}

static uint8_t
clock_map_is_sunlit (double pos_lat, double pos_long,
                         double sun_lat, double sun_long)
{
        double pos_vec[3];
        double sun_vec[3];
        double dot;

        /* twilight */
        double epsilon = 0.05;

        clock_map_compute_vector (pos_lat, pos_long, pos_vec);
        clock_map_compute_vector (sun_lat, sun_long, sun_vec);

        /* compute the dot product of the two */
        dot = pos_vec[0]*sun_vec[0] + pos_vec[1]*sun_vec[1]
                + pos_vec[2]*sun_vec[2];

        if (dot > epsilon) {
                return 0x00;
        }

        if (dot < -epsilon) {
                return 0xFF;
        }

        return (uint8_t)(-128 * ((dot / epsilon) - 1));
}

static pixel_t * pixel_at (bitmap_t * bitmap, int x, int y)
{
    return bitmap->pixels + bitmap->width * y + x;
}

static void
clock_map_render_shadow_pixbuf (bitmap_t *pixbuf)
{
        int x, y;
        int height, width;
        double sun_lat, sun_lon;
      //  time_t now = time (NULL);

        width = pixbuf->width;
        height = pixbuf->height;

        sun_position (now, &sun_lat, &sun_lon);
        
        printf("shade is: x=%i y=%i\n", pixbuf->width, pixbuf->height);

        for (y = 0; y < height; y++) {
                double lat = (height / 2.0 - y) / (height / 2.0) * 90.0;

                for (x = 0; x < width; x++) {
                        
                        double lon = (x - width / 2.0) / (width / 2.0) * 180.0;
						uint8_t shade;
		                shade = clock_map_is_sunlit (lat, lon, sun_lat, sun_lon);
                        
                        
                       	pixel_t * pixel = pixel_at ( pixbuf, x, y);
                       	pixel->gray = 0x00;
                       	pixel->alpha = shade;
          			//	pixel->red = shade;
            		//	pixel->green = shade;
            		//	pixel->blue = shade;
                }
        }
}

/* Given "bitmap", this returns the pixel of bitmap at the point 
   ("x", "y"). */


    
/* Write "bitmap" to a PNG file specified by "path"; returns 0 on
   success, non-zero on error. */

static int save_png_to_file (bitmap_t *bitmap, const char *path)
{
    FILE * fp;
    png_structp png_ptr = NULL;
    png_infop info_ptr = NULL;
    size_t x, y;
    png_byte ** row_pointers = NULL;
    /* "status" contains the return value of this function. At first
       it is set to a value which means 'failure'. When the routine
       has finished its work, it is set to a value which means
       'success'. */
    int status = -1;
    /* The following number is set by trial and error only. I cannot
       see where it it is documented in the libpng manual.
    */
    int channels = 2;
	int depth = 8;
	
    fp = fopen (path, "wb");
    if (! fp) {
        goto fopen_failed;
    }

    png_ptr = png_create_write_struct (PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    if (png_ptr == NULL) {
        goto png_create_write_struct_failed;
    }
    
    info_ptr = png_create_info_struct (png_ptr);
    if (info_ptr == NULL) {
        goto png_create_info_struct_failed;
    }
    
    /* Set up error handling. */

    if (setjmp (png_jmpbuf (png_ptr))) {
        goto png_failure;
    }
    
    /* Set image attributes. */

    png_set_IHDR (png_ptr,
                  info_ptr,
                  bitmap->width,
                  bitmap->height,
                  depth,
                  PNG_COLOR_TYPE_GRAY_ALPHA,
                  PNG_INTERLACE_NONE,
                  PNG_COMPRESSION_TYPE_DEFAULT,
                  PNG_FILTER_TYPE_DEFAULT);
    
    /* Initialize rows of PNG. */

    row_pointers = png_malloc (png_ptr, bitmap->height * sizeof (png_byte *));
    for (y = 0; y < bitmap->height; ++y) {
        png_byte *row = 
            png_malloc (png_ptr, sizeof (uint8_t) * bitmap->width * channels);
        row_pointers[y] = row;
        for (x = 0; x < bitmap->width; ++x) {
            pixel_t * pixel = pixel_at (bitmap, x, y);
            *row++ = pixel->gray;
            *row++ = pixel->alpha;
         //   *row++ = pixel->red;
         //   *row++ = pixel->green;
         //   *row++ = pixel->blue;
        }
    }
    
    /* Write the image data to "fp". */

    png_init_io (png_ptr, fp);
    png_set_rows (png_ptr, info_ptr, row_pointers);
    png_write_png (png_ptr, info_ptr, PNG_TRANSFORM_IDENTITY, NULL);

    /* The routine has successfully written the file, so we set
       "status" to a value which indicates success. */

    status = 0;
    
    for (y = 0; y < bitmap->height; y++) {
        png_free (png_ptr, row_pointers[y]);
    }
    png_free (png_ptr, row_pointers);
    
	 png_failure:
	 png_create_info_struct_failed:
		png_destroy_write_struct (&png_ptr, &info_ptr);
	 png_create_write_struct_failed:
		fclose (fp);
	 fopen_failed:
		return status;
}

int main (int argc, char *argv[])
{
  //GTimeVal timeval;
  double lat, lon;

 // gtk_init (&argc, &argv);

//  g_get_current_time (&timeval);

	
	bitmap_t shade;
	shade.width = 720;
    shade.height = 360;
    int i = 0;
    int sun = 0;
    int timevar = time(NULL);
    int keep = 0;
    int autoresize = 0;
    int verbose = 0;
    char daymap[255] = " ";
    char nightmap[255] = " ";
    char output[255] = " ";
    char sunfile[255] = " ";
    
    marker_t * markers;
    markers = NULL;
    int markercount = 0;
    
    
    for (i = 1; i < argc; i++)  /* Skip argv[0] (program name). */
    {
		if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0 )  /* Process optional arguments. */
		    {
		    	/* Print usage statement and exit (see below). */
		    	printf(" -d \tset dimensons to <width>x<height>\n");
		    	printf(" -df\tforce dimensons to <width>x<height>\n");
		    	printf(" -s \tprint the sun\n \tyou may parse the size of the marker");
		    	printf(" -sd \tdimensions for default sun\n");
		    	printf(" -sf \tfilename of the sun marker\n");
		    	printf(" -l \tcan be used more than once, add a location marker to <lat>,<lon>\n \t enter <lat> and <lon> as float; <lat> as N and <lon> as E\n \t use \"<lon>,<lat>\" for default marker\n \t use \"<lat>,<lon> <filename>\" for custom marker \n \t use\"<lat>,<lon> <size>\" for marker size \n");
		    	printf(" -m \tdaylight map to use, filename\n");
		    	printf(" -n \tnight map to use, filename\n");
		    	printf(" -o \toutput file\n");
		    	printf(" -k \tkeep automatically generated temporary files\n");
		    	printf(" -r \tauto resize maps to given dimension\n");
		    	printf(" -t \tset time for which the shade should be calculated\n");
		    	printf(" -h \n --help\topen this help option\n");
		    	return 0;
		    }
		    
		if ( strcmp(argv[i], "-v") == 0) verbose ++;
		
		if (strcmp(argv[i], "-t") == 0) // set time to calculate
		{

			if (i + 1 <= argc - 1)  /* There are enough arguments in argv. */
			{
				i++;
				if(verbose)	printf("setting time from %i to", (int)timevar);
				timevar = atoi(argv[i]);
				if(verbose) printf(" %i.\n", (int)timevar);
			}
			else
			{
				printf("not enough arguments passed to %s \n", argv[i]);
				return 1;
			}
		}
		
		if (strcmp(argv[i], "-r") == 0) // resize map automatically
		{
			if(verbose)	printf("resizing maps automatically\n");
			autoresize = 1;
		}
		
		if (strcmp(argv[i], "-s") == 0) // display sunspot and optional size
		{
			if(verbose)	printf("also printing sunspot\n");
			sun = 10;
			if (i+1 <= argc -1){
				if( strchr(argv[i+1], '-') == NULL){
					i++;
					sun = atoi(argv[i]);
					if( verbose) printf("set dimensions of sun to %i\n", sun );
				}
			}
		}
		
		if (strcmp(argv[i], "-k") == 0) // keep temp files
		{
			if(verbose) printf("keeping temporary files\n");
			keep = 1;
		}
		
		if (strcmp(argv[i], "-l") == 0) // add location marker
		{

            if (i + 1 <= argc - 1)  /* There are enough arguments in argv. */
            {
            	markercount ++;
                i++;
                char * pch = strtok(argv[i], ","); // grap position of x
                double dim[] = {0,0};
                int k = 0;
                while (pch != NULL && k <= 1) //read the numbers 
				{
					dim [k] = atof(pch);
					pch = strtok (NULL, "x");
					k++;
				}
				
				marker_t * tmpmarker;			// resize the old array

				tmpmarker = (marker_t*) realloc(markers, markercount * (sizeof(marker_t)));
				if(tmpmarker == NULL){
					free(markers);
					printf("error allocating new space for marker %i\n", markercount);
					return 2;
				}

				markers = tmpmarker;
				markers[markercount -1].lat = dim[0];
				markers[markercount -1].lon = dim[1];
				strcpy(markers[markercount -1].file, " ");
				markers[markercount -1].size = 3;
            	if (i+1 <= argc -1){
					if( strchr(argv[i+1], '.') != NULL){
						i++;
						strcpy (markers[markercount -1].file, argv[i]);
						if( verbose) printf("added marker to pos lat=%f lon=%f using picture file %s\n", 
												dim[0], dim[1], markers[markercount -1].file );
					}
					else{
						
						if( strchr(argv[i+1], '-') == NULL){
							i++;
							markers[markercount -1].size = atoi(argv[i]);
							if( verbose) printf("added marker to pos lat=%f lon=%f and size=%i\n", 
																			dim[0], dim[1], markers[markercount -1].size);
						}
						else if( verbose) printf("added marker to pos lat=%f lon=%f\n", dim[0], dim[1]);
            		}
				}
				else{
					if( verbose) printf("added marker to pos lat=%f lon=%f\n", dim[0], dim[1]);
            	}
            }
            else
            {
                /* Print usage statement and exit (see below). */
                printf("not enough arguments passed to %s \n", argv[i]);
                return 1;
            }
		}
		
		if (strcmp(argv[i], "-sf") == 0)  //select sunfile
		{
			if (i + 1 <= argc - 1)  /* There are enough arguments in argv. */
			{
				i++;
				strcpy( sunfile, argv[i]);
				if (verbose) printf("sunfile set to \"%s\"\n", sunfile);
			}
			else
			{
				printf("not enough arguments passed to %s \n", argv[i]);
				return 1;
			}
		}
		
		if (strcmp(argv[i], "-m") == 0)  // select daymap
		{

			if (i + 1 <= argc - 1)  /* There are enough arguments in argv. */
			{
				i++;
				strcpy( daymap, argv[i]);
				if (verbose) printf("daymap set to \"%s\"\n", daymap);
			}
			else
			{
				printf("not enough arguments passed to %s \n", argv[i]);
				return 1;
			}
		}
		
		if (strcmp(argv[i], "-n") == 0) // select nightmap
		{
			if (i + 1 <= argc - 1)  /* There are enough arguments in argv. */
			{
				i++;
				strcpy(nightmap, argv[i]);
				if (verbose) printf("nightmap set to \"%s\"\n", nightmap);
			}
			else
			{
				printf("not enough arguments passed to %s \n", argv[i]);
				return 1;
			}
		}
		
		if (strcmp(argv[i], "-o") == 0)  // select output
		{
			if (i + 1 <= argc - 1)  /* There are enough arguments in argv. */
			{
				i++;
				strcpy(output, argv[i]);
				if (verbose) printf("output set to \"%s\"\n", output);
			}
			else
			{
				printf("not enough arguments passed to %s \n", argv[i]);
				return 1;
			}
			
		}
		
		if (strcmp(argv[i], "-d") == 0)  // set dimensions
        {
            if (i + 1 <= argc - 1)  /* There are enough arguments in argv. */
            {
                i++;
                char * pch = strtok(argv[i], "x"); // grap position of x
                int dim[] = {0,0};
                int k = 0;
                while (pch != NULL && k <= 1)
				{
					dim [k] = atoi(pch);
					pch = strtok (NULL, "x");
					k++;
				}
               shade.width = dim[0];
               shade.height = dim[0]/2;
               if(shade.height != dim[0]){
               		printf("WARNING: aspect ratio hold, if you want to chang it use -df\n");
               }
              if( verbose) printf("dimensions of the shade x=%i y=%i\n", shade.width, shade.height);
            }
            else
            {
                /* Print usage statement and exit (see below). */
                printf("not enough arguments passed to %s \n", argv[i]);
                return 1;
            }
            if(verbose)	printf("aslso resizing maps\n");
			autoresize = 1;
        }
		
        if (strcmp(argv[i], "-df") == 0)  // set dimensions force
        {
            if (i + 1 <= argc - 1)  /* There are enough arguments in argv. */
            {
                i++;
                char * pch = strtok(argv[i], "x"); // grap position of x
                int dim[] = {0,0};
                int k = 0;
                while (pch != NULL && k <= 1)
				{
					dim [k] = atoi(pch);
					pch = strtok (NULL, "x");
					k++;
				}
               shade.width = dim[0];
               shade.height = dim[1];
               
              	if( verbose) printf("dimensions of the shade x=%i y=%i\n", dim[0], dim[1]);
            }
            else
            {
                /* Print usage statement and exit (see below). */
                printf("not enough arguments passed to %s \n", argv[i]);
                return 1;
            }
        }
        
        else
        {
            /* Process non-optional arguments here. */
        }
    }
   		// -------------------------------------------------------------------ende der komandozeleilen argumente ---------

	
	char imagemagick[255];// buffer for many obperations related to imagemagick
	
	if( strcmp(daymap, " ") == 0 ){ // set maps if not set
		strcpy(daymap, "daymap.jpg");
	}
	if( strcmp(nightmap, " ") == 0){
		strcpy(nightmap , "nightmap.jpg");
	}
	if( strcmp(output, " ") == 0){
		strcpy(output, "final.png");
	}
	
	if(!autoresize){				// set shade to fit daymap
		int twidth = 720, theight = twidth /2 ;
		getsize( daymap, &twidth, &theight );	// get image width and height
		if(twidth && theight){		// if values non 0
			shade.width = twidth;
			shade.height = theight;
			autoresize = 1;
			if (verbose) printf("shade rescaled to: x=%i y=%i\n", shade.width, shade.height);
		}
	}
	
	if(shade.width == 0 || shade.height == 0) {
    	printf( "some error occured, one or more dimensions equal 0\n");
    	return 1;
    }    
    shade.pixels = calloc (sizeof (pixel_t), shade.width * shade.height); // allocate space for shade
    
    now = timevar;
   // printf("now %i\n",now);
	sun_position (now, &lat, &lon);
	if(lat > 90) {lat -= 360;}
//	if(lat > 90) lat -= 180;
//	if(lon > 180) lon -= 180;
    if(verbose) { // some infos about now
		printf("current time :%d\nposition of the sun: %f°N %f°E\ncomputing shademap\n", (int)now, lat, lon); 
    }
    
	printf("computinge shade\n");
	clock_map_render_shadow_pixbuf (&shade);
	
	
	printf("saving shademap\n");
	save_png_to_file (& shade, "shade.png");


	printf("generaing whole map\n");	
	if( autoresize ){		//--------------------------------- resize day and nightmap to fit the shade dimensions
		int twidth = 720, theight = twidth /2 ;
		getsize( daymap, &twidth, &theight );	// get image width and height
		if (verbose) printf("dimensions of the daymap: x=%i y=%i\n", twidth, theight);
		if(twidth && theight) 	//--------------------------------- check if daymap has usable size
			if (theight != shade.height || twidth != shade.width){	//-------if not already right size
				if (verbose) printf("resizing daymap\n");
				char daycpy[255]; 
				strcpy( daycpy, daymap);
				char * daymapname = strtok(daycpy, "."); // grap position of x;
				char * appendix = strtok(NULL, "."); // grap position of x;
				//printf( "convert %s -resize %ix%i %s-r.%s \n", daymap, shade.width, shade.height, daymapname, appendix);
				sprintf(imagemagick, "convert %s -resize %ix%i! %s-r.%s ", daymap ,shade.width, shade.height, daymapname, appendix);
				if (verbose) printf("new dimensions of the daymap: x=%i y=%i\n", shade.width, shade.height);
				system( imagemagick );
				autoresize += 0b10;	// for auto clean
				sprintf(daymap, "%s-r.%s", daymapname, appendix);
			}
		
		getsize(nightmap, &twidth, &theight);
		if (verbose) printf("dimensions of the nightmap: x=%i y=%i\n", twidth, theight);
		if(twidth && theight)	//--------------------------------- check if nightmap has usable size
			if (twidth != shade.width || theight != shade.height){	//------ if not already in right size
				if (verbose) printf("resizing nightmap\n");
				char nightcpy[255]; 
				strcpy( nightcpy, nightmap);
				char * nightmapname = strtok(nightcpy, "."); // grap position of x;
				char * appendix = strtok(NULL, "."); // grap position of x;
				sprintf(imagemagick, "convert %s -resize %ix%i! %s-r.%s ", nightmap, shade.width, shade.height, nightmapname,appendix);
				if (verbose) printf("new dimensions of the nightmap: x=%i y=%i\n", shade.width, shade.height);
				system( imagemagick );
				autoresize += 0b100;	// for autoclean
				sprintf(nightmap,"%s-r.%s", nightmapname, appendix);
			}
	}
	
	
	if (verbose) printf("compositing maps\n");
	sprintf(imagemagick, "composite -compose Dst_In -alpha Set shade.png %s  night-shade.png", nightmap);
	system( imagemagick );
	
	sprintf(imagemagick, "composite -compose Atop night-shade.png %s day-night-shade.png", daymap);
	system( imagemagick );
	

	sprintf(imagemagick, "composite -dissolve 50 -alpha Set shade.png day-night-shade.png %s", output);
	system( imagemagick );
	

	
	
	if( sun >= 1){ 	// add sun spot
		int sposx = 0, sposy = 0;
		flatten(lon, lat, shade.width, shade.height, &sposx, &sposy);
		if( strcmp( sunfile , " ") != 0){
			int swidth = 10, sheight = 10;
			getsize( sunfile, &swidth, &sheight);
			if (verbose) printf("dimensions of the sun: x=%i y=%i\t\tyeah its really this small\n", swidth, sheight);
		
			sprintf(imagemagick, "composite -compose Atop %s -geometry +%i+%i %s %s", 
									sunfile ,sposx-swidth/2, sposy-sheight/2, output, output);
									
			system(imagemagick);
		}else{
			sprintf(imagemagick, " convert -size %ix%i xc:none -fill yellow -draw \'circle %f,%f 0,%f\'  sun.png", 
							sun, sun, (sun -1)/2., (sun -1)/2., (sun -1)/2.);
			system(imagemagick);
			
			sprintf(imagemagick, "composite -compose Atop sun.png -geometry +%f+%f %s %s", sposx-sun/2., sposy-sun/2., output, output);
			system(imagemagick);
			if (!keep) system(" rm sun.png"); // deleting here is easier
		}
	}
	
	if( markercount >= 1){
		int c;
		for(c = 0; c < markercount; c++){
			int mposx = 0, mposy = 0;
			flatten(markers[c].lon, markers[c].lat, shade.width, shade.height, &mposx, &mposy);
			
			if( strcmp( markers[c].file , " ") != 0){
				int mwidth = 3, mheight = 3;
				getsize( markers[c].file, &mwidth, &mheight);				
				if (verbose) printf("dimensions of the markerfile: x=%i y=%i\n", mwidth, mheight);
				
				sprintf(imagemagick, "composite -compose Atop %s -geometry +%i+%i %s %s", 
										markers[c].file ,mposx-mwidth/2, mposy-mheight/2, output, output);
				
				system(imagemagick);
			}else{
				sprintf(imagemagick, " convert -size %ix%i xc:none -fill red -draw \'circle %f,%f 0,%f\'  marker.png", 
							markers[c].size, markers[c].size, (markers[c].size-1)/2., (markers[c].size-1)/2., (markers[c].size-1)/2.);
				system(imagemagick);
				sprintf(imagemagick, "composite -compose Atop marker.png -geometry +%f+%f %s %s", mposx-markers[c].size/2., mposy-markers[c].size/2., output, output);
			
				system(imagemagick);
				if (!keep) system(" rm marker.png"); // deleting here is easier
			}
		}
	}
	
	if( !keep){	// if not wanted delete tem files
		system(" rm shade.png night-shade.png day-night-shade.png");

		if(autoresize & 0b100){
			sprintf(imagemagick, "rm %s", nightmap); // delete resized maps
			system(imagemagick);
			autoresize -= 0b100;
		}
		if(autoresize & 0b10){
			sprintf(imagemagick, "rm %s", daymap); // delete resized maps
			system(imagemagick);
			autoresize -= 0b10;
		}
	}
	free(markers); // cleanup
	printf("final map \"%s\" created\n", output);
	
  return 0;
}


// compiled with gcc -Wall shademap.c -g -O2 -lm -lpng -o daylightmap-1.out

// examele arguments 
 //./a.out -v -l 49.449167°,8.672222° 3 -r -d 5400 -s -sf sun_clipart.gif -n dnb_land_ocean_ice.2012.3600x1800.jpg  -m world.200404.3x5400x2700.jpg 

