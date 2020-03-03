//Sanket Manik Salunke  1001764897

#include "bitmap.h"

#include <getopt.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <errno.h>
#include <string.h>
#include<pthread.h>
#include<time.h>

int iteration_to_color( int i, int max );
int iterations_at_point( double x, double y, int max );

//void compute_image( struct bitmap *bm, double xmin, double xmax, double ymin, double ymax, int max , int nthread);            //argument added for nthreads
void *compute_image(void *p);

void show_help()
{
	printf("Use: mandel [options]\n");
	printf("Where options are:\n");
	printf("-m <max>    The maximum number of iterations per point. (default=1000)\n");
	printf("-x <coord>  X coordinate of image center point. (default=0)\n");
	printf("-y <coord>  Y coordinate of image center point. (default=0)\n");
	printf("-s <scale>  Scale of the image in Mandlebrot coordinates. (default=4)\n");
	printf("-W <pixels> Width of the image in pixels. (default=500)\n");
	printf("-H <pixels> Height of the image in pixels. (default=500)\n");
	printf("-o <file>   Set output file. (default=mandel.bmp)\n");
	printf("-h          Show this help text.\n");
	printf("\nSome examples are:\n");
	printf("mandel -x -0.5 -y -0.5 -s 0.2\n");
	printf("mandel -x -.38 -y -.665 -s .05 -m 100\n");
	printf("mandel -x 0.286932 -y 0.014287 -s .0005 -m 1000\n\n");
}

    
	struct params       //struct is created to pass the values
	{
		struct bitmap *bm;
		double xmin;
		double xmax;
		double ymin;
		double ymax; 
		int max;
		int start_image;
		int end_image;
		int nthread;   //nthtread default value
	};
	
	//struct params myparams = (struct params*);
	
	
	
int main( int argc, char *argv[] )
{
	char c;
	
	 //pthread_t mythread[nthread];                 //array defining number of threads
	 //struct params * p1;
	
	// These are the default configuration values used
	// if no command line arguments are given.

	const char *outfile = "mandel.bmp";
	double xcenter = 0;
	double ycenter = 0;
	double scale = 4;
	int    image_width = 500;
	int    image_height = 500;
	int max = 1000;
	int nthread= 1;    //default thread is 1
	
	int start_image=0;
	int end_image=0;
	
	//int i;
	

	
	// For each command line argument given,
	// override the appropriate configuration value.

	while((c = getopt(argc,argv,"x:y:s:W:H:m:n:o:h"))!=-1) {                       //nthreads from the user 
		switch(c) {
			case 'x':
				xcenter = atof(optarg);
				break;
			case 'y':
				ycenter = atof(optarg);
				break;
			case 's':
				scale = atof(optarg);
				break;
			case 'W':
				image_width = atoi(optarg);
				break;
			case 'H':
				image_height = atoi(optarg);
				break;
			case 'm':
				max = atoi(optarg);
				break;
			case 'n':
				nthread = atoi(optarg);
				break;
			case 'o':
				outfile = optarg;
				break;
			case 'h':
				show_help();
				exit(1);
				break;
		}
	}
	
	

	// Display the configuration of the image.
	printf("mandel: x=%lf y=%lf scale=%lf max=%d nthread= %d outfile=%s\n",xcenter,ycenter,scale,max,nthread,outfile);

	// Create a bitmap of the appropriate size.
	struct bitmap *bm = bitmap_create(image_width,image_height);
	
	struct params p1[nthread];
	
	//struct params* ptemp1 = (struct params*) p1;
	pthread_t mythread[nthread];                    //array defining number of threads
    int img= image_height/nthread;
	int i;
	
	double time_spent = 0.0;    //clock for check
	clock_t begin = clock();
	for (i = 0; i < nthread; i++)
	{
	
	p1[i].bm= bm;
	p1[i].xmin=xcenter-scale;
    p1[i].xmax=xcenter+scale;
	p1[i].ymin=ycenter-scale;
    p1[i].ymax=ycenter+scale; 
	p1[i].max=max;
	p1[i].start_image=start_image;
    p1[i].end_image=end_image;
	p1[i].nthread=nthread;
	

	// Fill it with a dark blue, for debugging
	//bitmap_reset(bm,MAKE_RGBA(0,0,255,0));

	// Compute the Mandelbrot image
	//compute_image(bm,xcenter-scale,xcenter+scale,ycenter-scale,ycenter+scale,max,nthread);
	
		p1[i].start_image= (i)*(img);
		p1[i].end_image = ((i+1)*(img));
		
		
		pthread_create(&mythread[i],NULL, compute_image,(void*) &p1[i]);     //creating a thread
		
	}
	
	for (i=0; i<nthread; i++ )
	{
		
		pthread_join(mythread[i], NULL);
	}
	
	clock_t end = clock();
	time_spent += (double)(end- begin) / CLOCKS_PER_SEC;
	
	printf("Time Spent %f seconds", time_spent);

	// Save the image in the stated file.
	if(!bitmap_save(bm,outfile)) {
		fprintf(stderr,"mandel: couldn't write to %s: %s\n",outfile,strerror(errno));
		return 1;
	}

	return 0;
}

/*
Compute an entire Mandelbrot image, writing each point to the given bitmap.
Scale the image to the range (xmin-xmax,ymin-ymax), limiting iterations to "max"
*/


void *compute_image(void* p)
{
	struct params* ptemp = (struct params*) p;
	
	
	int width = bitmap_width(ptemp->bm);
	int height = bitmap_height(ptemp->bm);

	// For every pixel in the image...

	//for(j=0;j<height;j++) {


    
	int i,j;

			
		for(j= ptemp->start_image; j< ptemp->end_image; j++)
        {
		for(i=0;i<width;i++) {

			// Determine the point in x,y space for that pixel.
			double x = ptemp->xmin + i*(ptemp->xmax-ptemp->xmin)/width;
			double y = ptemp->ymin + j*(ptemp->ymax-ptemp->ymin)/height;

			// Compute the iterations at that point.
			int iters = iterations_at_point(x,y,ptemp->max);

			// Set the pixel in the bitmap.
			bitmap_set(ptemp->bm,i,j,iters);
		}
	    }
		
		
		
    //}
}

/*
Return the number of iterations at point x, y
in the Mandelbrot space, up to a maximum of max.
*/

int iterations_at_point( double x, double y, int max )
{
	double x0 = x;
	double y0 = y;
	

	int iter = 0;

	while( (x*x + y*y <= 4) && iter < max ) {

		double xt = x*x - y*y + x0;
		double yt = 2*x*y + y0;

		x = xt;
		y = yt;

		iter++;
	}

	return iteration_to_color(iter,max);
}

/*
Convert a iteration number to an RGBA color.
Here, we just scale to gray with a maximum of imax.
Modify this function to make more interesting colors.
*/

int iteration_to_color( int i, int max )
{
	int gray = 255*i/max;
	return MAKE_RGBA(gray,gray,gray,0);
}
