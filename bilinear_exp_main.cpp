/*****************************************************************************/
// File: filtering_main.cpp
// Author: David Taubman
// Last Revised: 13 August, 2007
/*****************************************************************************/
// Copyright 2007, David Taubman, The University of New South Wales (UNSW)
/*****************************************************************************/

#include "io_bmp.h"
#include "image_comps.h"
#include<math.h>
#include<cmath>
#include "filtering.h"
#include<time.h>
/* ========================================================================= */
/*                 Implementation of `my_image_comp' functions               */
/* ========================================================================= */

/*****************************************************************************/
/*                  my_image_comp::perform_boundary_extension                */
/*****************************************************************************/




void my_image_comp::perform_boundary_extension()
{
  int r, c;
  
  // First extend upwards
  float *first_line = buf;
  for (r = 1; r <= border; r++)
	  for (c = 0; c < width; c++)
		  first_line[-r * stride + c] = first_line[(r - 1) * stride + c];//first_line[c];

  // Now extend downwards
  float *last_line = buf+(height-1)*stride;
  for (r = 1; r <= border; r++)
	  for (c = 0; c < width; c++)
		  last_line[r * stride + c] = last_line[-(r - 1) * stride + c];//[(height - 1) * -r * stride + c];    //     r * stride + c];//0;//last_line[c];

  // Now extend all rows to the left and to the right
  float *left_edge = buf-border*stride;
  float *right_edge = left_edge + width - 1;
  for (r=height+2*border; r > 0; r--, left_edge+=stride, right_edge+=stride)
    for (c=1; c <= border; c++)
      {
		left_edge[-c] =  left_edge[c];// left_edge[0];
		right_edge[c] =  right_edge[-c];// right_edge[0];
      }
}


/* ========================================================================= */
/*                              Global Functions                             */
/* ========================================================================= */

/*****************************************************************************/
/*                                apply_filter                               */
/*****************************************************************************/




/*****************************************************************************/
/*                                    main                                   */
/*****************************************************************************/

int  main(int argc, char *argv[])
{
	clock_t start_time = clock();
	
	//float seconds;
 

  if (argc != 4)
    {
      fprintf(stderr,"Usage: %s <in bmp file> <out bmp file>\n",argv[0]);
      return -1;
    }
  int filter_H = atoi(argv[3]);
  //printf("filter length %d ", filter_H);
 
  int err_code=0;
  try {
      // Read the input image
      bmp_in in;
      if ((err_code = bmp_in__open(&in,argv[1])) != 0)
        throw err_code;

	  
      int width = in.cols, height = in.rows;

	  int new_width = int(ceil(float(width) * 2.5F));
	  int new_height = int(ceil(float(height) * 2.5F));

      int n, num_comps = in.num_components;
      my_image_comp *input_comps = new my_image_comp[num_comps];
	  // Allocate storage for the filtered output
	  my_image_comp* output_comps = new my_image_comp[num_comps];
	  //my_image_comp* intermediate_comps = new my_image_comp[num_comps];
      
	  ///extent boundaries
	  for (n = 0; n < num_comps; n++) {
		  input_comps[n].init(height, width, filter_H); // Leave a border of filter extent size
	  }
        
      

      int r; // Declare row index
      io_byte *line = new io_byte[width*num_comps];
      for (r=height-1; r >= 0; r--)
        { // "r" holds the true row index we are reading, since the image is
          // stored upside down in the BMP file.
          if ((err_code = bmp_in__get_line(&in,line)) != 0)
            throw err_code;
          for (n=0; n < num_comps; n++)
            {
              io_byte *src = line+n; // Points to first sample of component n
              float *dst = input_comps[n].buf + r * input_comps[n].stride;
              for (int c=0; c < width; c++, src+=num_comps)
                dst[c] = (float) *src; // The cast to type "float" is not
                      // strictly required here, since bytes can always be
                      // converted to floats without any loss of information.
            }
        }
      bmp_in__close(&in);

      
	  for (n = 0; n < num_comps; n++) {

		  input_comps[n].perform_boundary_extension(); //extent ip image boundary
		 ///intermediate_comps[n].init(height, new_width, filter_H); // need a border equal to filter extent
		 // intermediate_comps[n].perform_boundary_extension();
		  output_comps[n].init(new_height, new_width, 0); // dont need a border for final output

		  //apply filters
		  
		  		 
	  }
	  
	  for (n = 0; n < num_comps; n++) {
		  bilinear(input_comps + n, output_comps + n);
		  //  horizontal_filter(input_comps + n, intermediate_comps + n, filter_H);
		 // intermediate_comps[n].perform_boundary_extension();
		 // vertical_filter(intermediate_comps + n, output_comps + n, filter_H);

	  }
	 

	 

        

      // Process the image, all in floating point (easy)
      //for (n=0; n < num_comps; n++)
        //input_comps[n].perform_boundary_extension();
      //for (n=0; n < num_comps; n++)
        //apply_filter(input_comps+n,output_comps+n, filter_H);

      // Write the image back out again
      bmp_out out;
      if ((err_code = bmp_out__open(&out,argv[2],new_width,new_height,num_comps)) != 0)
        throw err_code;
      for (r=new_height-1; r >= 0; r--)
        { // "r" holds the true row index we are writing, since the image is
          // written upside down in BMP files.
          for (n=0; n < num_comps; n++)
            {
              io_byte *dst = line+n; // Points to first sample of component n
              float *src = output_comps[n].buf + r * output_comps[n].stride;
			  for (int c = 0; c < new_width; c++, dst += num_comps)
			  {
				  if (src[c] > 255.0F)
					  src[c] = 255;
				  else if (src[c] <0.0F)
					  src[c] = 0;

				  *dst = (io_byte)src[c];
			  } // The cast to type "io_byte" is
                      // required here, since floats cannot generally be
                      // converted to bytes without loss of information.  The
                      // compiler will warn you of this if you remove the cast.
                      // There is in fact not the best way to do the
                      // conversion.  You should fix it up in the lab.
            }
          bmp_out__put_line(&out,line);
        }
      bmp_out__close(&out);
      delete[] line;
      delete[] input_comps;
      delete[] output_comps;
    }
  catch (int exc) {
      if (exc == IO_ERR_NO_FILE)
        fprintf(stderr,"Cannot open supplied input or output file.\n");
      else if (exc == IO_ERR_FILE_HEADER)
        fprintf(stderr,"Error encountered while parsing BMP file header.\n");
      else if (exc == IO_ERR_UNSUPPORTED)
        fprintf(stderr,"Input uses an unsupported BMP file format.\n  Current "
                "simple example supports only 8-bit and 24-bit data.\n");
      else if (exc == IO_ERR_FILE_TRUNC)
        fprintf(stderr,"Input or output file truncated unexpectedly.\n");
      else if (exc == IO_ERR_FILE_NOT_OPEN)
        fprintf(stderr,"Trying to access a file which is not open!(?)\n");
      return -1;
    }
  clock_t end_time = clock();
  float second = ((float)(end_time - start_time)) / CLOCKS_PER_SEC;
  printf("Computation Time %f \n ", second);
  return 0;
}
