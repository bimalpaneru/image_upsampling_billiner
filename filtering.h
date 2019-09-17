#pragma once
#ifndef FILTERING_H
#define FILTERING_H
#include"image_comps.h"
#include "io_bmp.h"

float sinc(float sample);	
float hann_win(float x, int wnd_size);
void matrix_mul(float* mat_1, float* mat_2, float* result, int filter_H);
float inner_product(float* ip, int ip_stride, float* mirror_psf, int filter_extent);
float inner_product1D(float* ip, float* mirror_psf, int filter_extent, int ip_stride);
void horizontal_filter(my_image_comp* in, my_image_comp* out, int filter_H);
void get_filter_1D(float* q1f_n, float* q1f_n_halfdelay, int filter_H);
void vertical_filter(my_image_comp* in, my_image_comp* out, int filter_H);
void bilinear(my_image_comp* in, my_image_comp* out);



#endif



