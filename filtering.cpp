#include"filtering.h"
#include<stdio.h>
#include<math.h>
#include<cmath>
#include "image_comps.h"

const float PI = 3.14F;
float sinc(float sample) {
	if (sample == 0)
		return 0.4F;
	else
		return (0.4F * sinf(PI * sample * 0.4F) / (PI * sample * 0.4F));

}

float hann_win(float x, int wnd_size) {
	float sum = 0.5F * (1.0F + cosf(PI * x / float(wnd_size)));
	return 0.5F * (1.0F + cosf(PI * x / float(wnd_size)));
}


void matrix_mul(float* mat_1, float* mat_2, float* result, int filter_H) {

	for (int x = 0; x <= 2 * filter_H; x++) {
		for (int y = 0; y <= 2 * filter_H; y++) {
			result[x * (2 * filter_H + 1) + y] = mat_1[x] * mat_2[y];

		}

	}

}

float inner_product(float* ip, int ip_stride, float* mirror_psf, int filter_extent) {
	float sum = 0.0F;

	for (int y = -filter_extent; y <= filter_extent; y++) {
		for (int x = -filter_extent; x <= filter_extent; x++) {
			ip[y * ip_stride + x];
			mirror_psf[y * (2 * filter_extent + 1) + x];
			sum += ip[y * ip_stride + x] * mirror_psf[y * (2 * filter_extent + 1) + x];
		}
	}
	return sum;
}


float inner_product1D(float* ip, float* mirror_psf, int filter_extent, int ip_stride) {
	float sum = 0.0F;
	for (int i = -filter_extent; i <= filter_extent; i++) {
		sum += ip [i * ip_stride] * mirror_psf[i];

	}
	return sum;

}




void horizontal_filter(my_image_comp* in, my_image_comp* out, int filter_H) {

	float* q1f_n = new float[2 * filter_H + 1];   //windowd sample of streatched sinc function for even positions
	float* q1f_n_halfdelay = new float[2 * filter_H + 1];//windowd sample of streatched sinc function for odd positions
	get_filter_1D(q1f_n, q1f_n_halfdelay, filter_H);
	
	//-----make the array point at the origin of the filter----//
	q1f_n = q1f_n + filter_H;
	q1f_n_halfdelay = q1f_n_halfdelay + filter_H;
	//-----make the array point at the origin of the filter----//

	for (int r = 0; r < out->height; r++) {

		for (int c = 0; c < out->width; c++) {

			float* op = out->buf + r * out->stride + c;

			//check for even numbered rows in intermediate image
			//doing inner product of every element of corrosponding row in original image 
			//and map it to intermediate image 
			//correspoinding row in original image  = 5*r/2 
			/*if (r % 2 == 0) {
				float* ip = in->buf + r * 5 * in->stride / 2 + c;
				*op = inner_product1D(ip, q1f_n, filter_H, 1);

			}*/

			if (c % 2 == 0) {
				float* ip = in->buf + r * in->stride + c * 5 / 2;
				*op = inner_product1D(ip, q1f_n, filter_H, 1);

			}



			//for odd numbered row in intermediate image
			//doing inner product of every element of corrosponding row in original image 
			//and map it to intermediate image 
			//correspoinding row in original  = 5*(r-1)/2 + 2
			/*else {
				float* ip = in->buf + (5 * (r - 1) / 2 + 2) * in->stride + c;
				*op = inner_product1D(ip, q1f_n, filter_H, 1);
			}*/

			else {
				float* ip = in->buf + r * in->stride + (5 * (c - 1) / 2 + 2);
			*op = inner_product1D(ip, q1f_n, filter_H, 1);
			}



		}

	}

}

void get_filter_1D(float* q1f_n, float* q1f_n_halfdelay, int filter_H) {
	float gain_q1f = 0.0F;
	float gain_q1f_halfdelay = 0.0F;
	for (int i = 0; i <= 2 * filter_H; i++) {

		q1f_n[i] = sinc(float(i) - filter_H) * hann_win((float(i) - filter_H), filter_H + 1);
		gain_q1f += q1f_n[i];
	}
	//q1f_n;
	//get the values of windowed sinc function delayed by half sample
	for (int i = 0; i <= 2 * filter_H; i++) {

		q1f_n_halfdelay[i] = sinc(float(i) - filter_H - 0.5F) * hann_win((float(i) - filter_H - 0.5F), filter_H + 1);
		gain_q1f_halfdelay += q1f_n_halfdelay[i];
	}

	for (int i = 0; i <= 2 * filter_H; i++) {
		q1f_n[i] = q1f_n[i] / gain_q1f;
		q1f_n_halfdelay[i] = q1f_n_halfdelay[i] / gain_q1f_halfdelay;
	}


}




void vertical_filter(my_image_comp* in, my_image_comp* out, int filter_H) {
	float* q1f_n = new float[2 * filter_H + 1];   //windowd sample of streatched sinc function for even positions
	float* q1f_n_halfdelay = new float[2 * filter_H + 1];//windowd sample of streatched sinc function for odd positions

	get_filter_1D(q1f_n, q1f_n_halfdelay, filter_H);


	//-----make the array point at the origin of the filter----//
	q1f_n = q1f_n + filter_H;
	q1f_n_halfdelay = q1f_n_halfdelay + filter_H;
	//-----make the array point at the origin of the filter----//

	for (int r = 0; r < out->height; r++) {

		for (int c = 0; c < out->width; c++) {

			float* op = out->buf + r * out->stride + c;
			
			/*if (c % 2 == 0) {
				float* ip = in->buf + r * in->stride + c * 5 / 2;
				*op = inner_product1D(ip, q1f_n, filter_H, in->stride);

			}*/
			if (r % 2 == 0) {
				float* ip = in->buf + r * 5 * in->stride / 2 + c;
				*op = inner_product1D(ip, q1f_n, filter_H, in->stride);

			}


			/*else {
				float* ip = in->buf + r * in->stride + 5 * (c - 1) / 2 + 2;
				*op = inner_product1D(ip, q1f_n_halfdelay, filter_H, in->stride);
			}*/

			else {
				float* ip = in->buf + (5 * (r - 1) / 2 + 2) * in->stride + c;
				*op = inner_product1D(ip, q1f_n, filter_H, in->stride);
			}

		}
	}
}





void apply_filter(my_image_comp* in, my_image_comp* out, int filter_H)
{

	//filter length H in 1D then 2H+1 is 1D filter length

	float* q1f_n = new float[2 * filter_H + 1];   //windowd sample of streatched sinc function for even positions
	float* q1f_n_halfdelay = new float[2 * filter_H + 1];//windowd sample of streatched sinc function for odd positions

	float* rowE_colE = new float[(2 * filter_H + 1) * (2 * filter_H + 1)];
	float* rowE_colO = new float[(2 * filter_H + 1) * (2 * filter_H + 1)];
	float* rowO_colO = new float[(2 * filter_H + 1) * (2 * filter_H + 1)];
	float* rowO_colE = new float[(2 * filter_H + 1) * (2 * filter_H + 1)];

	//get the values of windowed sinc function
	for (int i = 0; i <= 2 * filter_H; i++) {

		q1f_n[i] = sinc(float(i) - filter_H) * hann_win((float(i) - filter_H), filter_H + 1);

	}
	q1f_n;
	//get the values of windowed sinc function delayed by half sample
	for (int i = 0; i <= 2 * filter_H; i++) {

		q1f_n_halfdelay[i] = sinc(float(i) - filter_H - 0.5F) * hann_win((float(i) - filter_H - 0.5F), filter_H + 1);

	}
	q1f_n_halfdelay;

	//multiply to get filter matrices
	matrix_mul(q1f_n, q1f_n, rowE_colE, filter_H);
	matrix_mul(q1f_n, q1f_n_halfdelay, rowE_colO, filter_H);
	matrix_mul(q1f_n_halfdelay, q1f_n_halfdelay, rowO_colO, filter_H);
	matrix_mul(q1f_n_halfdelay, q1f_n, rowO_colE, filter_H);



	//------points at origin of  the filter matrix-----//
	rowE_colE = rowE_colE + (2 * filter_H + 1) * filter_H + filter_H;
	rowE_colO = rowE_colO + (2 * filter_H + 1) * filter_H + filter_H;

	rowO_colO = rowO_colO + (2 * filter_H + 1) * filter_H + filter_H;
	rowO_colE = rowO_colE + (2 * filter_H + 1) * filter_H + filter_H;
	//------points at origin of  the filter matrix-----//

	float s = 0.0F;
	filter_H;
	for (int x = 0; x <= 2 * filter_H; x++) {

		for (int y = 0; y <= 2 * filter_H; y++) {

			s += abs(rowE_colE[x * (2 * filter_H + 1)]);
			printf("%f *****", rowE_colE[x * (2 * filter_H + 1) + y]);

		}
		printf("\n");

	}

	//printf("\n Sum = %f *****", s); 


  // Check for consistent dimensions
	assert(in->border >= filter_H);
	//assert((out->height <= new_width) && (out->width <= in->width));

	// Perform the convolution
	for (int r = 0; r < out->height; r++) {

		for (int c = 0; c < out->width; c++) {

			float* op = out->buf + r * out->stride + c;

			if (r % 2 == 0 && c % 2 == 0) {
				float* ip = in->buf + 5 * r * in->stride / 2 + 5 * c / 2; //point inital position of the component buffer
				*op = inner_product(ip, in->stride, rowE_colE, filter_H);

			}

			else if (r % 2 == 0 && c % 2 != 0) {
				float* ip = in->buf + 5 * r / 2 * in->stride + 5 * c / 2 + 2;
				*op = inner_product(ip, in->stride, rowE_colO, filter_H);

			}


			else if (r % 2 != 0 && c % 2 != 0) {
				float* ip = in->buf + (5 * r / 2 + 2) * in->stride + (5 * c / 2 + 2);
				*op = inner_product(ip, in->stride, rowO_colO, filter_H);

			}

			else if (r % 2 != 0 && c % 2 == 0) {
				float* ip = in->buf + (5 * r / 2 + 2) * in->stride + 5 * c / 2;
				*op = inner_product(ip, in->stride, rowE_colO, filter_H);

			}


		}

	}

}

void bilinear(my_image_comp* in, my_image_comp* out) {
	float exp_ratio = 0.4F;
	float* a, * b, * c, * d;//, alpha, beta;
	float w, h;


	for (int y = 0; y < out->height; y++) {

		for (int x = 0; x < out->width; x++) {
			float* op = out->buf + y * out->stride + x;
			a = in->buf + int(y * exp_ratio) * in->stride + int(x * exp_ratio);
			b = a + 1;
			c = a + in->stride;
			d = c + 1;
			w = exp_ratio * x - int(exp_ratio * x);
			h = exp_ratio * y - int(exp_ratio * y);

			//alpha = *a + del * (*b - *a);
			//beta = *c + del * (*d - *c);

			*op = *a * (1 - w) * (1 - h) + *b * (w) * (1 - h) + *c * (h) * (1 - w) + *d * (w * h);

		}

	}




}