#include <cstring>
#include <cmath>
#include "stdafx.h"
#include "stdio.h"
#include <string>
#include <vector>
#include "stb_image.h"//open source utility to load jpeg format images
#include "stb_image_write.h"//open source utility to save rgb images
#include <Windows.h>
#include <MMSystem.h>
#pragma comment(lib, "winmm.lib")


#define CHANNEL_NUM 3//RBG


using namespace std;


int string_to_int(const string &a){
	int ret=0;
	sscanf_s(a.c_str(),"%d",&ret);
	return ret;
}

void load_image_greyscale(	const uint8_t *cur_image,
							uint8_t *binary_image,
							const int height,
							const int width);

void split_image_m_n_subimages(	uint8_t *binary_image,
								const int height,
								const int width,
								const int m,
								const int n,
								const int pixel_count,
								vector< vector<uint8_t> > &vec);

/*
From command line:<application> <jpeg image file name> <m> <n> <pixel count>
will load file and save the sub-image set to the current directory in both binary format
and as a png file. 
The threshold approach for the conversion from 3 channel (8-bit per channel) to single 
channel binary was to convert the pixel values to a single floating point value based 
on the //Rec. 601 luma formula, https://en.wikipedia.org/wiki/Luma_(video). Then the
threshold value was simply the midpoint between the min and max for the entire image.
While this is the most basic approach, it is used in the MATLAB implementation so 
can serve as a good starting approach. A more robust method would be the Otsu method 
available in Python; https://scikit-image.org/docs/dev/auto_examples/segmentation/plot_thresholding.html.
A 'from scratch' C++ Otsu implementation is overkill for this exercise, but could be done in a few hours.

Assumptions: 
A 'binary' image is an grayscale image where all values are either 0 or 1 (8-bit)
The first dimension (from left to right) is the image width (and corresponds to the sub image value 'm'),
and the second dimension is the image height (and corresponds to the sib image value 'n').

NOTE: Performance was not a consideration for the exercise. Converting some of the functions to CUDA
kernels would improve performance and removing some of the analysis loops would also improve performance.

Caveats: I mostly operate with CUDA, C, C++ and MATLAB so my coding style is a reflection of this work. While I do 
work with Python in the context of Tensorflow/Keras my personal opinion is that C++ with CUDA is usually
the best choice in terms of performance and stability. I am flexible but have learned from experience
that high-level abstraction often introduces more problematic issues than expected.

Oleg John Konings

*/


int main(int argc, char** argv){

	if(argc!=5){
		printf("\nApplication requires 4 inputs after appplication name; jpeg file name, m, n, pixel count.\n");
		exit(1);

	}


	printf("\nfile name image to load=%s \n",argv[1]);
	const int mm=string_to_int(argv[2]);
	const int nn=string_to_int(argv[3]);
	const int pixel_count=string_to_int(argv[4]);

	int width=0, height=0, bpp=0;

    uint8_t* rgb_image = stbi_load(argv[1], &width, &height, &bpp, 3);

	printf("\nwidth=%d,height=%d,bpp=%d, m=%d, n=%d, pixel_count=%d\n",width,height,bpp,mm,nn,pixel_count);
	const size_t elem_binary_image=(size_t)(height)*(size_t)(width);
	const size_t bytes_binary_image=elem_binary_image*sizeof(uint8_t);
	const size_t elem_sub_image=(size_t)(mm)*(size_t)(nn);
	const size_t bytes_sub_image=elem_sub_image*sizeof(uint8_t);

	uint8_t *H_binary_image=(uint8_t*)malloc(bytes_binary_image);

	load_image_greyscale(	rgb_image,
							H_binary_image,
							height,
							width);	

	vector< vector<uint8_t> > vec(0);

	split_image_m_n_subimages(	H_binary_image,
								height,
								width,
								mm,
								nn,
								pixel_count,
								vec);
	int num_images=(int)vec.size();
	printf("\nNumber of sub_images=%d, writing as binary image files and png format to current directory..\n",num_images);

	char name_buffer[1009];
	const int quality=50;
	FILE *pFile=NULL;
	uint8_t temp_rgb_val=0;
	uint8_t *temp_rgb_image=(uint8_t*)malloc(elem_sub_image*CHANNEL_NUM*sizeof(uint8_t));
	for(int i=0;i<vec.size();i++){
		sprintf_s(name_buffer,"%s_sub_image_%d.bin",argv[1],i);
		pFile=fopen(name_buffer,"wb");
		fwrite((uint8_t*)(&vec[i][0]),sizeof(uint8_t),elem_sub_image,pFile);
		fclose(pFile);
		sprintf_s(name_buffer,"%s_sub_image_%d.png",argv[1],i);
		for(int jj=0,jjj=0;jj<(mm*nn);jj++,jjj+=3){
			temp_rgb_val= vec[i][jj];
			temp_rgb_val= (temp_rgb_val>(uint8_t)0) ? (uint8_t)255:(uint8_t)0;
			temp_rgb_image[jjj]=temp_rgb_val;
			temp_rgb_image[jjj+1]=temp_rgb_val;
			temp_rgb_image[jjj+2]=temp_rgb_val;
		}
		stbi_write_png(name_buffer,mm,nn,CHANNEL_NUM,temp_rgb_image,mm*CHANNEL_NUM);

	}


    stbi_image_free(rgb_image);
	if(H_binary_image)free(H_binary_image);
	if(temp_rgb_image)free(temp_rgb_image);
	return 0;
}


void load_image_greyscale(	const uint8_t *cur_image,
							uint8_t *binary_image,
							const int height,
							const int width){

	memset(binary_image,0,height*width*sizeof(uint8_t));
	const size_t image_elem=(size_t)(height)*(size_t)(width)*3ULL;
	float R=0.0f, G=0.0f,B=0.0f,Y=0.0f;
	double mean_val=0.;
	float max_val=-1.0f, min_val=(float)(1<<29);

	int hist[256];
	memset(hist,0,sizeof(hist));

	for(size_t i=0,b_i=0;i<image_elem;i+=3ULL,b_i++){
		R=(float)(cur_image[i]);
		G=(float)(cur_image[i+1]);
		B=(float)(cur_image[i+2]);
		Y=0.2989f*R+0.5870f*G+0.1140f*B;//Rec. 601 luma formula, https://en.wikipedia.org/wiki/Luma_(video)
		mean_val+=(double)(Y);
		max_val= (Y>max_val) ? Y:max_val;
		min_val= (Y<min_val) ? Y:min_val;
		hist[min(255,(int)(Y))]++;

	}
	mean_val/=(double)(height*width);
	int max_hist=-1, idx=-1;
	for(int i=0;i<256;i++){
		if(hist[i]>max_hist){
			max_hist=hist[i];
			idx=i;
		}

	}
	printf("\nmax Y=%f, min Y=%f, mean Y=%lf, mode Y=%d, frequency mode=%d\n",max_val,min_val,mean_val,idx,max_hist);
	double stdDev=0.,variance=0.;
	for(size_t i=0,b_i=0;i<image_elem;i+=3ULL,b_i++){
		R=(float)(cur_image[i]);
		G=(float)(cur_image[i+1]);
		B=(float)(cur_image[i+2]);
		Y=0.2989f*R+0.5870f*G+0.1140f*B;//Rec. 601 luma formula, https://en.wikipedia.org/wiki/Luma_(video)
		variance+=pow(((double)(Y))-mean_val,2.0);
	}
	variance/=(double)(height*width);
	stdDev=sqrt(variance);
	float threshold_value=(max_val-min_val)/2.0f+min_val;//default threshold value will be halfway between the min and max of Y value for all image pixels

	printf("\nstdDev=%lf, threshold=%f\n",stdDev,threshold_value);

	for(size_t i=0,b_i=0;i<image_elem;i+=3ULL,b_i++){
		R=(float)(cur_image[i]);
		G=(float)(cur_image[i+1]);
		B=(float)(cur_image[i+2]);
		Y=0.2989f*R+0.5870f*G+0.1140f*B;//Rec. 601 luma formula, https://en.wikipedia.org/wiki/Luma_(video)

		binary_image[b_i]= (Y<threshold_value) ? 0:1;

	}

}

void split_image_m_n_subimages(	uint8_t *binary_image,
								const int height,
								const int width,
								const int m,
								const int n,
								const int pixel_count,
								vector< vector<uint8_t> > &vec){

	const size_t image_elem=(size_t)(height)*(size_t)(width);

	if(height%n || width%m || (m<=0) || (n<=0) ){
		printf("\nError! image height/width not evenly divisible by subimage size!\n");
		return;
	}
	int cur_pixel_count=0;
	vec.clear();
	vector<int> count_vec(0);


	for(int i=0;i<height;i+=n){

		for(int j=0;j<width;j+=m){

			cur_pixel_count=0;
			vector<uint8_t> sub_vec((m*n),0);
			int k=0;

			for(int ii=0;ii<n;ii++)for(int jj=0;jj<m;jj++){
				uint8_t cur=(uint8_t)(binary_image[(i+ii)*width+(j+jj)]);

				if(0==cur){
					sub_vec[k++]=cur;//not really necessary since vector should be all zeros upon init, but keeps the implementation simple

				}else{
					
					if(cur_pixel_count<pixel_count){//if the current non-zero pixel count is less than the max then add to current sub image
						sub_vec[k++]=cur;
						cur_pixel_count++;

					}else{//if the current count is at max (going one by one so if this branch is hit you are at the max exactly)
						vec.push_back(sub_vec);
						count_vec.push_back(cur_pixel_count);
						sub_vec.assign((m*n),0);
						sub_vec[k++]=cur;//keep the position so that the sub images can be added to get the full image with the total pixel count 
						cur_pixel_count=1;

					}
				}				
			}

			vec.push_back(sub_vec);//this will be be last sub-image
			count_vec.push_back(cur_pixel_count);//keeping track of pixel counts per sub image

		}
	}

	for(int i=0;i<count_vec.size();i++){
		printf("\nsub image %d has a pixel count of %d\n",i,count_vec[i]);
	}
}







