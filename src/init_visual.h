#pragma once
#ifndef INIT_VISUAL_H
#define INIT_VISUAL_H

#include "Library.h"

const double inf = 1e80;
const double small_epsilon = 0.0001;
const double debug_small_epsilon = 1e-10;
const double pi = 3.14159265358979323846;
const double earth_radius = 6371000; //6371000 meters

struct SLAM
{
	double q_SquareNorm;
	//vector<int> L_ell;
	//vector<int> U_ell;
	double W_L_ell;
	double W_U_ell;
	double*A_L_ell;
	double*A_U_ell;
	double S_L_ell;
	double S_U_ell;
	double k;
	vector<double*> query_list;
	vector<double> result_list;
	double W_q;
	//double R_q_card;
	double*A_R_q;
	double S_R_q;
	double max_KDE;

	//double**outMatrix_SLAM_obj;
};

struct SWS
{
	double max_KDE;
	double*q;
	double*sliding_window;
	double*sliding_window_L; //Used in triangular kernel
	double*sliding_window_R; //Used in triangular kernel
};

struct statistics
{
	int n; //number of points in a dataset (after filtering)
	double bandwidth_s; //spatial bandwidth value
	double bandwidth_t; //temporal bandwidth value
	double time_q; //Used in online STKDV
	char*dataFileName_JSON; //Stores (x, y, time, w) or (x, y, w)
	char*dataFileName_CSV; //Stores (x, y, time, w) or (x, y, w)
	vector<double*> base_dataMatrix; //This base_dataMatrix should not be deleted
	vector<double> base_weightVector; //This base_weightVector should not be deleted
	vector<double*> featureVector; //feature vector of all data points (after the filter operation)
	vector<double> weightVector; //weight value for all data points
	double**queryVector;

	string outputFileName; //output file Name
	//dim = 2 (Used for KDV and online STKDV) 
	//dim = 3 (Used for batch-based STKDV and batch-based bandwidth exploration in KDV)
	int dim;

	double x_L; double x_U;
	double y_L; double y_U;
	double t_L; double t_U; //Used in STKDV (batch-based KDV)
	double incr_x; double incr_y; double incr_t;

	//0: Uniform kernel 
	//1: Epanechnikov kernel
	//2: Quartic kernel
	//3: Triangular kernel
	int kernel_s_type; //spatial kernel 
	int kernel_t_type; //temporal kernel

	//Bucket-based algorithm
	vector< vector<int> > B_L;
	vector< vector<int> > B_U;
	double*C_L; double*C_U;
	double**v_L; double**v_U;
	double*H_L; double*H_U;
	double*C_R_q;
	double**v_R_q;
	double*H_R_q;

	//(1) KDV 
	//(2) online STKDV  
	//(3) batch-based STKDV
	int KDV_type;

	//Multithreading of SLAM
	int num_threads;

	//Used in SLAM (SLAM_SORT)
	int static_coord;
	int dynamic_coord;
	int static_pixel_size;
	int dynamic_pixel_size;
	vector<SLAM> SLAM_vec;

	//Used in SWS, EDWIN_single, and EDWIN_multiple
	double*q;
	double**sorted_featureVector;

	//Used in SWS
	double*sorted_weightVector;
	vector<SWS> SWS_vec;

	//Used in EDWIN_single and EDWIN_multiple
	vector<double**> S_plane_vec;
	vector<double> sorted_fV_timestamp_vec;
	int start_t_win_pos;
	int end_t_win_pos;

	//Used in EDWIN_single and EDWIN_multiple
	double max_EDWIN_KDE;

	//Used in EDWIN_single
	double cur_time;

	//Used in EDWIN_multiple
	vector<double**> S_D_plane_vec;
	vector<double**> S_I_plane_vec;
	int start_S_D_t_win_pos;
	int end_S_D_t_win_pos;
	int start_S_I_t_win_pos;
	int end_S_I_t_win_pos;

	//output visualization
	int row_pixels;	int col_pixels; int t_pixels;
	double**outMatrix;
	double***outCube;
	//int out_r; int out_c; int out_t;
};

//void GPS_to_x_y(double longitude, double latitude, double& x, double& y, statistics& stat);
//void x_y_to_GPS(double x, double y, double& longitude, double& latitude, statistics& stat);
double computeSqNorm(double*q, int dim);
double inner_product(double*q, double*p, int dim);
double sq_euclid_dist(double*q, double*p, int dim);
void initQuery(statistics& stat); //KDV
void update_incr_values(statistics& stat); //Batch-based STKDV

#endif