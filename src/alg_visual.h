#pragma once
#ifndef ALG_VISUAL_H
#define ALG_VISUAL_H

#include "baseline.h"
#include "SLAM.h"
#include "SWS.h"
#include "EDWIN_otf.h"
#include "EDWIN_multiple.h"

class alg_visual
{
public:
	void load_parameters(int argc, char**argv);
	void filter_datasets();
	void init_visual();
	void visual_Algorithm();
	void matrix_normalization(double max_KDE);
	void cube_normalization(double max_KDE);
	string saveMatrix_toString_CSV();
	string saveCube_toString_CSV();
	string saveMatrix_toString();
	string saveCube_toString();
	string compute(int argc, char**argv);
	void clear_memory();
	void clear_basic_memory();

	//Used for testing
	void obtain_L_U();
	void output_File();
	//void output_matrix();
	//void output_cube();
	//void obtain_long_lat_t_L_U();
	//void compute_test(int argc, char**argv);

	//Call it once
	void load_datasets(char**argv);
	void load_datasets_CSV(char**argv);

private:
	statistics stat;
};

#endif