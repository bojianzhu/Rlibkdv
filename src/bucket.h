#pragma once
#ifndef BUCKET_H
#define BUCKET_H

#include "init_visual.h"

struct bound_entry
{
	int id;
	double bound_value;
	bool is_LB;

	bool operator < (const bound_entry& another) const {
		return bound_value < another.bound_value;
	}
};

void init_Bucket(statistics& stat);
void clear_Bucket(statistics& stat);
void erase_Bucket(statistics& stat);
void bucket_algorithm_row(statistics& stat, double k, int y_index, vector<double**>& stat_plane_vec);
void bucket_algorithm(statistics& stat, vector<double**>& stat_plane_vec);

#endif
