#pragma once
#ifndef SLAM_H
#define SLAM_H

#include "init_visual.h"
#include "bucket.h"

/*struct bound_entry
{
	int id;
	double bound_value;
	bool is_LB;

	bool operator < (const bound_entry& another) const {
		return bound_value < another.bound_value;
	}
};*/

void SLAM_visual(statistics& stat);

#endif