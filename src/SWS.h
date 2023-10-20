#pragma once
#ifndef SWS_H
#define SWS_H

//#include "init_visual.h"
#include "baseline.h"

struct index_time_pair
{
	int index;
	double time;
	bool operator<(const index_time_pair& pair) const { return time < pair.time; }
};

struct win_status
{
	//double s_bandwidth;
	//double t_bandwidth;

	//prev window
	double start_t_win_val_prev;
	double end_t_win_val_prev;

	//cur window
	double start_t_win_val;
	double end_t_win_val;
	int start_t_win_pos;
	int	end_t_win_pos;
	int center_t_win_pos; //Used in triangular kernel
};

void SWS_visual(statistics& stat);

#endif