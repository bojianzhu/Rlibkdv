#include "baseline.h"

double SCAN_2D(double*q, statistics& stat)
{
	double temp_value;
	double sq_dist_value;
	//double dist_value;
	double incr_value = 0;

	for (int i = 0; i < stat.n; i++)
	{
		sq_dist_value = sq_euclid_dist(q, stat.featureVector[i], stat.dim);

		if (sq_dist_value > stat.bandwidth_s*stat.bandwidth_s)
			continue;

		if (stat.kernel_s_type == 0) //Uniform kernel
			incr_value += stat.weightVector[i] * (1.0 - (1.0 / stat.bandwidth_s));
		if (stat.kernel_s_type == 1) //Epanechnikov kernel
			incr_value += stat.weightVector[i] * (1.0 - (1.0 / (stat.bandwidth_s*stat.bandwidth_s))*sq_dist_value);
		if (stat.kernel_s_type == 2) //Quartic kernel
		{
			temp_value = 1.0 - (1.0 / (stat.bandwidth_s*stat.bandwidth_s))*sq_dist_value;
			incr_value += stat.weightVector[i] * (temp_value * temp_value);
		}
	}

	return incr_value;
}

void SCAN_visual(statistics& stat)
{
	double*q;
	for (int c = 0; c < stat.col_pixels; c++)
	{
		for (int r = 0; r < stat.row_pixels; r++)
		{
			q = stat.queryVector[r*stat.col_pixels + c];
			stat.outMatrix[r][c] = SCAN_2D(q, stat);
		}
	}
}

void SCAN_otf_STKDV_visual(statistics& stat)
{
	double*q = new double[stat.dim];

	q[2] = stat.cur_time;
	for (int r = 0; r < stat.row_pixels; r++)
	{
		q[0] = stat.x_L + r * stat.incr_x;
		for (int c = 0; c < stat.col_pixels; c++)
		{
			q[1] = stat.y_L + c * stat.incr_y;
			stat.outMatrix[r][c] = 0;

			for (int i = 0; i < stat.n; i++)
				stat.outMatrix[r][c] += stat.weightVector[i] * spatial_kernel(q, stat.featureVector[i], stat) * temporal_kernel(q, stat.featureVector[i], stat);
		}
	}
}

void SCAN_batch_STKDV_visual(statistics& stat)
{
	double*q = new double[stat.dim];

	for (int r = 0; r < stat.row_pixels; r++)
	{
		q[0] = stat.x_L + r * stat.incr_x;
		for (int c = 0; c < stat.col_pixels; c++)
		{
			q[1] = stat.y_L + c * stat.incr_y;
			for (int t = 0; t < stat.t_pixels; t++)
			{
				q[2] = stat.t_L + t * stat.incr_t;
				stat.outCube[r][c][t] = 0;
				for (int i = 0; i < stat.n; i++)
					stat.outCube[r][c][t] += stat.weightVector[i] * spatial_kernel(q, stat.featureVector[i], stat) * temporal_kernel(q, stat.featureVector[i], stat);
			}
		}
	}
}

double spatial_kernel(double*q, double*p, statistics& stat)
{
	double value;

	if (stat.kernel_s_type == 1) //Epanechnikov kernel
	{
		value = 1 - ((q[0] - p[0])*(q[0] - p[0]) + (q[1] - p[1])*(q[1] - p[1])) / (stat.bandwidth_s*stat.bandwidth_s);

		if (value < 0)
			return 0;
		else
			return value;
	}

	if (stat.kernel_s_type == 2) //Quartic kernel
	{
		value = 1 - ((q[0] - p[0])*(q[0] - p[0]) + (q[1] - p[1])*(q[1] - p[1])) / (stat.bandwidth_s*stat.bandwidth_s);

		if (value < 0)
			return 0;
		else
			return value * value;
	}

	if (stat.kernel_s_type == 3) //Triangular kernel
	{
		value = 1 - fabs(sqrt((q[0] - p[0])*(q[0] - p[0]) + (q[1] - p[1])*(q[1] - p[1]))) / stat.bandwidth_s;

		if (value < 0)
			return 0;
		else
			return value;
	}

	return -inf;
}

double temporal_kernel(double*q, double*p, statistics& stat)
{
	double value;

	if (stat.kernel_t_type == 1) //Epanechnikov kernel
	{
		value = 1 - ((q[2] - p[2])*(q[2] - p[2])) / (stat.bandwidth_t*stat.bandwidth_t);
		if (value < 0)
			return 0;
		else
			return value;
	}

	if (stat.kernel_t_type == 2) //Quartic kernel
	{
		value = 1 - ((q[2] - p[2])*(q[2] - p[2])) / (stat.bandwidth_t*stat.bandwidth_t);
		if (value < 0)
			return 0;
		else
			return value * value;
	}

	if (stat.kernel_t_type == 3) //Triangular kernel
	{
		value = 1 - fabs(q[2] - p[2]) / stat.bandwidth_t;
		if (value < 0)
			return 0;
		else
			return value;
	}

	return -inf;
}