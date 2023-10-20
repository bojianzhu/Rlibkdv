#include "init_visual.h"

double computeSqNorm(double*q, int dim)
{
	double sqNorm = 0;
	for (int d = 0; d < dim; d++)
		sqNorm += q[d] * q[d];

	return sqNorm;
}

double inner_product(double*q, double*p, int dim)
{
	double ip = 0;

	for (int d = 0; d < dim; d++)
		ip += q[d] * p[d];

	return ip;
}

double sq_euclid_dist(double*q, double*p, int dim)
{
	double dist = 0;
	for (int d = 0; d < dim; d++)
		dist += (q[d] - p[d])*(q[d] - p[d]);
	return dist;
}

void initQuery(statistics& stat)
{
	int total_q = stat.row_pixels*stat.col_pixels;
	double x_coord;
	double y_coord;
	stat.queryVector = new double*[total_q];

	if (stat.row_pixels != 1 || stat.col_pixels != 1)
	{
		stat.incr_x = (stat.x_U - stat.x_L) / (stat.row_pixels - 1);
		stat.incr_y = (stat.y_U - stat.y_L) / (stat.col_pixels - 1);
	}

	if (stat.row_pixels == 1)
		stat.incr_x = 0;
	if (stat.col_pixels == 1)
		stat.incr_y = 0;

	for (int q = 0; q < total_q; q++)
		stat.queryVector[q] = new double[stat.dim];

	for (int r = 0; r < stat.row_pixels; r++)
	{
		x_coord = stat.x_L + r * stat.incr_x;
		for (int c = 0; c < stat.col_pixels; c++)
		{
			y_coord = stat.y_L + c * stat.incr_y;
			stat.queryVector[r*stat.col_pixels + c][0] = x_coord;
			stat.queryVector[r*stat.col_pixels + c][1] = y_coord;
		}
	}
}

void update_incr_values(statistics& stat)
{
	if (stat.row_pixels == 0 || stat.col_pixels == 0 || stat.t_pixels == 0)
	{
		//cout << "Not valid input dimensions!" << endl;
		//exit(0);
	}

	stat.incr_x = (stat.x_U - stat.x_L) / stat.row_pixels;
	stat.incr_y = (stat.y_U - stat.y_L) / stat.col_pixels;
	stat.incr_t = (stat.t_U - stat.t_L) / stat.t_pixels;
}
