#include "EDWIN_otf.h"

void sort_FeatureVector(statistics& stat)
{
	vector<index_time_pair> pair_vector;
	index_time_pair pair;
	int cur_index;

	stat.sorted_featureVector = new double*[stat.n];
	for (int i = 0; i < stat.n; i++)
	{
		pair.index = i;
		pair.time = stat.featureVector[i][2];
		pair_vector.push_back(pair);
		stat.sorted_featureVector[i] = new double[3];
	}

	sort(pair_vector.begin(), pair_vector.end());

	for (int i = 0; i < stat.n; i++)
	{
		cur_index = pair_vector[i].index;
		for (int d = 0; d < 3; d++)
			stat.sorted_featureVector[i][d] = stat.featureVector[cur_index][d];
	}
}

void init_EDWIN_otf(statistics& stat)
{
	double**temp_Plane;
	for (int u = 0; u <= 2; u++)
	{
	  temp_Plane = new double*[stat.row_pixels];
		stat.S_plane_vec.push_back(temp_Plane);
		//stat.S_plane_vec[u] = new double*[stat.row_pixels];
		for (int x_index = 0; x_index < stat.row_pixels; x_index++)
		{
			stat.S_plane_vec[u][x_index] = new double[stat.col_pixels];
			for (int y_index = 0; y_index < stat.col_pixels; y_index++)
				stat.S_plane_vec[u][x_index][y_index] = 0;
		}
	}

	stat.q = new double[2];
	sort_FeatureVector(stat);
	for (int i = 0; i < stat.n; i++)
		stat.sorted_fV_timestamp_vec.push_back(stat.sorted_featureVector[i][2]);
}

void clear_EDWIN_otf(statistics& stat)
{
	for (int u = 0; u <= 2; u++)
		for (int x_index = 0; x_index < stat.row_pixels; x_index++)
			for (int y_index = 0; y_index < stat.col_pixels; y_index++)
				stat.S_plane_vec[u][x_index][y_index] = 0;
}

void EDWIN_otf_visual(statistics& stat)
{
	vector<double>::iterator left, right;
	double coeff_S_W0, coeff_S_W1, coeff_S_W2;
	stat.max_EDWIN_KDE = -inf;

	init_EDWIN_otf(stat);

	//Find the lower and upper bound of the time window
	left = lower_bound(stat.sorted_fV_timestamp_vec.begin(),
		stat.sorted_fV_timestamp_vec.end(), stat.cur_time - stat.bandwidth_t);
	right = upper_bound(stat.sorted_fV_timestamp_vec.begin(),
		stat.sorted_fV_timestamp_vec.end(), stat.cur_time + stat.bandwidth_t);

	stat.start_t_win_pos = left - stat.sorted_fV_timestamp_vec.begin();
	stat.end_t_win_pos = min((int)(right - stat.sorted_fV_timestamp_vec.begin()), stat.n) - 1;

	coeff_S_W0 = 1.0 - (stat.cur_time * stat.cur_time) / (stat.bandwidth_t*stat.bandwidth_t);
	coeff_S_W1 = 2.0 * stat.cur_time / (stat.bandwidth_t*stat.bandwidth_t);
	coeff_S_W2 = -1.0 / (stat.bandwidth_t*stat.bandwidth_t);

	init_Bucket(stat);
	bucket_algorithm(stat, stat.S_plane_vec);
	erase_Bucket(stat);

	for (int x_index = 0; x_index < stat.row_pixels; x_index++)
	{
		for (int y_index = 0; y_index < stat.col_pixels; y_index++)
		{
			stat.outMatrix[x_index][y_index] = coeff_S_W0 * stat.S_plane_vec[0][x_index][y_index]
				+ coeff_S_W1 * stat.S_plane_vec[1][x_index][y_index]
				+ coeff_S_W2 * stat.S_plane_vec[2][x_index][y_index];

			stat.max_EDWIN_KDE = max(stat.max_EDWIN_KDE, stat.outMatrix[x_index][y_index]);
		}
	}


	//clear_EDWIN_otf(stat);
}
