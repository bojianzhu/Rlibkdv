#include "EDWIN_multiple.h"

void init_EDWIN_multiple(statistics& stat)
{
	double**temp_Plane;

	for (int u = 0; u < 3; u++)
	{
	  temp_Plane = new double*[stat.row_pixels];
		stat.S_plane_vec.push_back(temp_Plane);
		temp_Plane = new double*[stat.row_pixels];
		stat.S_D_plane_vec.push_back(temp_Plane);
		temp_Plane = new double*[stat.row_pixels];
		stat.S_I_plane_vec.push_back(temp_Plane);

		//stat.S_D_plane_vec[u] = new double*[stat.row_pixels];
    //stat.S_I_plane_vec[u] = new double*[stat.row_pixels];
		//stat.S_plane_vec[u] = new double*[stat.row_pixels];
		for (int x_index = 0; x_index < stat.row_pixels; x_index++)
		{
			stat.S_D_plane_vec[u][x_index] = new double[stat.col_pixels];
			stat.S_I_plane_vec[u][x_index] = new double[stat.col_pixels];
			stat.S_plane_vec[u][x_index] = new double[stat.col_pixels];
		}
	}

	stat.q = new double[2];
	sort_FeatureVector(stat);
	for (int i = 0; i < stat.n; i++)
		stat.sorted_fV_timestamp_vec.push_back(stat.sorted_featureVector[i][2]);
}

void clear_EDWIN_multiple(statistics& stat)
{
	for (int u = 0; u <= 2; u++)
	{
		for (int x_index = 0; x_index < stat.row_pixels; x_index++)
		{
			for (int y_index = 0; y_index < stat.col_pixels; y_index++)
			{
				stat.S_D_plane_vec[u][x_index][y_index] = 0;
				stat.S_I_plane_vec[u][x_index][y_index] = 0;
			}
		}
	}
}

void EDWIN_multiple(statistics& stat)
{
	vector<double>::iterator left, right;
	double coeff_S_W0, coeff_S_W1, coeff_S_W2;
	double cur_time;
	double D_largest_time;
	double I_largest_time;
	stat.max_EDWIN_KDE = -inf;
	//int prev_start_t_win;
	//int prev_end_t_win;
	//for (int i = 0; i < stat.n; i++)
	//	cout << stat.sorted_featureVector[i][0] << " " << stat.sorted_featureVector[i][1] << " " << stat.sorted_featureVector[i][2] << endl;

	init_EDWIN_multiple(stat);

	init_Bucket(stat);
	for (int t_index = 0; t_index < stat.t_pixels; t_index++)
	{
		cur_time = stat.t_L + t_index * stat.incr_t;

		if (t_index > 0)
		{
			D_largest_time = cur_time - stat.bandwidth_t;
			I_largest_time = cur_time + stat.bandwidth_t;

			stat.start_S_D_t_win_pos = stat.start_t_win_pos;
			stat.end_S_D_t_win_pos = stat.start_t_win_pos - 1;
			stat.start_S_I_t_win_pos = min(stat.end_t_win_pos + 1, stat.n - 1);
			stat.end_S_I_t_win_pos = stat.start_S_I_t_win_pos - 1;

			for (int i = stat.start_S_D_t_win_pos; i < stat.n; i++)
			{
				if (stat.sorted_fV_timestamp_vec[i] < D_largest_time)
					stat.end_S_D_t_win_pos = i;
				else
					break;
			}

			for (int i = stat.start_S_I_t_win_pos; i < stat.n; i++)
			{
				if (stat.sorted_fV_timestamp_vec[i] < I_largest_time)
					stat.end_S_I_t_win_pos = i;
				else
					break;
			}

			stat.start_t_win_pos = stat.start_S_D_t_win_pos; stat.end_t_win_pos = stat.end_S_D_t_win_pos;
			bucket_algorithm(stat, stat.S_D_plane_vec);
			stat.start_t_win_pos = stat.start_S_I_t_win_pos; stat.end_t_win_pos = stat.end_S_I_t_win_pos;
			bucket_algorithm(stat, stat.S_I_plane_vec);

			for (int u = 0; u <= 2; u++)
				for (int x_index = 0; x_index < stat.row_pixels; x_index++)
					for (int y_index = 0; y_index < stat.col_pixels; y_index++)
						stat.S_plane_vec[u][x_index][y_index] += (stat.S_I_plane_vec[u][x_index][y_index] - stat.S_D_plane_vec[u][x_index][y_index]);


			stat.start_t_win_pos = stat.end_S_D_t_win_pos + 1;
			stat.end_t_win_pos = stat.end_S_I_t_win_pos;
		}
		else
		{
			//Find the lower and upper bound of the time window
			left = lower_bound(stat.sorted_fV_timestamp_vec.begin(),
				stat.sorted_fV_timestamp_vec.end(), cur_time - stat.bandwidth_t);
			right = upper_bound(stat.sorted_fV_timestamp_vec.begin(),
				stat.sorted_fV_timestamp_vec.end(), cur_time + stat.bandwidth_t);

			stat.start_t_win_pos = left - stat.sorted_fV_timestamp_vec.begin();
			stat.end_t_win_pos = min((int)(right - stat.sorted_fV_timestamp_vec.begin()), stat.n) - 1;

			bucket_algorithm(stat, stat.S_plane_vec);
		}

		coeff_S_W0 = 1.0 - (cur_time * cur_time) / (stat.bandwidth_t*stat.bandwidth_t);
		coeff_S_W1 = 2.0 * cur_time / (stat.bandwidth_t*stat.bandwidth_t);
		coeff_S_W2 = -1.0 / (stat.bandwidth_t*stat.bandwidth_t);

		for (int x_index = 0; x_index < stat.row_pixels; x_index++)
		{
			for (int y_index = 0; y_index < stat.col_pixels; y_index++)
			{
				stat.outCube[x_index][y_index][t_index] = coeff_S_W0 * stat.S_plane_vec[0][x_index][y_index]
					+ coeff_S_W1 * stat.S_plane_vec[1][x_index][y_index]
					+ coeff_S_W2 * stat.S_plane_vec[2][x_index][y_index];

				stat.max_EDWIN_KDE = max(stat.outCube[x_index][y_index][t_index], stat.max_EDWIN_KDE);
			}
		}

		clear_EDWIN_multiple(stat);
	}

	erase_Bucket(stat);
}
