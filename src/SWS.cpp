#include "SWS.h"

void obtain_q(int x_index, int y_index, int t_index, statistics& stat, SWS& sws_obj)
{
	sws_obj.q[0] = stat.x_L + x_index * stat.incr_x;
	sws_obj.q[1] = stat.y_L + y_index * stat.incr_y;
	sws_obj.q[2] = stat.t_L + t_index * stat.incr_t;
}

void init_SWS(statistics& stat)
{
	vector<index_time_pair> pair_vector;
	index_time_pair pair;
	int cur_index;
	SWS temp_SWS_obj;

	stat.sorted_featureVector = new double*[stat.n];
	stat.sorted_weightVector = new double[stat.n];
	for (int i = 0; i < stat.n; i++)
	{
		pair.index = i;
		pair.time = stat.featureVector[i][2];
		pair_vector.push_back(pair);
		stat.sorted_featureVector[i] = new double[stat.dim];
	}

	sort(pair_vector.begin(), pair_vector.end());

	for (int i = 0; i < stat.n; i++)
	{
		cur_index = pair_vector[i].index;
		for (int d = 0; d < stat.dim; d++)
			stat.sorted_featureVector[i][d] = stat.featureVector[cur_index][d];

		stat.sorted_weightVector[i] = stat.weightVector[cur_index];
	}

	//Init for SWS parallelization
	for (int tid = 0; tid < stat.num_threads; tid++)
	{
		stat.SWS_vec.push_back(temp_SWS_obj);
		stat.SWS_vec[tid].max_KDE = -inf;

		stat.SWS_vec[tid].q = new double[stat.dim];
		if (stat.kernel_t_type == 1) //Epanechnikov kernel
			stat.SWS_vec[tid].sliding_window = new double[3];
		if (stat.kernel_t_type == 2) //Quartic kernel
			stat.SWS_vec[tid].sliding_window = new double[5];

		//code here for triangular kernel
	}
}

double compute_init_window_density(statistics& stat, SWS& sws_obj, win_status& win)
{
	win.start_t_win_val = sws_obj.q[2] - stat.bandwidth_t;
	win.end_t_win_val = sws_obj.q[2] + stat.bandwidth_t;
	win.start_t_win_pos = 0;
	double win_density_value;
	double kernel_s_value;
	double t_p_pow;
	bool isStart = false;
	bool isEnd = false;
	//Used in Quartic kernel
	double gamma_t_pow_2, gamma_t_pow_4;
	double t_q_pow_2, t_q_pow_3, t_q_pow_4;

	//initialization
	if (stat.kernel_t_type == 1) //Epanechnikov kernel
	{
		for (int w = 0; w < 3; w++)
			sws_obj.sliding_window[w] = 0;
	}
	if (stat.kernel_t_type == 2) //Quartic kernel
	{
		for (int w = 0; w < 5; w++)
			sws_obj.sliding_window[w] = 0;
	}

	for (int i = 0; i < stat.n; i++)
	{
		if (isStart == false)
		{
			if (stat.sorted_featureVector[i][2] > win.start_t_win_val) //This is the first point that is bigger than the starting point of the interval
			{
				isStart = true;
				win.start_t_win_pos = i;
			}
		}

		if (isEnd == false)
		{
			if (stat.sorted_featureVector[i][2] > win.end_t_win_val)
			{
				isEnd = true;
				win.end_t_win_pos = i - 1;
			}
			else
			{
				if (isStart == true) //This point is inside the interval
				{
					kernel_s_value = spatial_kernel(sws_obj.q, stat.sorted_featureVector[i], stat);
					t_p_pow = 1;

					sws_obj.sliding_window[0] += stat.sorted_weightVector[i] * kernel_s_value;
					if (stat.kernel_t_type == 1) //Epanechnikov kernel
					{
						for (int w = 1; w < 3; w++)
						{
							t_p_pow *= stat.sorted_featureVector[i][2];
							sws_obj.sliding_window[w] += stat.sorted_weightVector[i] * t_p_pow * kernel_s_value;
						}
					}
					if (stat.kernel_t_type == 2) //Quartic kernel
					{
						for (int w = 1; w < 5; w++)
						{
							t_p_pow *= stat.sorted_featureVector[i][2];
							sws_obj.sliding_window[w] += stat.sorted_weightVector[i] * t_p_pow * kernel_s_value;
						}
					}
				}

				if (i == stat.n - 1) //The last point
					win.end_t_win_pos = stat.n - 1;
			}
		}

		if (isStart == true && isEnd == true)
			break;
	}

	if (stat.kernel_t_type == 1) //Epanechnikov kernel
		win_density_value = (1 - (sws_obj.q[2] * sws_obj.q[2]) / (stat.bandwidth_t*stat.bandwidth_t))*sws_obj.sliding_window[0]
		+ 2 * (sws_obj.q[2] * sws_obj.sliding_window[1]) / (stat.bandwidth_t*stat.bandwidth_t)
		- sws_obj.sliding_window[2] / (stat.bandwidth_t*stat.bandwidth_t);

	if (stat.kernel_t_type == 2) //Quartic kernel
	{
		gamma_t_pow_2 = 1 / (stat.bandwidth_t*stat.bandwidth_t);
		gamma_t_pow_4 = gamma_t_pow_2 * gamma_t_pow_2;
		t_q_pow_2 = sws_obj.q[2] * sws_obj.q[2];
		t_q_pow_3 = t_q_pow_2 * sws_obj.q[2];
		t_q_pow_4 = t_q_pow_3 * sws_obj.q[2];
		win_density_value = (1 - 2 * gamma_t_pow_2*t_q_pow_2 + gamma_t_pow_4 * t_q_pow_4)*sws_obj.sliding_window[0]
			+ (4 * gamma_t_pow_2*sws_obj.q[2] - 4 * gamma_t_pow_4*t_q_pow_3)*sws_obj.sliding_window[1]
			+ (6 * gamma_t_pow_4*t_q_pow_2 - 2 * gamma_t_pow_2)*sws_obj.sliding_window[2]
			- 4 * gamma_t_pow_4*sws_obj.q[2] * sws_obj.sliding_window[3] + gamma_t_pow_4 * sws_obj.sliding_window[4];
	}

	return win_density_value;
}

void update_sliding_window(statistics& stat, SWS& sws_obj, vector<int>& index_set, bool is_positive)
{
	double weight;
	double kernel_s_value;
	double t_p_pow;
	int id;

	if (is_positive == true)
		weight = 1;
	else
		weight = -1;

	for (int i = 0; i < (int)index_set.size(); i++)
	{
		id = index_set[i];

		kernel_s_value = spatial_kernel(sws_obj.q, stat.sorted_featureVector[id], stat);
		sws_obj.sliding_window[0] += weight * stat.sorted_weightVector[id] * kernel_s_value;

		if (stat.kernel_t_type == 1) //Epanechnikov kernel
		{
			t_p_pow = 1;
			for (int w = 1; w < 3; w++)
			{
				t_p_pow *= stat.sorted_featureVector[id][2];
				sws_obj.sliding_window[w] += weight * stat.sorted_weightVector[id] * t_p_pow * kernel_s_value;
			}
		}

		if (stat.kernel_t_type == 2) //Quartic kernel
		{
			t_p_pow = 1;
			for (int w = 1; w < 5; w++)
			{
				t_p_pow *= stat.sorted_featureVector[id][2];
				sws_obj.sliding_window[w] += weight * stat.sorted_weightVector[id] * t_p_pow * kernel_s_value;
			}
		}
	}
}

double incr_update_window_density(statistics& stat, SWS& sws_obj, win_status& win)
{
	double win_density_value;
	vector<int> D_set;
	vector<int> I_set;
	int cur_index;
	bool isStart = false;
	bool isEnd = false;
	//Used in Quartic kernel
	double gamma_t_pow_2, gamma_t_pow_4;
	double t_q_pow_2, t_q_pow_3, t_q_pow_4;

	win.start_t_win_val_prev = win.start_t_win_val;
	win.end_t_win_val_prev = win.end_t_win_val;
	win.start_t_win_val = sws_obj.q[2] - stat.bandwidth_t;
	win.end_t_win_val = sws_obj.q[2] + stat.bandwidth_t;

	cur_index = win.start_t_win_pos;
	for (int i = cur_index; i < stat.n; i++)
	{
		if (isStart == false)
		{
			if (stat.sorted_featureVector[i][2] > win.start_t_win_val)
			{
				win.start_t_win_pos = i;
				isStart = true;
			}
		}

		if (isStart == true)
			break;

		if (stat.sorted_featureVector[i][2] <= min(win.end_t_win_val_prev, win.start_t_win_val))
			D_set.push_back(i);
	}

	cur_index = win.end_t_win_pos;
	for (int i = cur_index; i < stat.n; i++)
	{
		if (isEnd == false)
		{
			if (stat.sorted_featureVector[i][2] > win.end_t_win_val)
			{
				win.end_t_win_pos = i - 1;
				isEnd = true;
			}
		}

		if (isEnd == true)
			break;

		if (stat.sorted_featureVector[i][2] > max(win.end_t_win_val_prev, win.start_t_win_val))
			I_set.push_back(i);
	}

	update_sliding_window(stat, sws_obj, D_set, false);
	update_sliding_window(stat, sws_obj, I_set, true);

	if (stat.kernel_t_type == 1) //Epanechnikov kernel
		win_density_value = (1 - (sws_obj.q[2] * sws_obj.q[2]) / (stat.bandwidth_t*stat.bandwidth_t))*sws_obj.sliding_window[0]
		+ 2 * (sws_obj.q[2] * sws_obj.sliding_window[1]) / (stat.bandwidth_t*stat.bandwidth_t)
		- sws_obj.sliding_window[2] / (stat.bandwidth_t*stat.bandwidth_t);

	if (stat.kernel_t_type == 2) //Quartic kernel
	{
		gamma_t_pow_2 = 1 / (stat.bandwidth_t*stat.bandwidth_t);
		gamma_t_pow_4 = gamma_t_pow_2 * gamma_t_pow_2;
		t_q_pow_2 = sws_obj.q[2] * sws_obj.q[2];
		t_q_pow_3 = t_q_pow_2 * sws_obj.q[2];
		t_q_pow_4 = t_q_pow_3 * sws_obj.q[2];
		win_density_value = (1 - 2 * gamma_t_pow_2*t_q_pow_2 + gamma_t_pow_4 * t_q_pow_4)*sws_obj.sliding_window[0]
			+ (4 * gamma_t_pow_2*sws_obj.q[2] - 4 * gamma_t_pow_4*t_q_pow_3)*sws_obj.sliding_window[1]
			+ (6 * gamma_t_pow_4*t_q_pow_2 - 2 * gamma_t_pow_2)*sws_obj.sliding_window[2]
			- 4 * gamma_t_pow_4*sws_obj.q[2] * sws_obj.sliding_window[3] + gamma_t_pow_4 * sws_obj.sliding_window[4];
	}

	return win_density_value;
}

void SWS_algorithm(statistics& stat, int tid)
{
	int num_entries = stat.row_pixels*stat.col_pixels;
	int x_index, y_index;
	win_status win;

	for (int e = tid; e < num_entries; e = e + stat.num_threads)
	{
		x_index = (int)floor((double)e / stat.col_pixels);
		y_index = e - x_index * stat.col_pixels;

		obtain_q(x_index, y_index, 0, stat, stat.SWS_vec[tid]);
		if (stat.kernel_t_type == 1 || stat.kernel_t_type == 2) //Epanechnikov and quartic kernels
		{
			stat.outCube[x_index][y_index][0] = compute_init_window_density(stat, stat.SWS_vec[tid], win);
			stat.SWS_vec[tid].max_KDE = max(stat.SWS_vec[tid].max_KDE, stat.outCube[x_index][y_index][0]);
		}

		//Incremental update algorithm (Section 3.2 in our paper)
		for (int t_index = 1; t_index < stat.t_pixels; t_index++)
		{
			stat.SWS_vec[tid].q[2] = stat.t_L + t_index * stat.incr_t;

			if (stat.kernel_t_type == 1 || stat.kernel_t_type == 2) //Epanechnikov and quartic kernels
			{
				stat.outCube[x_index][y_index][t_index] = incr_update_window_density(stat, stat.SWS_vec[tid], win);
				stat.SWS_vec[tid].max_KDE = max(stat.SWS_vec[tid].max_KDE, stat.outCube[x_index][y_index][t_index]);
			}
		}
	}
}

void SWS_visual(statistics& stat)
{
	thread*th_vec = new thread[stat.num_threads];
	init_SWS(stat);

	for (int tid = 0; tid < stat.num_threads; tid++)
		th_vec[tid] = thread(SWS_algorithm, ref(stat), tid);

	for (int tid = 0; tid < stat.num_threads; tid++)
		th_vec[tid].join();
}