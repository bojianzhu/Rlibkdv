#include "SLAM.h"

void envelope_point_set(statistics& stat, vector<int>& E_k, SLAM& slam_obj)
{
	for (int i = 0; i < stat.n; i++)
	{
		if (fabs(stat.featureVector[i][stat.static_coord] - slam_obj.k) < stat.bandwidth_s)
			E_k.push_back(i);
	}
}

void bound_list(statistics& stat, vector<int>& E_k, vector<bound_entry>& List, SLAM& slam_obj)
{
	double pm_value;
	int id;

	bound_entry b_entry_ell;
	bound_entry b_entry_u;

	for (int i = 0; i < (int)E_k.size(); i++)
	{
		id = E_k[i];
		pm_value = sqrt(stat.bandwidth_s*stat.bandwidth_s - (slam_obj.k - stat.featureVector[id][stat.static_coord])*(slam_obj.k - stat.featureVector[id][stat.static_coord]));

		b_entry_ell.id = id;
		b_entry_ell.bound_value = stat.featureVector[id][stat.dynamic_coord] - pm_value;
		b_entry_ell.is_LB = true;

		b_entry_u.id = id;
		b_entry_u.bound_value = stat.featureVector[id][stat.dynamic_coord] + pm_value;
		b_entry_u.is_LB = false;

		List.push_back(b_entry_ell);
		List.push_back(b_entry_u);
	}
}

void init_SLAM(statistics& stat)
{
	SLAM temp_SLAM_obj;
	for (int th = 0; th < stat.num_threads; th++)
	{
		stat.SLAM_vec.push_back(temp_SLAM_obj);
		stat.SLAM_vec[th].W_L_ell = 0;
		stat.SLAM_vec[th].W_U_ell = 0;
		stat.SLAM_vec[th].A_L_ell = new double[stat.dim];
		stat.SLAM_vec[th].A_U_ell = new double[stat.dim];
		stat.SLAM_vec[th].S_L_ell = 0;
		stat.SLAM_vec[th].S_U_ell = 0;
		stat.SLAM_vec[th].W_q = 0;
		//stat.SLAM_vec[th].R_q_card = 0;
		stat.SLAM_vec[th].A_R_q = new double[stat.dim];
		stat.SLAM_vec[th].S_R_q = 0;
		stat.SLAM_vec[th].max_KDE = -inf;

		for (int q_id = 0; q_id < stat.dynamic_pixel_size; q_id++)
		{
			double*query = new double[stat.dim];
			stat.SLAM_vec[th].query_list.push_back(query);
			stat.SLAM_vec[th].result_list.push_back(0);
		}

		for (int d = 0; d < stat.dim; d++)
		{
			stat.SLAM_vec[th].A_L_ell[d] = 0;
			stat.SLAM_vec[th].A_U_ell[d] = 0;
			stat.SLAM_vec[th].A_R_q[d] = 0;
		}
	}
}

void clear_SLAM(statistics& stat, SLAM& slam_obj)
{
	for (int d = 0; d < stat.dim; d++)
	{
		slam_obj.A_L_ell[d] = 0;
		slam_obj.A_U_ell[d] = 0;
	}

	slam_obj.S_L_ell = 0;
	slam_obj.S_U_ell = 0;
	slam_obj.W_L_ell = 0;
	slam_obj.W_U_ell = 0;
	//slam_obj.L_ell.clear();
	//slam_obj.U_ell.clear();
}

void SLAM_SORT(statistics& stat, SLAM& slam_obj)
{
	vector<int> E_k;
	vector<bound_entry> List;
	int i_q, i_b;
	double ip;
	double p_SquareNorm;
	bool p_finish = false;

	envelope_point_set(stat, E_k, slam_obj);
	bound_list(stat, E_k, List, slam_obj);

	sort(List.begin(), List.end());
	i_q = 0;
	i_b = 0;

	if (E_k.size() == 0) //degenerate case 1
	{
		for (int i_q = 0; i_q < stat.dynamic_pixel_size; i_q++)
			slam_obj.result_list[i_q] = 0;
		return;
	}

	if (stat.dynamic_pixel_size == 0) //degenerate case 2
		return;

	while (i_q < stat.dynamic_pixel_size)
	{
		if (p_finish == true || slam_obj.query_list[i_q][stat.dynamic_coord] <= List[i_b].bound_value)
		{
			slam_obj.q_SquareNorm = computeSqNorm(slam_obj.query_list[i_q], stat.dim);
			slam_obj.W_q = slam_obj.W_L_ell - slam_obj.W_U_ell;
			//slam_obj.R_q_card = slam_obj.L_ell.size() - slam_obj.U_ell.size();
			for (int d = 0; d < stat.dim; d++)
				slam_obj.A_R_q[d] = slam_obj.A_L_ell[d] - slam_obj.A_U_ell[d];
			slam_obj.S_R_q = slam_obj.S_L_ell - slam_obj.S_U_ell;

			ip = inner_product(slam_obj.query_list[i_q], slam_obj.A_R_q, stat.dim);
			slam_obj.result_list[i_q] = slam_obj.W_q - (1.0 / (stat.bandwidth_s*stat.bandwidth_s))
				*(slam_obj.W_q*slam_obj.q_SquareNorm - 2 * ip + slam_obj.S_R_q);

			i_q++;
		}
		else
		{
			p_SquareNorm = 0;
			if (List[i_b].is_LB == true)
			{
				slam_obj.W_L_ell += stat.weightVector[List[i_b].id];
				//slam_obj.L_ell.push_back(List[i_b].id);
				for (int d = 0; d < stat.dim; d++)
				{
					slam_obj.A_L_ell[d] += stat.weightVector[List[i_b].id] * stat.featureVector[List[i_b].id][d];
					p_SquareNorm += stat.featureVector[List[i_b].id][d] * stat.featureVector[List[i_b].id][d];
				}

				slam_obj.S_L_ell += stat.weightVector[List[i_b].id] * p_SquareNorm;
			}
			else
			{
				slam_obj.W_U_ell += stat.weightVector[List[i_b].id];
				//slam_obj.U_ell.push_back(List[i_b].id);
				for (int d = 0; d < stat.dim; d++)
				{
					slam_obj.A_U_ell[d] += stat.weightVector[List[i_b].id] * stat.featureVector[List[i_b].id][d];
					p_SquareNorm += stat.featureVector[List[i_b].id][d] * stat.featureVector[List[i_b].id][d];
				}

				slam_obj.S_U_ell += stat.weightVector[List[i_b].id] * p_SquareNorm;
			}

			i_b++;
			if (i_b >= 2 * (int)E_k.size())
				p_finish = true;
		}
	}
}

void SLAM_scan_x(statistics& stat, int thread_num)
{
	for (int c = thread_num; c < stat.col_pixels; c = c + stat.num_threads)
	{
		for (int r = 0; r < stat.row_pixels; r++)
		{
			stat.SLAM_vec[thread_num].query_list[r][0] = stat.queryVector[r*stat.col_pixels + c][0];
			stat.SLAM_vec[thread_num].query_list[r][1] = stat.queryVector[r*stat.col_pixels + c][1];
		}

		stat.SLAM_vec[thread_num].k = stat.SLAM_vec[thread_num].query_list[0][1]; //set the parameter k
		SLAM_SORT(stat, stat.SLAM_vec[thread_num]);

		for (int r = 0; r < stat.row_pixels; r++)
		{
			stat.outMatrix[r][c] = stat.SLAM_vec[thread_num].result_list[r];
			if (stat.outMatrix[r][c] > stat.SLAM_vec[thread_num].max_KDE)
				stat.SLAM_vec[thread_num].max_KDE = stat.outMatrix[r][c];
		}

		clear_SLAM(stat, stat.SLAM_vec[thread_num]);
	}
}

void SLAM_visual(statistics& stat)
{
	stat.dynamic_pixel_size = stat.row_pixels;
	stat.static_pixel_size = stat.col_pixels;
	stat.static_coord = 1;
	stat.dynamic_coord = 0;
	thread*th_vec = new thread[stat.num_threads];

	init_SLAM(stat);
	//SLAM_scan_x(stat);

	for (int tid = 0; tid < stat.num_threads; tid++)
		th_vec[tid] = thread(SLAM_scan_x, ref(stat), tid);

	for (int tid = 0; tid < stat.num_threads; tid++)
		th_vec[tid].join();
}