#include "bucket.h"

void envelope_point_set(statistics& stat, double k, vector<int>& E_k)
{
	for (int i = stat.start_t_win_pos; i <= stat.end_t_win_pos; i++)
	{
		if (fabs(stat.sorted_featureVector[i][1] - k) < stat.bandwidth_s)
			E_k.push_back(i);
	}
}

void bound_list(statistics& stat, double k, vector<int>& E_k, vector<bound_entry>& List)
{
	double pm_value;
	int id;

	bound_entry b_entry_ell;
	bound_entry b_entry_u;

	for (int i = 0; i < (int)E_k.size(); i++)
	{
		id = E_k[i];
		pm_value = sqrt(stat.bandwidth_s*stat.bandwidth_s - (k - stat.sorted_featureVector[id][1])*(k - stat.sorted_featureVector[id][1]));

		b_entry_ell.id = id;
		b_entry_ell.bound_value = stat.sorted_featureVector[id][0] - pm_value;
		b_entry_ell.is_LB = true;

		b_entry_u.id = id;
		b_entry_u.bound_value = stat.sorted_featureVector[id][0] + pm_value;
		b_entry_u.is_LB = false;

		List.push_back(b_entry_ell);
		List.push_back(b_entry_u);
	}
}

void init_Bucket(statistics& stat)
{
	vector<int> idList;

	stat.C_L = new double[3]; stat.C_U = new double[3]; stat.C_R_q = new double[3];
	stat.v_L = new double*[3]; stat.v_U = new double*[3]; stat.v_R_q = new double*[3];
	stat.H_L = new double[3]; stat.H_U = new double[3]; stat.H_R_q = new double[3];

	for (int u = 0; u <= 2; u++)
	{
		stat.C_L[u] = 0; stat.C_U[u] = 0; stat.C_R_q[u] = 0;
		stat.v_L[u] = new double[2]; stat.v_U[u] = new double[2]; stat.v_R_q[u] = new double[2];
		stat.H_L[u] = 0; stat.H_U[u] = 0; stat.H_R_q[u] = 0;

		for (int d = 0; d < 2; d++)
		{
			stat.v_L[u][d] = 0;
			stat.v_U[u][d] = 0;
			stat.v_R_q[u][d] = 0;
		}
	}

	for (int x_index = 0; x_index <= stat.row_pixels; x_index++)
	{
		stat.B_L.push_back(idList);
		stat.B_U.push_back(idList);
	}
}

void clear_Bucket(statistics& stat)
{
	for (int u = 0; u <= 2; u++)
	{
		stat.C_L[u] = 0; stat.C_U[u] = 0; stat.C_R_q[u] = 0;
		stat.H_L[u] = 0; stat.H_U[u] = 0; stat.H_R_q[u] = 0;

		for (int d = 0; d < 2; d++)
		{
			stat.v_L[u][d] = 0;
			stat.v_U[u][d] = 0;
			stat.v_R_q[u][d] = 0;
		}
	}

	for (int x_index = 0; x_index <= stat.row_pixels; x_index++)
	{
		stat.B_L[x_index].clear();
		stat.B_U[x_index].clear();
	}
}

void erase_Bucket(statistics& stat)
{
	for (int u = 0; u <= 2; u++)
	{
		delete[] stat.v_L[u];
		delete[] stat.v_U[u];
		delete[] stat.v_R_q[u];
	}

	delete[] stat.C_L; delete[] stat.C_U; delete[] stat.C_R_q;
	delete[] stat.H_L; delete[] stat.H_U; delete[] stat.H_R_q;

	stat.B_L.clear();
	stat.B_U.clear();
}

void bucket_algorithm_row(statistics& stat, double k, int y_index, vector<double**>& stat_plane_vec)
{
	vector<int> E_k;
	vector<bound_entry> List;
	int id;
	int i_l, i_u;
	double t_p_u;
	double p_square_norm;
	double bandwidth_square;

	envelope_point_set(stat, k, E_k);
	bound_list(stat, k, E_k, List);

	for (int i = 0; i < (int)List.size(); i++)
	{
		id = List[i].id;
		if (List[i].is_LB == true)
		{
			i_l = (int)max(ceil((List[i].bound_value - stat.x_L) / stat.incr_x), 0.0);
			i_l = (int)min(i_l, stat.row_pixels);
			stat.B_L[i_l].push_back(id);
		}
		if (List[i].is_LB == false)
		{
			i_u = (int)min(ceil((List[i].bound_value - stat.x_L) / stat.incr_x), (double)stat.row_pixels);
			i_u = (int)max(i_u, 0);
			stat.B_U[i_u].push_back(id);
		}
	}

	for (int x_index = 0; x_index < stat.row_pixels; x_index++)
	{
		for (int e = 0; e < (int)stat.B_L[x_index].size(); e++)
		{
			id = stat.B_L[x_index][e];

			t_p_u = 1;
			for (int u = 0; u <= 2; u++)
			{
				stat.C_L[u] += t_p_u;

				p_square_norm = 0;
				for (int d = 0; d < 2; d++)
				{
					stat.v_L[u][d] += t_p_u * stat.sorted_featureVector[id][d];
					p_square_norm += stat.sorted_featureVector[id][d] * stat.sorted_featureVector[id][d];
				}
				stat.H_L[u] += t_p_u * p_square_norm;

				t_p_u *= stat.sorted_featureVector[id][2];
			}
		}

		for (int e = 0; e < (int)stat.B_U[x_index].size(); e++)
		{
			id = stat.B_U[x_index][e];

			t_p_u = 1;
			for (int u = 0; u <= 2; u++)
			{
				stat.C_U[u] += t_p_u;

				p_square_norm = 0;
				for (int d = 0; d < 2; d++)
				{
					stat.v_U[u][d] += t_p_u * stat.sorted_featureVector[id][d];
					p_square_norm += stat.sorted_featureVector[id][d] * stat.sorted_featureVector[id][d];
				}
				stat.H_U[u] += t_p_u * p_square_norm;

				t_p_u *= stat.sorted_featureVector[id][2];
			}
		}

		for (int u = 0; u <= 2; u++)
		{
			stat.C_R_q[u] = stat.C_L[u] - stat.C_U[u];
			stat.H_R_q[u] = stat.H_L[u] - stat.H_U[u];
			for (int d = 0; d < 2; d++)
				stat.v_R_q[u][d] = stat.v_L[u][d] - stat.v_U[u][d];

			bandwidth_square = stat.bandwidth_s*stat.bandwidth_s;
			//Code here (Compute S_I^{(u)}(q))
			stat.q[0] = stat.x_L + x_index * stat.incr_x;
			stat.q[1] = k;

			stat_plane_vec[u][x_index][y_index] = (1.0 - inner_product(stat.q, stat.q, 2) / bandwidth_square)*stat.C_R_q[u]
				+ (2.0 / bandwidth_square) * inner_product(stat.q, stat.v_R_q[u], 2) - stat.H_R_q[u] / bandwidth_square;
		}
	}
}

void bucket_algorithm(statistics& stat, vector<double**>& stat_plane_vec)
{
	//init_Bucket(stat);
	for (int y_index = 0; y_index < stat.col_pixels; y_index++)
	{
		bucket_algorithm_row(stat, stat.y_L + stat.incr_y*y_index, y_index, stat_plane_vec);
		clear_Bucket(stat);
	}

	//erase_Bucket(stat);
}