#include "alg_visual.h"

void alg_visual::load_parameters(int argc, char**argv)
{
	//Testing Mode (KDV_type = 1)
	/*stat.num_threads = 8;
	stat.x_L = 1000;
	stat.x_U = 30000;
	stat.y_L = 1000;
	stat.y_U = 50000;
	stat.outputFileName = (char*)"./Results/Atlanta_SLAM_Testing_T8_960_800";
	stat.row_pixels = 960;
	stat.col_pixels = 800;
	stat.kernel_s_type = 1;
	stat.bandwidth_s = 1000;*/

	/*stat.num_threads = 4;
	stat.row_pixels = 128;
	stat.col_pixels = 128;
	stat.kernel_s_type = 1;
	stat.outputFileName = (char*)"./Results/Atlanta_SLAM_T4_v2";
	stat.bandwidth_s = 1000;
	obtain_L_U();*/

	//Testing Mode (KDV_type = 2)
	/*stat.outputFileName = (char*)"./Results/data_ST_SCAN_otf";
	stat.num_threads = 1;
	stat.x_L = 260000;
	stat.x_U = 300000;
	stat.y_L = 10000;
	stat.y_U = 30000;
	stat.row_pixels = 32;
	stat.col_pixels = 32;
	stat.kernel_s_type = 1;
	stat.bandwidth_s = 1000;*/
	//obtain_L_U();

	//Testing Mode (KDV_type = 3)
	/*stat.outputFileName = (char*)"./Results/case_ST_EDWIN_kt";
	stat.num_threads = 1;
	stat.x_L = 113.5222;
	stat.x_U = 113.6022;
	stat.y_L = 22.1036;
	stat.y_U = 22.2195;
	stat.row_pixels = 256;
	stat.col_pixels = 256;
	stat.kernel_s_type = 1;
	stat.bandwidth_s = 0.0001;*/

	stat.num_threads = atoi(argv[3]);
	stat.x_L = atof(argv[4]);
	stat.x_U = atof(argv[5]);
	stat.y_L = atof(argv[6]);
	stat.y_U = atof(argv[7]);
	stat.row_pixels = atoi(argv[8]);
	stat.col_pixels = atoi(argv[9]);
	stat.kernel_s_type = atoi(argv[10]);
	stat.bandwidth_s = atof(argv[11]);

	if (stat.KDV_type == 1) //KDV
		stat.dim = 2;

	if (stat.KDV_type == 2) //online STKDV
	{
		//Testing Mode (KDV_type = 2)
		/*stat.t_L = 100;
		stat.t_U = 2800;
		stat.kernel_t_type = 1;
		stat.bandwidth_t = 14;
		stat.cur_time = 2500;*/

		stat.dim = 3;
		stat.t_L = atof(argv[12]);
		stat.t_U = atof(argv[13]);
		stat.kernel_t_type = atoi(argv[14]);
		stat.bandwidth_t = atof(argv[15]);
		stat.cur_time = atof(argv[16]);
	}
	if (stat.KDV_type == 3) //batch-based STKDV
	{
		//Testing Mode (KDV_type = 3)
		/*stat.t_L = 1;
		stat.t_U = 30;
		stat.t_pixels = 30;
		stat.kernel_t_type = 1;
		stat.bandwidth_t = 1;*/

		stat.dim = 3;
		stat.t_L = atof(argv[12]);
		stat.t_U = atof(argv[13]);
		stat.t_pixels = atoi(argv[14]);
		stat.kernel_t_type = atoi(argv[15]);
		stat.bandwidth_t = atof(argv[16]);
	}
}

void alg_visual::filter_datasets()
{
	int ori_n = stat.base_dataMatrix.size();
	int n = 0;
	double x_value, y_value, weight;

	for (int i = 0; i < ori_n; i++)
	{
		x_value = stat.base_dataMatrix[i][0];
		y_value = stat.base_dataMatrix[i][1];
		weight = stat.base_weightVector[i];

		//filter those data points that are out of range
		if (x_value < stat.x_L - stat.bandwidth_s || x_value > stat.x_U + stat.bandwidth_s
			|| y_value < stat.y_L - stat.bandwidth_s || y_value > stat.y_U + stat.bandwidth_s)
			continue;

		stat.featureVector.push_back(new double[stat.dim]);
		stat.weightVector.push_back(weight);
		stat.featureVector[n][0] = x_value;
		stat.featureVector[n][1] = y_value;

		if (stat.KDV_type == 2 || stat.KDV_type == 3)
			stat.featureVector[n][2] = stat.base_dataMatrix[i][2];

		n++;
	}

	stat.n = n;
}

void alg_visual::init_visual()
{
	if (stat.KDV_type == 1 || stat.KDV_type == 2) //KDV or online STKDV
	{
		initQuery(stat);

		stat.outMatrix = new double*[stat.row_pixels];
		for (int r = 0; r < stat.row_pixels; r++)
			stat.outMatrix[r] = new double[stat.col_pixels];

		//do the preprocessing here (if any)
	}

	if (stat.KDV_type == 3) //batch-based STKDV
	{
		update_incr_values(stat);
		stat.outCube = new double**[stat.row_pixels];
		for (int r = 0; r < stat.row_pixels; r++)
			stat.outCube[r] = new double*[stat.col_pixels];

		for (int r = 0; r < stat.row_pixels; r++)
			for (int c = 0; c < stat.col_pixels; c++)
				stat.outCube[r][c] = new double[stat.t_pixels];

		//do the preprocessing here (if any)
	}
}

void alg_visual::visual_Algorithm()
{
	if (stat.KDV_type == 1) //KDV (Epanechnikov kernel only)
		SLAM_visual(stat);
		//SCAN_visual(stat);

	if (stat.KDV_type == 2) //STKDV (Epanechnikov kernel only)
		EDWIN_otf_visual(stat);

	if (stat.KDV_type == 3) //STKDV (Epanechnikov kernel only)
		EDWIN_multiple(stat);
		//SWS_visual(stat); //STKDV (Epanechnikov and Quartic kernels only)
		//SCAN_batch_STKDV_visual(stat);
}

void alg_visual::matrix_normalization(double max_KDE)
{
	for (int r = 0; r < stat.row_pixels; r++)
		for (int c = 0; c < stat.col_pixels; c++)
			stat.outMatrix[r][c] *= 255.0 / max_KDE;
}

void alg_visual::cube_normalization(double max_KDE)
{
	for (int r = 0; r < stat.row_pixels; r++)
		for (int c = 0; c < stat.col_pixels; c++)
			for (int t = 0; t < stat.t_pixels; t++)
				stat.outCube[r][c][t] *= 255.0 / max_KDE;
}

string alg_visual::saveMatrix_toString_CSV()
{
	double max_KDE = -inf;
	double x, y;
	stringstream outString_ss;

	if (stat.KDV_type == 1)
	{
		for (int tid = 0; tid < stat.num_threads; tid++)
			max_KDE = max(max_KDE, stat.SLAM_vec[tid].max_KDE);
	}
	else //KDV_type = 2
		max_KDE = stat.max_EDWIN_KDE;

	matrix_normalization(max_KDE);

	for (int r = 0; r < stat.row_pixels; r++)
	{
		for (int c = 0; c < stat.col_pixels; c++)
		{
			if (stat.outMatrix[r][c] < small_epsilon)
				continue;

			x = stat.queryVector[r*stat.col_pixels + c][0];
			y = stat.queryVector[r*stat.col_pixels + c][1];

			outString_ss << setprecision(10) << x << "," << y << "," << stat.outMatrix[r][c] << endl;
		}
	}

	clear_memory();
	return outString_ss.str();
}

string alg_visual::saveCube_toString_CSV()
{
	double max_KDE = -inf;
	double x, y, time;
	stringstream outString_ss;

	//for (int tid = 0; tid < stat.num_threads; tid++)
	//	max_KDE = max(max_KDE, stat.SWS_vec[tid].max_KDE);
	max_KDE = stat.max_EDWIN_KDE;

	cube_normalization(max_KDE);

	for (int r = 0; r < stat.row_pixels; r++)
	{
		x = stat.x_L + r * stat.incr_x;
		for (int c = 0; c < stat.col_pixels; c++)
		{
			y = stat.y_L + c * stat.incr_y;
			for (int t = 0; t < stat.t_pixels; t++)
			{
				time = stat.t_L + t * stat.incr_t;
				if (stat.outCube[r][c][t] < small_epsilon)
					continue;

				outString_ss << setprecision(10) << x << "," << y << "," << time << "," << stat.outCube[r][c][t] << endl;
			}
		}
	}

	// cout << outString_ss.str() << endl;
	clear_memory();
	return outString_ss.str();
}

string alg_visual::saveMatrix_toString()
{
	int counter = 0;
	double max_KDE = -inf;
	double x, y;
	stringstream outString_ss;

	if (stat.KDV_type == 1)
	{
		for (int tid = 0; tid < stat.num_threads; tid++)
			max_KDE = max(max_KDE, stat.SLAM_vec[tid].max_KDE);
	}
	else //KDV_type = 2
		max_KDE = stat.max_EDWIN_KDE;

	matrix_normalization(max_KDE);

	outString_ss << "[";

	for (int r = 0; r < stat.row_pixels; r++)
	{
		for (int c = 0; c < stat.col_pixels; c++)
		{
			if (stat.outMatrix[r][c] < small_epsilon)
				continue;
			else
			{
				if (counter != 0)
					outString_ss << "," << endl;
				counter++;
			}

			x = stat.queryVector[r*stat.col_pixels + c][0];
			y = stat.queryVector[r*stat.col_pixels + c][1];

			outString_ss << setprecision(10) << "{\"x\": " << x << ", \"y\": " << y << ", \"value\": " << stat.outMatrix[r][c] << "}";
		}
	}
	outString_ss << "]";

	clear_memory();
	return outString_ss.str();
}

string alg_visual::saveCube_toString()
{
	int counter = 0;
	double max_KDE = -inf;
	double x, y, time;
	stringstream outString_ss;

	//for (int tid = 0; tid < stat.num_threads; tid++)
	//	max_KDE = max(max_KDE, stat.SWS_vec[tid].max_KDE);
	max_KDE = stat.max_EDWIN_KDE;

	cube_normalization(max_KDE);

	outString_ss << "[";
	for (int r = 0; r < stat.row_pixels; r++)
	{
		x = stat.x_L + r * stat.incr_x;
		for (int c = 0; c < stat.col_pixels; c++)
		{
			y = stat.y_L + c * stat.incr_y;
			for (int t = 0; t < stat.t_pixels; t++)
			{
				time = stat.t_L + t * stat.incr_t;
				if (stat.outCube[r][c][t] < small_epsilon)
					continue;
				else
				{
					if (counter != 0)
						outString_ss << "," << endl;
					counter++;
				}

				outString_ss << setprecision(10) << "{\"x\": " << x << ", \"y\": " << y << ", \"time\": " << time << ", \"value\": " << stat.outCube[r][c][t] << "}";
			}
		}
	}
	outString_ss << "]";

	clear_memory();
	return outString_ss.str();
}

string alg_visual::compute(int argc, char**argv)
{
	string outString;
	load_parameters(argc, argv);
	filter_datasets();
	init_visual();
	visual_Algorithm();
	//output_File(); //debug mode
	//exit(0); //debug mode

	if (stat.KDV_type == 1 || stat.KDV_type == 2)
		return saveMatrix_toString_CSV();
	if (stat.KDV_type == 3)
		return saveCube_toString_CSV();

	return "";
}

void alg_visual::clear_memory()
{
	int total_q = stat.row_pixels*stat.col_pixels;
	int ori_n = stat.base_dataMatrix.size();

	for (int i = 0; i < ori_n; i++)
		delete[] stat.featureVector[i];
	stat.featureVector.clear();
	stat.weightVector.clear();

	if (stat.KDV_type == 1 || stat.KDV_type == 2)
	{
		for (int q = 0; q < total_q; q++)
			delete[] stat.queryVector[q];
		delete[] stat.queryVector;

		for (int r = 0; r < stat.row_pixels; r++)
			delete[] stat.outMatrix[r];
		delete[] stat.outMatrix;
	}

	if (stat.KDV_type == 1) //KDV
	{
		for (int th = 0; th < stat.num_threads; th++)
		{
			delete[] stat.SLAM_vec[th].A_L_ell;
			delete[] stat.SLAM_vec[th].A_U_ell;
			delete[] stat.SLAM_vec[th].A_R_q;

			for (int q_id = 0; q_id < stat.dynamic_pixel_size; q_id++)
				delete[] stat.SLAM_vec[th].query_list[q_id];

			stat.SLAM_vec[th].query_list.clear();
			stat.SLAM_vec[th].result_list.clear();
		}

		stat.SLAM_vec.clear();
	}

	if (stat.KDV_type == 2) //Exploratory STKDV
	{
		for (int u = 0; u <= 2; u++)
			for (int x_index = 0; x_index < stat.row_pixels; x_index++)
				delete[] stat.S_plane_vec[u][x_index];

		for (int u = 0; u <= 2; u++)
			delete[] stat.S_plane_vec[u];

		stat.S_plane_vec.clear();

		delete[] stat.q;

		for (int i = 0; i < stat.n; i++)
			delete[] stat.sorted_featureVector[i];

		delete[] stat.sorted_featureVector;
		stat.sorted_fV_timestamp_vec.clear();
	}

	if (stat.KDV_type == 3) //Batch-based STKDV
	{
		for (int r = 0; r < stat.row_pixels; r++)
			for (int c = 0; c < stat.col_pixels; c++)
				delete[] stat.outCube[r][c];

		for (int r = 0; r < stat.row_pixels; r++)
			delete[] stat.outCube[r];

		delete[] stat.outCube;

		for (int i = 0; i < stat.n; i++)
			delete[] stat.sorted_featureVector[i];
		delete[] stat.sorted_featureVector;

		//the EDWIN_multiple method
		delete[] stat.q;
		for (int u = 0; u < 3; u++)
		{
			for (int x_index = 0; x_index < stat.row_pixels; x_index++)
			{
				delete[] stat.S_plane_vec[u][x_index];
				delete[] stat.S_D_plane_vec[u][x_index];
				delete[] stat.S_I_plane_vec[u][x_index];
			}
		}

		for (int u = 0; u < 3; u++)
		{
			delete[] stat.S_plane_vec[u];
			delete[] stat.S_D_plane_vec[u];
			delete[] stat.S_I_plane_vec[u];
		}

		stat.S_plane_vec.clear();
		stat.S_D_plane_vec.clear();
		stat.S_I_plane_vec.clear();
		stat.sorted_fV_timestamp_vec.clear();

		//SWS method
		/*for (int tid = 0; tid < stat.num_threads; tid++)
		{
			delete[] stat.SWS_vec[tid].q;
			delete[] stat.SWS_vec[tid].sliding_window;
		}

		stat.SWS_vec.clear();*/
	}
}

void alg_visual::clear_basic_memory()
{
	int ori_n = stat.base_dataMatrix.size();
	for (int i = 0; i < ori_n; i++)
		delete[] stat.base_dataMatrix[i];

	stat.base_dataMatrix.clear();
	stat.base_weightVector.clear();
}

void alg_visual::obtain_L_U()
{
	stat.x_L = inf; stat.y_L = inf; //Used for testing
	stat.x_U = -inf; stat.y_U = -inf; //Used for testing

	if (stat.KDV_type == 2 || stat.KDV_type == 3)
	{
		stat.t_L = inf;
		stat.t_U = -inf;
	}

	int ori_n = stat.base_dataMatrix.size();

	for (int i = 0; i < ori_n; i++)
	{
		if (stat.base_dataMatrix[i][0] < stat.x_L)
			stat.x_L = stat.base_dataMatrix[i][0];
		if (stat.base_dataMatrix[i][0] > stat.x_U)
			stat.x_U = stat.base_dataMatrix[i][0];
		if (stat.base_dataMatrix[i][1] < stat.y_L)
			stat.y_L = stat.base_dataMatrix[i][1];
		if (stat.base_dataMatrix[i][1] > stat.y_U)
			stat.y_U = stat.base_dataMatrix[i][1];

		if (stat.KDV_type == 2 || stat.KDV_type == 3)
		{
			if (stat.base_dataMatrix[i][2] < stat.t_L)
				stat.t_L = stat.base_dataMatrix[i][2];
			if (stat.base_dataMatrix[i][2] > stat.t_U)
				stat.t_U = stat.base_dataMatrix[i][2];
		}
	}
}

void alg_visual::output_File()
{
	fstream outFile;
	double x, y, time;

	outFile.open(stat.outputFileName, ios::in | ios::out | ios::trunc);
	if (outFile.is_open() == false)
		//cout << "Cannot open output file!" << endl;

	for (int r = 0; r < stat.row_pixels; r++)
	{
		x = stat.x_L + r * stat.incr_x;
		for (int c = 0; c < stat.col_pixels; c++)
		{
			y = stat.y_L + c * stat.incr_y;

			if (stat.KDV_type == 1 || stat.KDV_type == 2)
				outFile << stat.outMatrix[r][c] << endl;
				//outFile << x << " " << y << " " << stat.outMatrix[r][c] << endl;

			if (stat.KDV_type == 3)
			{
				for (int t = 0; t < stat.t_pixels; t++)
				{
					time = stat.t_L + t * stat.incr_t;
					outFile << stat.outCube[r][c][t] << endl;
					//outFile << x << " " << y << " " << time << " " << stat.outCube[r][c][t] << endl;
				}
			}
		}
	}

	outFile.close();
}

void alg_visual::load_datasets(char**argv)
{
	fstream dataFile_json;
	string lineString;
	char* token;
	int ori_n = 0;
	double x, y, t, weight;

	//stat.dataFileName_JSON = (char*)"../../../datasets/Atlanta/Atlanta_spatial_shift.json";
	//stat.KDV_type = 1;

	stat.dataFileName_JSON = argv[1];
	stat.KDV_type = atoi(argv[2]);

	dataFile_json.open(stat.dataFileName_JSON);
	if (dataFile_json.is_open() == false)
	{
		//cout << "Cannot Open File!" << endl;
		//exit(1);
	}

	while (getline(dataFile_json, lineString))
	{
		if (lineString == "")
			break;

		token = strtok((char*)lineString.c_str(), " :,}"); token = strtok(NULL, " :,}");
		x = atof(token);
		token = strtok(NULL, " :,}"); token = strtok(NULL, " :,}");
		y = atof(token);
		stat.base_dataMatrix.push_back(new double[3]);
		stat.base_dataMatrix[ori_n][0] = x;
		stat.base_dataMatrix[ori_n][1] = y;

		if (stat.KDV_type == 2 || stat.KDV_type == 3) //Online STKDV or batch-based STKDV
		{
			token = strtok(NULL, " :,}"); token = strtok(NULL, " :,}");
			t = atof(token);
			stat.base_dataMatrix[ori_n][2] = t;
		}

		token = strtok(NULL, " :,}"); token = strtok(NULL, " :,}");
		weight = atof(token);
		stat.base_weightVector.push_back(weight);

		ori_n++;
	}

	dataFile_json.close();
}

void alg_visual::load_datasets_CSV(char**argv)
{
	fstream dataFile_CSV;
	string lineString;
	char* token;
	int ori_n = 0;
	double x, y, t, weight;

	//Testing Mode
	//stat.dataFileName_CSV = (char*)"../../../Datasets/cases.csv";
	//stat.KDV_type = 3;

	string data_str = argv[1];
	stat.KDV_type = atoi(argv[2]);

	// dataFile_CSV.open(stat.dataFileName_CSV, ios::in | ios::out);
	// if (dataFile_CSV.is_open() == false)
	// {
	// 	cout << "Cannot Open File!" << endl;
	// 	exit(1);
	// }
    istringstream iss(data_str);
	getline(iss, lineString);
	while (getline(iss, lineString))
	{
		if (lineString == "")
			break;

		token = strtok((char*)lineString.c_str(), " ,");
		x = atof(token);
		token = strtok(NULL, " ,");
		y = atof(token);
		stat.base_dataMatrix.push_back(new double[3]);
		stat.base_dataMatrix[ori_n][0] = x;
		stat.base_dataMatrix[ori_n][1] = y;

		if (stat.KDV_type == 2 || stat.KDV_type == 3) //Online STKDV or batch-based STKDV
		{
			token = strtok(NULL, " ,");
			t = atof(token);
			stat.base_dataMatrix[ori_n][2] = t;
		}

		//token = strtok(NULL, " ,");
		//weight = atof(token);
		weight = 1;
		stat.base_weightVector.push_back(weight);

		ori_n++;
	}

	dataFile_CSV.close();
}
