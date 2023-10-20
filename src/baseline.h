#pragma once
#ifndef BASELINE_H
#define BASELINE_H

#include "init_visual.h"

double SCAN_2D(double*q, statistics& stat);
void SCAN_visual(statistics& stat);
void SCAN_otf_STKDV_visual(statistics& stat);
void SCAN_batch_STKDV_visual(statistics& stat);
double spatial_kernel(double*q, double*p, statistics& stat);
double temporal_kernel(double*q, double*p, statistics& stat);

#endif