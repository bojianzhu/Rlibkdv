#pragma once
#ifndef EDWIN_OTF_H
#define EDWIN_OTF_H

#include "init_visual.h"
#include "SWS.h"
#include "bucket.h"

void sort_FeatureVector(statistics& stat);
void init_EDWIN_otf(statistics& stat);
void clear_EDWIN_otf(statistics& stat);
void EDWIN_otf_visual(statistics& stat);

#endif