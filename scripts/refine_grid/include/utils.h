#ifndef _UTILS_H_
#define _UTILS_H_

#include <iostream>
#include <cmath>
#include <string>

const std::string LINE_1 = "===========================================================================================================================";
const std::string LINE_2 = "---------------------------------------------------------------------------------------------------------------------------";

inline double calc_norm (double x1, double y1, double z1, double x2, double y2, double z2) {
    return sqrt(pow((x1-x2),2) + pow((y1-y2),2) + pow((z1-z2),2));
}

void usage (const char pName[]);
void print_progress (int iter, int max_iter);

#endif