#pragma once

#include <cmath>
#include <functional>
#include <vector>

#include <ALFI/dist.h>

using FuncName = std::string;
using Lambda = std::function<double(double)>;
using Function = std::pair<FuncName, Lambda>;
using Interval = std::pair<double, double>;

/**
	@brief Collection of functions with their interpolation intervals.

	Each entry consists of:
	- Function: pair (`std::pair`) containing:
	  - Name (`std::string`).
	  - Lambda function (`std::function<double(double)>`).
	- Interval: range (`std::pair<double, double>`).
*/
const std::vector<std::pair<Function, Interval>> funcs_and_ints = {
	{{"exp", [](double x) { return std::exp(x); }}, {-3, 3}},
	{{"sin", [](double x) { return std::sin(x); }}, {-M_PI, M_PI}},
	{{"cos", [](double x) { return std::cos(x); }}, {-M_PI, M_PI}},
	{{"f1", [](double x) { return std::abs(x) + x/2 - x*x; }}, {-1, 1}},
	{{"f2", [](double x) { return -3*std::sin(10*x) + 10*sin(std::abs(x) + x/2); }}, {-10, 10}},
};

const std::vector<alfi::dist::Type> dists = {
	alfi::dist::Type::UNIFORM,
	alfi::dist::Type::CHEBYSHEV,
	alfi::dist::Type::CHEBYSHEV_2,
};

const int64_t nn = 10000;