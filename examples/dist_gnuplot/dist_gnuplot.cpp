#include <gnuplot-iostream.h>

#include <ALFI.h>

struct PointLine {
	std::string title;
	std::vector<double> x_points;
};

void plot_data(const std::vector<PointLine>& lines) {
	static const std::vector<std::string> colors = {
		"blue",        // blue        #0000ff =   0   0 255
		"red",         // red         #ff0000 = 255   0   0
		"green",       // green       #00ff00 =   0 255   0
		"orange",      // orange      #ffa500 = 255 165   0
		"purple",      // purple      #c080ff = 192 128 255
		"brown",       // brown       #a52a2a = 165  42  42
		"magenta",     // magenta     #ff00ff = 255   0 255
		"dark-blue",   // dark-blue   #00008b =   0   0 139
		"coral",       // coral       #ff7f50 = 255 127  80
		"dark-violet", // dark-violet #9400d3 = 148   0 211
		"salmon",      // salmon      #fa8072 = 250 128 114
		"dark-red",    // dark-red    #8b0000 = 139   0   0
		"web-blue",    // web-blue    #0080ff =   0 128 255
		"web-green",   // web-green   #00c000 =   0 192   0
		"dark-orange"  // dark-orange #c04000 = 192  64   0
	};

	const double upper_y = 0;
	const double lower_y = -static_cast<double>(lines.size()) - 1;

	Gnuplot gp;

	gp << "set terminal wxt size 1100,600\n";

	gp << "set xrange [0:1]\n";
	gp << "set yrange [" << lower_y << ":" << upper_y << "]\n";

	gp << "plot ";
	for (size_t i = 0; i < lines.size(); ++i) {
		if (i > 0)
			gp << ", ";
		gp << "'-' with points pt 7 ps 1.5 lc rgb '"
			<< colors[i % colors.size()] << "' title '"
			<< lines[i].title << "'";
	}
	gp << '\n';

	double current_layer = upper_y - 1;
	for (const auto& line : lines) {
		for (const auto& x : line.x_points) {
			gp << x << " " << current_layer << "\n";
		}
		gp << "e\n";
		--current_layer;
	}
}

int main() {
	const size_t n = 9;
	const double a = 0;
	const double b = 1;

	const double ratio = 2;
	const double logistic_steepness = 8;
	const double erf_steepness = 3;

	plot_data({
		{"Uniform", alfi::dist::uniform(n, a, b)},
		{"Quadratic", alfi::dist::quadratic(n, a, b)},
		{"Cubic", alfi::dist::cubic(n, a, b)},
		{"Chebyshev", alfi::dist::chebyshev(n, a, b)},
		{"Stretched Chebyshev", alfi::dist::chebyshev_stretched(n, a, b)},
		{"Chebyshev Second Kind", alfi::dist::chebyshev_2(n, a, b)},
		{"Chebyshev Ellipse", alfi::dist::chebyshev_ellipse(n, a, b, ratio)},
		{"Stretched Chebyshev Ellipse", alfi::dist::chebyshev_ellipse_stretched(n, a, b, ratio)},
		{"Chebyshev Ellipse Second Kind", alfi::dist::chebyshev_ellipse_2(n, a, b, ratio)},
		{"Logistic", alfi::dist::logistic(n, a, b, logistic_steepness)},
		{"Stretched Logistic", alfi::dist::logistic_stretched(n, a, b, logistic_steepness)},
		{"Error function", alfi::dist::erf(n, a, b, erf_steepness)},
		{"Stretched Error function", alfi::dist::erf_stretched(n, a, b, erf_steepness)},
	});

	return 0;
}