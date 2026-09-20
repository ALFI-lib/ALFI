#pragma once

#include <iostream>
#include <cmath>

#include "../config.h"

namespace alfi::poly {
	template <typename Number = DefaultNumber, template <typename, typename...> class Container = DefaultContainer>
	Container<Number> lagrange_vals(const Container<Number>& X, const Container<Number>& Y, const Container<Number>& xx) {
		const auto nn = xx.size();

		if (X.size() != Y.size()) {
			std::cerr << "Error in function " << __FUNCTION__
					  << ": Vectors X (of size " << X.size()
					  << ") and Y (of size " << Y.size()
					  << ") are not the same size. Returning an array of NaNs..." << std::endl;
			Container<Number> yy(nn);
			std::fill(yy.begin(), yy.end(), NAN);
			return yy;
		}

		if (X.empty()) {
			std::cerr << "Error in function " << __FUNCTION__
					  << ": Vectors X and Y are empty. Cannot interpolate. Returning an array of NaNs..." << std::endl;
			Container<Number> yy(nn);
			std::fill(yy.begin(), yy.end(), NAN);
			return yy;
		}

		const auto N = X.size();

		Container<Number> yy(nn);

#if defined(_OPENMP) && !defined(ALFI_DISABLE_OPENMP)
#pragma omp parallel for
#endif
		for (SizeT k = 0; k < nn; ++k) {
			yy[k] = 0;
			for (SizeT i = 0; i < N; ++i) {
				Number l = 1;
				for (SizeT j = 0; j < N; ++j) {
					if (i != j) {
						l *= (xx[k] - X[j]) / (X[i] - X[j]);
					}
				}
				yy[k] += Y[i] * l;
			}
		}

		return yy;
	}
}