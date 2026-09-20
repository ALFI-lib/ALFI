#pragma once

#include <iostream>
#include <cmath>

#include "../config.h"

namespace alfi::poly {
	template <typename Number = DefaultNumber, template <typename, typename...> class Container = DefaultContainer>
	Container<Number> lagrange(const Container<Number>& X, const Container<Number>& Y) {
		if (X.size() != Y.size()) {
			std::cerr << "Error in function " << __FUNCTION__
					  << ": Vectors X (of size " << X.size()
					  << ") and Y (of size " << Y.size()
					  << ") are not the same size. Returning {NAN}..." << std::endl;
			return {NAN};
		}

		if (X.empty()) {
			std::cerr << "Error in function " << __FUNCTION__
					  << ": Vectors X and Y are empty. Cannot interpolate. Returning {NAN}..." << std::endl;
			return {NAN};
		}

		const auto N = X.size();

		Container<Number> P(N);
		std::fill(P.begin(), P.end(), 0);

		Container<Number> l(N);

		for (SizeT k = 0; k < N; ++k) {
			l.resize(1);
			l[0] = 1;
			for (SizeT j = 0; j < N; ++j) {
				if (j != k) {
					// l = conv(l, [1/(X[k]-X[j]), -X[j]/(X[k]-X[j])]);
					l.resize(l.size() + 1);
					l[l.size()-1] = 0;
					for (SizeT i = l.size() - 1; i > 0; --i) {
						l[i] = (l[i] - X[j] * l[i-1]) / (X[k] - X[j]);
					}
					l[0] /= X[k] - X[j];
				}
			}
			for (SizeT i = 0; i < l.size(); ++i) {
				P[i] += Y[k] * l[i];
			}
		}

		return P;
	}
}