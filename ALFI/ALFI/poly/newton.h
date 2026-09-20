#pragma once

#include <iostream>
#include <cmath>

#include "../config.h"

namespace alfi::poly {
	template <typename Number = DefaultNumber, template <typename, typename...> class Container = DefaultContainer>
	Container<Number> newton(const Container<Number>& X, const Container<Number>& Y) {
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

		Container<Number> F = Y;

		for (SizeT i = 1; i < N; ++i) {
			for (SizeT j = N - 1; j >= i; --j) {
				F[j] = (F[j] - F[j-1]) / (X[j] - X[j-i]);
			}
		}

		Container<Number> P(N);
		std::fill(P.begin(), P.end() - 1, 0);
		P[P.size()-1] = F[0];

		Container<Number> f(N);

		for (SizeT i = 0; i < N - 1; ++i) {
			f.resize(1);
			f[0] = F[i+1];
			for (SizeT j = 0; j <= i; ++j) {
				// f = conv(f, [1, -X[j]]);
				f.resize(f.size() + 1);
				f[f.size()-1] = 0;
				for (SizeT k = f.size() - 1; k > 0; --k) {
					f[k] -= X[j] * f[k-1];
				}
			}
			// P = conv(P, [0, 1]); P = P + f;
			// P is buffered, so no need in conv
			for (SizeT j = 0; j < f.size(); ++j) {
				P[N-1-i-1+j] += f[j];
			}
		}

		return P;
	}
}