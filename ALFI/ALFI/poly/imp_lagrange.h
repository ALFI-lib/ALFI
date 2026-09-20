#pragma once

#include <iostream>
#include <cmath>

#include "../config.h"

namespace alfi::poly {
	template <typename Number = DefaultNumber, template <typename, typename...> class Container = DefaultContainer>
	Container<Number> imp_lagrange(const Container<Number>& X, const Container<Number>& Y) {
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

		Container<Number> w_rev(N);
		std::fill(w_rev.begin(), w_rev.end(), 1);
		for (SizeT i = 0; i < N; ++i) {
			for (SizeT j = 0; j < N; ++j) {
				if (i != j) {
					w_rev[i] *= (X[i] - X[j]);
				}
			}
		}

		Container<Number> P(N);
		std::fill(P.begin(), P.end(), 0);

		Container<Number> l(N + 1);

		for (SizeT k = 0; k < N; ++k) {
			l.resize(1);
			l[0] = 1;
			for (SizeT j = 0; j < N; ++j) {
				// l = conv(l, [1, -X[j]]);
				l.resize(l.size() + 1);
				l[l.size()-1] = 0;
				for (SizeT i = l.size() - 1; i > 0; --i) {
					l[i] -= X[j] * l[i-1];
				}
			}
		}

		for (SizeT k = 0; k < N; ++k) {
			Container<Number> l_cur(N);
			l_cur[0] = l[0];
			for (SizeT i = 1; i < N; ++i) {
				l_cur[i] = l[i] + l_cur[i-1] * X[k];
			}
			for (SizeT i = 0; i < N; ++i) {
				P[i] += Y[k] * l_cur[i] / w_rev[k];
			}
		}

		return P;
	}
}