#pragma once

#include <iostream>

#include "config.h"

namespace alfi::poly {

	template <typename Number = DefaultNumber, typename Container = DefaultContainer<Number>>
	Number val(const Container& coeffs, Number x) {
		Number result = 0;
		for (const Number& c : coeffs) {
			result = x * result + c;
		}
		return result;
	}

	template <typename Number = DefaultNumber, typename Container = DefaultContainer<Number>>
	Container val(const Container& coeffs, const Container& X) {
		Container result(X.size(), 0);
		for (const Number& c : coeffs) {
			for (size_t i = 0; i < X.size(); ++i) {
				result[i] = X[i] * result[i] + c;
			}
		}
		return result;
	}

	template <typename Number = DefaultNumber, typename Container = DefaultContainer<Number>>
	Number val(const Container& coeffs, Number x, Number x0) {
		Number result = 0;
		for (const Number& c : coeffs) {
			result = (x - x0) * result + c;
		}
		return result;
	}

	template <typename Number = DefaultNumber, typename Container = DefaultContainer<Number>>
	Container val(const Container& coeffs, const Container& X, Number x0) {
		Container result(X.size(), 0);
		for (const Number& c : coeffs) {
			for (size_t i = 0; i < X.size(); ++i) {
				result[i] = (X[i] - x0) * result[i] + c;
			}
		}
		return result;
	}

	template <typename Number = DefaultNumber, typename Container = DefaultContainer<Number>>
	Container lagrange(const Container& X, const Container& Y, Number x0 = 0) {
		if (X.size() != Y.size()) {
			std::cerr << "Error in function " << __FUNCTION__
					  << ": Vectors X (of size " << X.size()
					  << ") and Y (of size " << Y.size()
					  << ") are not the same size. Returning an empty array..." << std::endl;
			return {};
		}

		if (X.empty()) {
			std::cerr << "Error in function " << __FUNCTION__
					  << ": Vectors X and Y are empty. Cannot interpolate. Returning an empty array..." << std::endl;
			return {};	
		}

		const size_t n = X.size();

		Container L(n, 0);

		Container l(n);

		for (size_t k = 0; k < n; ++k) {
			l.resize(1);
			l[0] = 1;
			for (size_t j = 0; j < n; ++j) {
				if (j != k) {
					// A = conv(A, [1/(X[k]-X[j]), -(X[j] - x0)/(X[k]-X[j])]);
					l.resize(l.size() + 1);
					l[l.size()-1] = 0;
					for (size_t i = l.size() - 1; i > 0; --i)
						l[i] = (l[i] - (X[j] - x0) * l[i-1]) / (X[k] - X[j]);
					l[0] /= X[k] - X[j];
				}
			}
			for (size_t i = 0; i < l.size(); ++i) {
				L[i] += Y[k] * l[i];
			}
		}

		return L;
	}

	template <typename Number = DefaultNumber, typename Container = DefaultContainer<Number>>
	Container newton(const Container& X, const Container& Y, Number x0 = 0) {
		if (X.size() != Y.size()) {
			std::cerr << "Error in function " << __FUNCTION__
					  << ": Vectors X (of size " << X.size()
					  << ") and Y (of size " << Y.size()
					  << ") are not the same size. Returning an empty array..." << std::endl;
			return {};
		}

		if (X.empty()) {
			std::cerr << "Error in function " << __FUNCTION__
					  << ": Vectors X and Y are empty. Cannot interpolate. Returning an empty array..." << std::endl;
			return {};	
		}

		const size_t n = X.size();

		Container F = Y;

		for (size_t i = 1; i < n; ++i)
			for (size_t j = n - 1; j >= i; --j)
				F[j] = (F[j] - F[j-1]) / (X[j] - X[j-i]);

		Container N(n, 0);
		N[N.size()-1] = F[0];

		Container f(n);
		for (size_t i = 0; i < n - 1; ++i) {
			f.resize(1);
			f[0] = F[i+1];
			for (size_t j = 0; j <= i; ++j) {
				// f = conv(f, [1, -(X[j] - x0)]);
				f.resize(f.size() + 1);
				f[f.size()-1] = 0;
				for (size_t k = f.size() - 1; k > 0; --k)
					f[k] -= (X[j] - x0) * f[k-1];
			}
			// N = conv(N, [0, 1]); N = N + f;
			// N is buffered, so no need in conv
			for (size_t j = 0; j < f.size(); ++j)
				N[n-1-i-1+j] += f[j];
		}

		return N;
	}
}