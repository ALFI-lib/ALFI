#pragma once

#include <iostream>
#include <cmath>

#include "config.h"
#include "util/numeric.h"

namespace alfi::poly {
	template <typename Number = DefaultNumber>
	std::enable_if_t<!traits::has_size<Number>::value, Number>
	val(const auto& coeffs, Number x) {
		Number result = 0;
		for (const Number& c : coeffs) {
			result = x * result + c;
		}
		return result;
	}

	template <typename Number = DefaultNumber, template <typename> class Container = DefaultContainer, typename ContainerX = Container<Number>>
	std::enable_if_t<traits::has_size<ContainerX>::value, Container<Number>>
	val(const auto& coeffs, const ContainerX& X) {
		Container<Number> result(X.size(), 0);
#if defined(_OPENMP) && !defined(ALFI_DISABLE_OPENMP)
#pragma omp parallel for
#endif
		for (SizeT i = 0; i < X.size(); ++i) {
			for (const Number& c : coeffs) {
				result[i] = X[i] * result[i] + c;
			}
		}
		return result;
	}

	template <typename Number = DefaultNumber, template <typename> class Container = DefaultContainer>
	Container<Number> lagrange(const auto& X, const auto& Y) {
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

	template <typename Number = DefaultNumber, template <typename> class Container = DefaultContainer>
	Container<Number> lagrange_vals(const auto& X, const auto& Y, const auto& xx) {
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

	template <typename Number = DefaultNumber, template <typename> class Container = DefaultContainer>
	Container<Number> imp_lagrange(const auto& X, const auto& Y) {
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

	template <typename Number = DefaultNumber, template <typename> class Container = DefaultContainer>
	Container<Number> imp_lagrange_vals(const auto& X, const auto& Y, const auto& xx, Number epsilon = std::numeric_limits<Number>::epsilon()) {
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

		Container<Number> w_rev(N);
		std::fill(w_rev.begin(), w_rev.end(), 1);
		for (SizeT i = 0; i < N; ++i) {
			for (SizeT j = 0; j < N; ++j) {
				if (i != j) {
					w_rev[i] *= (X[i] - X[j]);
				}
			}
		}

		Container<Number> yy(nn);

#if defined(_OPENMP) && !defined(ALFI_DISABLE_OPENMP)
#pragma omp parallel for
#endif
		for (SizeT k = 0; k < nn; ++k) {
			Number l = 1;
			for (SizeT i = 0; i < N; ++i) {
				l *= (xx[k] - X[i]);
			}
			yy[k] = 0;
			for (SizeT i = 0; i < N; ++i) {
				if (util::numeric::are_equal(xx[k], X[i], epsilon)) {
					yy[k] = Y[i];
					break;
				} else {
					yy[k] += Y[i] * l / (w_rev[i] * (xx[k] - X[i]));
				}
			}
		}

		return yy;
	}

	template <typename Number = DefaultNumber, template <typename> class Container = DefaultContainer>
	Container<Number> newton(const auto& X, const auto& Y) {
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

	template <typename Number = DefaultNumber, template <typename> class Container = DefaultContainer>
	Container<Number> newton_vals(const auto& X, const auto& Y, const auto& xx) {
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

		Container<Number> F = Y;

		for (SizeT i = 1; i < N; ++i) {
			for (SizeT j = N - 1; j >= i; --j) {
				F[j] = (F[j] - F[j-1]) / (X[j] - X[j-i]);
			}
		}

		Container<Number> yy(nn);

#if defined(_OPENMP) && !defined(ALFI_DISABLE_OPENMP)
#pragma omp parallel for
#endif
		for (SizeT k = 0; k < nn; ++k) {
			yy[k] = F[N-1];
			for (SizeT iter = 0; iter < N - 1; ++iter) {
				const auto i = N - 2 - iter;
				yy[k] = yy[k] * (xx[k] - X[i]) + F[i];
			}
		}

		return yy;
	}
}