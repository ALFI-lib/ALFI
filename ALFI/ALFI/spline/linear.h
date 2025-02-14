#pragma once

#include <algorithm>
#include <iostream>
#include <cmath>

#include "../config.h"
#include "../util/misc.h"

namespace alfi::spline {
	template <typename Number = DefaultNumber, template <typename> class Container = DefaultContainer>
	class LinearSpline {
	public:
		static Container<Number> compute_coeffs(const auto& X, const auto& Y) {
			if (X.size() != Y.size()) {
				std::cerr << "Error in function " << __FUNCTION__
						  << ": Vectors X (of size " << X.size()
						  << ") and Y (of size " << Y.size()
						  << ") are not the same size. Returning an empty array..." << std::endl;
				return {};
			}

			const auto n = X.size();

			if (n <= 2) {
				if (n == 0) {
					return {};
				} else if (n == 1) {
					return {Y[0]};
				} else {
					return {(Y[1] - Y[0]) / (X[1] - X[0]), Y[0]};
				}
			}

			Container<Number> coeffs;

			coeffs.resize(2 * (n - 1));

			for (SizeT i = 0; i < n - 1; ++i) {
				coeffs[2*i] = (Y[i+1] - Y[i]) / (X[i+1] - X[i]);
				coeffs[2*i+1] = Y[i];
			}

			return coeffs;
		}

		LinearSpline() = default;

		template <typename ContainerXType>
		LinearSpline(ContainerXType&& X, const auto& Y) {
			construct(std::forward<ContainerXType>(X), Y);
		}

		LinearSpline(const LinearSpline& other) = default;
		LinearSpline(LinearSpline&& other) noexcept = default;

		LinearSpline& operator=(const LinearSpline& other) = default;
		LinearSpline& operator=(LinearSpline&& other) noexcept = default;

		template <typename ContainerXType>
		void construct(ContainerXType&& X, const auto& Y) {
			if (X.size() != Y.size()) {
				std::cerr << "Error in function " << __FUNCTION__
						  << ": Vectors X (of size " << X.size()
						  << ") and Y (of size " << Y.size()
						  << ") are not the same size. Doing nothing..." << std::endl;
				return;
			}
			auto coeffs = compute_coeffs(X, Y);
			if (coeffs.empty() && !X.empty()) {
				std::cerr << "Error in function " << __FUNCTION__
						  << ": failed to construct coefficients. Not changing the object..." << std::endl;
				return;
			}
			_X = std::forward<ContainerXType>(X);
			_coeffs = std::move(coeffs);
		}

		Number eval(Number x) const {
			return eval(x, std::distance(_X.begin(), util::misc::first_leq_or_begin(_X.begin(), _X.end(), x)));
		}
		Number eval(Number x, SizeT segment) const {
			if (_coeffs.empty()) {
				return NAN;
			} else if (_coeffs.size() == 1) {
				return _coeffs[0];
			}
			segment = std::clamp(segment, static_cast<SizeT>(0), static_cast<SizeT>(_X.size() - 2));
			x = x - _X[segment];
			return _coeffs[2*segment] * x + _coeffs[2*segment+1];
		}

		Container<Number> eval(const Container<Number>& xx, bool sorted = true) const {
			Container<Number> result(xx.size());
			if (sorted) {
				for (SizeT i = 0, i_x = 0; i < xx.size(); ++i) {
					const Number x = xx[i];
					while (i_x + 1 < _X.size() && x >= _X[i_x+1])
						++i_x;
					result[i] = eval(x, i_x);
				}
			} else {
				for (SizeT i = 0; i < xx.size(); ++i) {
					result[i] = eval(xx[i]);
				}
			}
			return result;
		}

		Number operator()(Number x) const {
			return eval(x);
		}
		Container<Number> operator()(const Container<Number>& xx) const {
			return eval(xx);
		}

		const Container<Number>& X() const & {
			return _X;
		}
		Container<Number>&& X() && {
			return std::move(_X);
		}

		const Container<Number>& coeffs() const & {
			return _coeffs;
		}
		Container<Number>&& coeffs() && {
			return std::move(_coeffs);
		}

		static std::pair<SizeT, SizeT> segment(SizeT index) {
			return {2*index, 2*(index+1)};
		}

	private:
		Container<Number> _X = {};
		Container<Number> _coeffs = {};
	};
}