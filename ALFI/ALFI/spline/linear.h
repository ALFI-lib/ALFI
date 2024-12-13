#pragma once

#include <iostream>
#include <cmath>

#include "../config.h"
#include "../util.h"

namespace alfi::spline {
	template <typename Number = DefaultNumber, typename Container = DefaultContainer<Number>>
	class LinearSpline {
	public:
		LinearSpline() = default;

		template <typename ContainerXType>
		LinearSpline(ContainerXType&& X, const Container& Y) {
			construct(std::forward<ContainerXType>(X), Y);
		}

		LinearSpline(const LinearSpline& other) = default;
		LinearSpline(LinearSpline&& other) noexcept = default;

		LinearSpline& operator=(const LinearSpline& other) = default;
		LinearSpline& operator=(LinearSpline&& other) noexcept = default;

		template <typename ContainerXType>
		void construct(ContainerXType&& X, const Container& Y) {
			if (X.size() != Y.size()) {
				std::cerr << "Error in function " << __FUNCTION__
						  << ": Vectors X (of size " << X.size()
						  << ") and Y (of size " << Y.size()
						  << ") are not the same size. Doing nothing..." << std::endl;
				return;
			}

			const auto guard = util::SimpleScopeGuard([&]() {
				_X = std::forward<ContainerXType>(X);
			});

			const size_t n = X.size();

			if (n <= 2) {
				if (n == 0) {
					_coeffs = {};
				} else if (n == 1) {
					_coeffs = {Y[0]};
				} else {
					_coeffs = {(Y[1] - Y[0]) / (X[1] - X[0]), Y[0]};
				}
				return;
			}

			_coeffs.resize(2 * (n - 1));

			for (size_t i = 0; i < n - 1; ++i) {
				_coeffs[2*i] = (Y[i+1] - Y[i]) / (X[i+1] - X[i]);
				_coeffs[2*i+1] = Y[i];
			}
		}

		Number eval(Number x) const {
			return eval(x, std::distance(_X.begin(), util::first_leq_or_begin(_X.begin(), _X.end(), x)));
		}
		Number eval(Number x, size_t segment_index) const {
			if (_X.empty())
				return NAN;
			segment_index = std::max(std::min(segment_index, _X.size() - 2), static_cast<size_t>(0));
			x = x - _X[segment_index];
			return _coeffs[2*segment_index] * x + _coeffs[2*segment_index+1];
		}

		Container eval(const Container& xx, bool sorted = true) const {
			Container result(xx.size());
			if (sorted) {
				for (size_t i = 0, i_x = 0; i < xx.size(); ++i) {
					const Number x = xx[i];
					while (i_x < _X.size() && x >= _X[i_x+1])
						++i_x;
					result[i] = eval(x, i_x);
				}
			} else {
				for (size_t i = 0; i < xx.size(); ++i) {
					result[i] = eval(xx[i]);
				}
			}
			return result;
		}

		Number operator()(Number x) const {
			return eval(x);
		}
		Container operator()(const Container& xx) const {
			return eval(xx);
		}

		const Container& X() const & {
			return _X;
		}
		Container&& X() && {
			return std::move(_X);
		}

		const Container& coeffs() const & {
			return _coeffs;
		}
		Container&& coeffs() && {
			return std::move(_coeffs);
		}

		static std::pair<size_t, size_t> segment(size_t index) {
			return {2*index, 2*(index+1)};
		}

	private:
		Container _X = {};
		Container _coeffs = {};
	};
}