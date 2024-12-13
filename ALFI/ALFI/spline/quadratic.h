#pragma once

#include <iostream>
#include <cmath>

#include "../config.h"
#include "../util.h"
#include "util.h"

namespace alfi::spline {
	template <typename Number = DefaultNumber, typename Container = DefaultContainer<Number>>
	class QuadraticSpline {
	public:
		enum class Type {
			/**
			 * The second derivative at the first point is equal to zero.
			 */
			NATURAL_START,
			/**
			 * The second derivative at the last point is equal to zero.
			 */
			NATURAL_END,
			/**
			 * A parabola is constructed through the first three points. Then each subsequent segment is built sequentially.\n
			 * This is equivalent to the condition of continuity of the second derivative at the second point.\n
			 * In this way, the second point "ceases to be a knot point".
			 */
			NOT_A_KNOT_START,
			/**
			 * A parabola is constructed through the last three points. Then each previous segment is built sequentially.\n
			 * This is equivalent to the condition of continuity of the second derivative at the second-to-last point.\n
			 * In this way, the second-to-last point "ceases to be a knot point".
			 */
			NOT_A_KNOT_END,

			Default = NATURAL_START,
		};

		QuadraticSpline() = default;

		template <typename ContainerXType>
		QuadraticSpline(ContainerXType&& X, const Container& Y, Type type = Type::Default) {
			construct(std::forward<ContainerXType>(X), Y, type);
		}

		QuadraticSpline(const QuadraticSpline& other) = default;
		QuadraticSpline(QuadraticSpline&& other) noexcept = default;

		QuadraticSpline& operator=(const QuadraticSpline& other) = default;
		QuadraticSpline& operator=(QuadraticSpline&& other) noexcept = default;

		template <typename ContainerXType>
		void construct(ContainerXType&& X, const Container& Y, Type type = Type::Default) {
			if (X.size() != Y.size()) {
				std::cerr << "Error in function " << __FUNCTION__
						  << ": Vectors X (of size " << X.size()
						  << ") and Y (of size " << Y.size()
						  << ") are not the same size. Doing nothing..." << std::endl;
				return;
			}

			_type = type;
			const auto guard = alfi::util::SimpleScopeGuard([&]() {
				_X = std::forward<ContainerXType>(X);
			});

			const size_t n = X.size();

			if (n <= 1) {
				_coeffs = util::simple_spline<Number,Container>(X, Y, 2);
				return;
			}

			_coeffs.resize(3 * (n - 1));

			if (type == Type::NATURAL_START) {
				if (n <= 2) {
					_coeffs = util::simple_spline<Number,Container>(X, Y, 2);
					return;
				}
				_coeffs[0] = 0; // a1 = 0
				_coeffs[1] = (Y[1] - Y[0]) / (X[1] - X[0]); // b1
				_coeffs[2] = Y[0]; // c1
				for (size_t i = 1; i < n - 1; ++i) {
					_coeffs[3*i+1] = 2 * (Y[i] - Y[i-1]) / (X[i] - X[i-1]) - _coeffs[3*(i-1)+1];
					_coeffs[3*i] = ((Y[i+1] - Y[i]) / (X[i+1] - X[i]) - _coeffs[3*i+1]) / (X[i+1] - X[i]);
					// another way: _coeffs[3*i] = (_coeffs[3*(i+1)+1] - _coeffs[3*i+1]) / 2 / (X[i+1] - X[i]);
					_coeffs[3*i+2] = Y[i];
				}
			} else if (type == Type::NATURAL_END) {
				if (n <= 2) {
					_coeffs = util::simple_spline<Number,Container>(X, Y, 2);
					return;
				}
				_coeffs[3*(n-2)] = 0; // an-1 = 0
				_coeffs[3*(n-2)+1] = (Y[n-1] - Y[n-2]) / (X[n-1] - X[n-2]); // bn-1
				_coeffs[3*(n-2)+2] = Y[n-2]; // cn-1
				for (size_t iter = 0; iter <= n - 3; ++iter) {
					const size_t i = n - 3 - iter;
					_coeffs[3*i+1] = 2 * (Y[i+1] - Y[i]) / (X[i+1] - X[i]) - _coeffs[3*(i+1)+1];
					// another way: _coeffs[3*i] = ((Y[i+1] - Y[i]) / (X[i+1] - X[i]) - _coeffs[3*i+1]) / (X[i+1] - X[i]);
					_coeffs[3*i] = (_coeffs[3*(i+1)+1] - _coeffs[3*i+1]) / 2 / (X[i+1] - X[i]);
					_coeffs[3*i+2] = Y[i];
				}
			} else if (type == Type::NOT_A_KNOT_START) {
				if (n <= 3) {
					_coeffs = util::simple_spline<Number,Container>(X, Y, 2);
					return;
				}
				{ // first three points
					const Number h1 = X[1] - X[0], h2 = X[2] - X[1];
					const Number d1 = (Y[1] - Y[0]) / h1, d2 = (Y[2] - Y[1]) / h2;
					_coeffs[0] = (d2 - d1) / (h1 + h2); // a1
					_coeffs[1] = d1 - h1 * _coeffs[0];  // b1
					_coeffs[2] = Y[0];                  // c1
					_coeffs[3] = _coeffs[0];            // a2
					_coeffs[4] = d1 + h1 * _coeffs[0];  // b2
					_coeffs[5] = Y[1];                  // c2
				}
				for (size_t i = 2; i < n - 1; ++i) {
					const Number hi = X[i+1] - X[i], hi_1 = X[i] - X[i-1];
					const Number d = 2 * _coeffs[3*(i-1)] * hi_1 + _coeffs[3*(i-1)+1];
					_coeffs[3*i] = ((Y[i+1] - Y[i])/hi - d) / hi;
					_coeffs[3*i+1] = d;
					_coeffs[3*i+2] = Y[i];
				}
			} else if (type == Type::NOT_A_KNOT_END) {
				if (n <= 3) {
					_coeffs = util::simple_spline<Number,Container>(X, Y, 2);
					return;
				}
				{ // last three points
					const Number hn_2 = X[n-2] - X[n-3], hn_1 = X[n-1] - X[n-2];
					const Number dn_2 = (Y[n-2] - Y[n-3]) / hn_2, dn_1 = (Y[n-1] - Y[n-2]) / hn_1;
					_coeffs[3*(n-3)] = (dn_1 - dn_2) / (hn_2 + hn_1);    // an-2
					_coeffs[3*(n-3)+1] = dn_2 - hn_2 * _coeffs[3*(n-3)]; // bn-2
					_coeffs[3*(n-3)+2] = Y[n-3];                         // cn-2
					_coeffs[3*(n-2)] = _coeffs[3*(n-3)];                 // an-1
					_coeffs[3*(n-2)+1] = dn_2 + hn_2 * _coeffs[3*(n-3)]; // bn-1
					_coeffs[3*(n-2)+2] = Y[n-2];                         // cn-1
				}
				for (size_t iter = 0; iter <= n - 4; ++iter) {
					const size_t i = n - 4 - iter;
					const Number hi = X[i+1] - X[i];
					_coeffs[3*i] = (_coeffs[3*(i+1)+1] - (Y[i+1] - Y[i])/hi) / hi;
					_coeffs[3*i+1] = _coeffs[3*(i+1)+1] - 2*_coeffs[3*i]*hi;
					// another way: _coeffs[3*i+1] = (Y[i+1] - Y[i])/hi - _coeffs[3*i]*hi;
					_coeffs[3*i+2] = Y[i];
				}
			} else {
				std::cerr << "Error in function " << __FUNCTION__
						  << ": unknown type (" << static_cast<int>(type)
						  << "). Setting _coeffs to an empty array..." << std::endl;
				_coeffs = {};
			}
		}

		Number eval(Number x) const {
			return eval(x, std::distance(_X.begin(), alfi::util::first_leq_or_begin(_X.begin(), _X.end(), x)));
		}
		Number eval(Number x, size_t segment_index) const {
			if (_X.empty())
				return NAN;
			segment_index = std::max(std::min(segment_index, _X.size() - 2), static_cast<size_t>(0));
			x = x - _X[segment_index];
			return (_coeffs[3*segment_index] * x + _coeffs[3*segment_index+1]) * x + _coeffs[3*segment_index+2];
		}

		Container eval(const Container& evalX, bool sorted = true) const {
			Container result(evalX.size());
			if (sorted) {
				for (size_t i = 0, i_x = 0; i < evalX.size(); ++i) {
					const Number evalx = evalX[i];
					while (i_x < _X.size() && evalx >= _X[i_x + 1])
						++i_x;
					result[i] = eval(evalx, i_x);
				}
			} else {
				for (size_t i = 0; i < evalX.size(); ++i) {
					result[i] = eval(evalX[i]);
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

		Type type() const {
			return _type;
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
			return {3*index, 3*(index+1)};
		}

	private:
		Type _type = Type::Default;
		Container _X = {};
		Container _coeffs = {};
	};
}