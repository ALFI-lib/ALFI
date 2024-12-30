#pragma once

#include <iostream>
#include <cmath>

#include "../config.h"
#include "../util/arrays.h"
#include "../util/misc.h"
#include "../util/spline.h"

namespace alfi::spline {
	template <typename Number = DefaultNumber, template <typename> class Container = DefaultContainer>
	class CubicSpline {
	public:
		enum class Type {
			/**
			 * Second derivatives at the end points are equal to zero.\n
			 * Has minimum curvature of all interpolating functions and is most stable.\n
			 */
			NATURAL,
			/**
			 * The third derivative is continuous at the second and second-to-last points.
			 * In this way, the second and second-to-last points "cease to be knot points".
			 */
			NOT_A_KNOT,
			/**
			 * The first and second derivatives at the end points are equal.
			 */
			PERIODIC,
			/**
			 * A cubic curve is constructed through the first four points. Then each subsequent segment is built sequentially.\n
			 * This is equivalent to the condition of continuity of the third derivative at the second and third points.\n
			 * In this way, the second and third points "cease to be knot points".
			 */
			NOT_A_KNOT_START,
			/**
			 * A cubic curve is constructed through the last four points. Then each previous segment is built sequentially.\n
			 * This is equivalent to the condition of continuity of the third derivative at the third-to-last and second-to-last points.\n
			 * In this way, the third-to-last and second-to-last points "cease to be knot points".
			 */
			NOT_A_KNOT_END,
			/**
			 * The arithmetic mean_array between NOT_A_KNOT_START and NOT_A_KNOT_END.
			 */
			SEMI_NOT_A_KNOT,

			Default = NATURAL,
		};

		static Container<Number> compute_coeffs(const auto& X, const auto& Y, Type type = Type::Default) {
			if (X.size() != Y.size()) {
				std::cerr << "Error in function " << __FUNCTION__
						  << ": Vectors X (of size " << X.size()
						  << ") and Y (of size " << Y.size()
						  << ") are not the same size. Returning an empty array..." << std::endl;
				return {};
			}

			const SizeT n = X.size();

			if (n <= 1) {
				return util::spline::simple_spline<Number,Container>(X, Y, 3);
			}

			Container<Number> coeffs(4 * (n - 1));

			if (type == Type::NATURAL) {
				if (n <= 2) {
					return util::spline::simple_spline<Number,Container>(X, Y, 3);
				}

				Container<Number> dX(n - 1), dY(n - 1);
				for (SizeT i = 0; i < n - 1; ++i) {
					dX[i] = X[i+1] - X[i];
					dY[i] = Y[i+1] - Y[i];
				}

				Container<Number> B(n);
				B[0] = B[n-1] = 0;
				// Tridiagonal matrix algorithm
				Container<Number> diagonal(n - 2), right(n - 2);
				for (SizeT i = 0; i < n - 2; ++i) {
					diagonal[i] = 2 * (dX[i] + dX[i+1]);
					right[i] = 3 * (dY[i+1] / dX[i+1] - dY[i] / dX[i]);
				}
				// Forward pass
				for (SizeT i = 1; i < n - 2; ++i) {
					const auto m = dX[i] / diagonal[i-1];
					diagonal[i] -= m * dX[i];
					right[i] -= m * right[i-1];
				}
				// Backward pass
				B[n-2] = right[n-3] / diagonal[n-3];
				for (SizeT i = n - 3; i >= 1; --i) {
					B[i] = (right[i-1] - dX[i] * B[i+1]) / diagonal[i-1];
				}

				Container<Number> A(n - 1);
				for (SizeT i = 0; i < A.size(); ++i) {
					A[i] = (B[i+1] - B[i]) / (3*dX[i]);
				}

				Container<Number> C(n - 1);
				for (SizeT i = 0; i < C.size(); ++i) {
					C[i] = dY[i] / dX[i] - dX[i] * ((2*B[i] + B[i+1]) / 3);
				}

				// const auto D = Y; - no need to create extra array

				for (SizeT i = 0, index = 0; i < n - 1; ++i) {
					coeffs[index++] = A[i];
					coeffs[index++] = B[i];
					coeffs[index++] = C[i];
					coeffs[index++] = Y[i]; // = D[i];
				}
			} else if (type == Type::NOT_A_KNOT) {
				if (n <= 4) {
					return util::spline::simple_spline<Number,Container>(X, Y, 3);
				}

				Container<Number> dX(n - 1), dY(n - 1);
				for (SizeT i = 0; i < n - 1; ++i) {
					dX[i] = X[i+1] - X[i];
					dY[i] = Y[i+1] - Y[i];
				}

				const auto big_denominator = dX[n-2] * (dX[n-3] + dX[n-2]) * (2*dX[n-3] + dX[n-2]);

				Container<Number> B(n - 1);
				// Tridiagonal matrix algorithm
				Container<Number> diagonal(n - 3), right(n - 3);
				diagonal[0] = ((dX[0] + dX[1]) * (dX[0] + 2 * dX[1])) / dX[1];
				for (SizeT i = 1; i < n - 4; ++i) {
					diagonal[i] = 2 * (dX[i] + dX[i+1]);
				}
				diagonal[n-4] = 2 * dX[n-4] + 3 * dX[n-3] * (1 - dX[n-3] / (2 * dX[n-3] + dX[n-2]));
				for (SizeT i = 0; i < n - 4; ++i) {
					right[i] = 3 * (dY[i+1] / dX[i+1] - dY[i] / dX[i]);
				}
				right[n-4] = 3 * (dY[n-3]/dX[n-3] - dY[n-4]/dX[n-4]) - (3*dX[n-3]*(dY[n-2]*dX[n-3] - dY[n-3]*dX[n-2])) / big_denominator;
				// Forward pass
				diagonal[1] -= dX[1] / diagonal[0] * (dX[1] - (dX[0] * dX[0]) / dX[1]);
				right[1] -= dX[1] / diagonal[0] * right[0];
				for (SizeT i = 2; i < n - 3; ++i) {
					const auto m = dX[i] / diagonal[i-1];
					diagonal[i] -= m * dX[i];
					right[i] -= m * right[i-1];
				}
				// Backward pass
				B[n-3] = right[n-4] / diagonal[n-4];
				for (SizeT i = n - 4; i >= 2; --i) {
					B[i] = (right[i-1] - dX[i] * B[i+1]) / diagonal[i-1];
				}
				B[1] = (right[0] - (dX[1] - (dX[0] * dX[0]) / dX[1]) * B[2]) / diagonal[0];

				B[n-2] = (dX[n-2]*(B[n-3]*(dX[n-2]*dX[n-2] - dX[n-3]*dX[n-3]) - 3*dY[n-3]) + 3*dX[n-3]*dY[n-2]) / big_denominator;

				Container<Number> A(n - 1);
				for (SizeT i = 1; i < A.size() - 1; ++i) {
					A[i] = (B[i+1] - B[i]) / (3*dX[i]);
				}
				A[0] = A[1];
				A[n-2] = A[n-3];

				B[0] = B[1] - 3 * A[0] * dX[0];

				Container<Number> C(n - 1);
				for (SizeT i = 0; i < C.size() - 1; ++i) {
					C[i] = dY[i] / dX[i] - dX[i] * ((2*B[i] + B[i+1]) / 3);
				}
				C[n-2] = dY[n-2] / dX[n-2] - A[n-2] * dX[n-2] * dX[n-2] - B[n-2] * dX[n-2];

				// const auto D = Y; - no need to create extra array

				for (SizeT i = 0, index = 0; i < n - 1; ++i) {
					coeffs[index++] = A[i];
					coeffs[index++] = B[i];
					coeffs[index++] = C[i];
					coeffs[index++] = Y[i]; // = D[i];
				}
			} else if (type == Type::PERIODIC) {
				if (n <= 2) {
					return util::spline::simple_spline<Number,Container>(X, Y, 3);
				}

				Container<Number> dX(n - 1), dY(n - 1);
				for (SizeT i = 0; i < n - 1; ++i) {
					dX[i] = X[i+1] - X[i];
					dY[i] = Y[i+1] - Y[i];
				}

				Container<Number> B(n);
				if (n == 3) {
					B[0] = 3 * (dY[0] / dX[0] - dY[1] / dX[1]) / (dX[0] + dX[1]);
					B[1] = -B[0];
					B[2] = B[0];
				} else {
					Container<Number> diag(n - 1), right(n - 1);
					for (SizeT i = 0; i < n - 2; ++i) {
						diag[i] = 2 * (dX[i] + dX[i+1]);
						right[i] = 3 * (dY[i+1] / dX[i+1] - dY[i] / dX[i]);
					}
					diag[n-2] = 2 * (dX[n-2] + dX[0]);
					right[n-2] = 3 * (dY[0]/dX[0] - dY[n-2]/dX[n-2]);

					Container<Number> last_row(n - 2), last_col(n - 2);
					last_row[0] = last_col[0] = dX[0];
					for (SizeT i = 1; i < n - 3; ++i) {
						last_row[i] = last_col[i] = 0;
					}
					last_row[n-3] = last_col[n-3] = dX[n-2];

					for (SizeT i = 0; i < n - 3; ++i) {
						const auto m1 = dX[i+1] / diag[i];
						diag[i+1] -= m1 * dX[i+1];
						last_col[i+1] -= m1 * last_col[i];
						right[i+1] -= m1 * right[i];
						const auto m2 = last_row[i] / diag[i];
						last_row[i+1] -= m2 * dX[i+1];
						diag[n-2] -= m2 * last_col[i];
						right[n-2] -= m2 * right[i];
					}
					diag[n-2] -= last_row[n-3] / diag[n-3] * last_col[n-3];
					right[n-2] -= last_row[n-3] / diag[n-3] * right[n-3];

					B[n-1] = right[n-2] / diag[n-2];
					B[n-2] = (right[n-3] - last_col[n-3] * B[n-1]) / diag[n-3];
					for (SizeT i = n - 3; i >= 1; --i) {
						B[i] = (right[i-1] - dX[i]*B[i+1] - last_col[i-1]*B[n-1]) / diag[i-1];
					}
					B[0] = B[n-1];
				}

				Container<Number> A(n - 1);
				for (SizeT i = 0; i < A.size(); ++i) {
					A[i] = (B[i+1] - B[i]) / (3*dX[i]);
				}

				Container<Number> C(n - 1);
				for (SizeT i = 0; i < C.size(); ++i) {
					C[i] = dY[i] / dX[i] - dX[i] * ((2*B[i] + B[i+1]) / 3);
				}

				for (SizeT i = 0, index = 0; i < n - 1; ++i) {
					coeffs[index++] = A[i];
					coeffs[index++] = B[i];
					coeffs[index++] = C[i];
					coeffs[index++] = Y[i]; // = D[i];
				}
			} else if (type == Type::NOT_A_KNOT_START) {
				if (n <= 4) {
					return util::spline::simple_spline<Number,Container>(X, Y, 3);
				}

				Container<Number> dX(n - 1), dY(n - 1);
				for (SizeT i = 0; i < n - 1; ++i) {
					dX[i] = X[i+1] - X[i];
					dY[i] = Y[i+1] - Y[i];
				}

				{ // first four points
					const auto h1 = dX[0], h2 = dX[1];
					const auto h12 = X[2]-X[0], h13 = X[3]-X[0];
					const auto d1 = Y[1]-Y[0], d12 = Y[2]-Y[0], d13 = Y[3]-Y[0];

					const auto abc = util::linalg::linsolve(
						{{std::pow(h1, 3), std::pow(h1, 2), h1},
							{std::pow(h12, 3), std::pow(h12, 2), h12},
							{std::pow(h13, 3), std::pow(h13, 2), h13}},
						{d1, d12, d13});
					if (abc.empty()) {
						std::cerr << "Error in function " << __FUNCTION__ << ": Could not compute. Returning an empty array..." << std::endl;
						return {};
					}

					const auto a = abc[0], b = abc[1], c = abc[2];

					coeffs[0] = a; // a1
					coeffs[1] = b; // b1
					coeffs[2] = c; // c1
					coeffs[3] = Y[0]; // d1
					coeffs[4] = a; // a2
					coeffs[5] = 3 * a * h1 + b; // b2
					coeffs[6] = 3 * a * std::pow(h1, 2) + 2 * b * h1 + c; // c2
					coeffs[7] = Y[1]; // d2
					coeffs[8] = a; // a3
					coeffs[9] = 3 * a * h2 + coeffs[5]; // b3
					coeffs[10] = 3 * a * std::pow(h2, 2) + 2 * coeffs[5] * h2 + coeffs[6]; // c3
					coeffs[11] = Y[2]; // d3
				}

				for (SizeT i = 3; i < n - 1; ++i) {
					coeffs[4*i+3] = Y[i];
					coeffs[4*i+2] = 3 * coeffs[4*(i-1)] * pow(dX[i-1], 2) + 2 * coeffs[4*(i-1)+1] * dX[i-1] + coeffs[4*(i-1)+2];
					coeffs[4*i+1] = 3 * coeffs[4*(i-1)] * dX[i-1] + coeffs[4*(i-1)+1];
					coeffs[4*i] = (dY[i]/dX[i] - coeffs[4*i+1] * dX[i] - coeffs[4*i+2]) / pow(dX[i], 2);
				}
			} else if (type == Type::NOT_A_KNOT_END) {
				if (n <= 4) {
					return util::spline::simple_spline<Number,Container>(X, Y, 3);
				}

				Container<Number> dX(n - 1), dY(n - 1);
				for (SizeT i = 0; i < n - 1; ++i) {
					dX[i] = X[i+1] - X[i];
					dY[i] = Y[i+1] - Y[i];
				}

				{ // last four points
					const auto h1 = dX[n-4], h2 = dX[n-3];
					const auto h12 = X[n-2]-X[n-4], h13 = X[n-1]-X[n-4];
					const auto d1 = dY[n-4], d12 = Y[n-2] - Y[n-4], d13 = Y[n-1] - Y[n-4];

					const auto abc = util::linalg::linsolve(
						{{std::pow(h1, 3), std::pow(h1, 2), h1},
							{std::pow(h12, 3), std::pow(h12, 2), h12},
							{std::pow(h13, 3), std::pow(h13, 2), h13}},
						{d1, d12, d13});
					if (abc.empty()) {
						std::cerr << "Error in function " << __FUNCTION__ << ": Could not compute. Returning an empty array..." << std::endl;
						return {};
					}

					const auto a = abc[0], b = abc[1], c = abc[2];

					coeffs[4*(n-4)] = a; // an-3
					coeffs[4*(n-4)+1] = b; // bn-3
					coeffs[4*(n-4)+2] = c; // cn-3
					coeffs[4*(n-4)+3] = Y[n-4]; // dn-3
					coeffs[4*(n-3)] = a; // an-2
					coeffs[4*(n-3)+1] = 3 * a * h1 + b; // bn-2
					coeffs[4*(n-3)+2] = 3 * a * std::pow(h1, 2) + 2 * b * h1 + c; // cn-2
					coeffs[4*(n-3)+3] = Y[n-3]; // dn-2
					coeffs[4*(n-2)] = a; // an-1
					coeffs[4*(n-2)+1] = 3 * a * h2 + coeffs[4*(n-3)+1]; // bn-1
					coeffs[4*(n-2)+2] = 3 * a * std::pow(h2, 2) + 2 * coeffs[4*(n-3)+1] * h2 + coeffs[4*(n-3)+2]; // cn-1
					coeffs[4*(n-2)+3] = Y[n-2]; // dn-1
				}

				for (SizeT iter = 0; iter <= n - 5; ++iter) {
					const auto i = n - 5 - iter;
					coeffs[4*i] = (coeffs[4*(i+1)+1] - coeffs[4*(i+1)+2]/dX[i] + dY[i]/pow(dX[i], 2)) / dX[i];
					coeffs[4*i+1] = coeffs[4*(i+1)+1] - 3 * coeffs[4*i] * dX[i];
					coeffs[4*i+2] = dY[i]/dX[i] - coeffs[4*i]*pow(dX[i], 2) - coeffs[4*i+1]*dX[i];
					coeffs[4*i+3] = Y[i];
				}
			} else if (type == Type::SEMI_NOT_A_KNOT) {
				return util::arrays::mean(compute_coeffs(X, Y, Type::NOT_A_KNOT_START), compute_coeffs(X, Y, Type::NOT_A_KNOT_END));
			} else {
				std::cerr << "Error in function " << __FUNCTION__
						  << ": unknown type (" << static_cast<int>(type)
						  << "). Returning an empty array..." << std::endl;
				return {};
			}

			return coeffs;
		}

		CubicSpline() = default;

		template <typename ContainerXType>
		CubicSpline(ContainerXType&& X, const auto& Y, Type type = Type::Default) {
			construct(std::forward<ContainerXType>(X), Y, type);
		}

		CubicSpline(const CubicSpline& other) = default;
		CubicSpline(CubicSpline&& other) noexcept = default;

		CubicSpline& operator=(const CubicSpline& other) = default;
		CubicSpline& operator=(CubicSpline&& other) noexcept = default;

		template <typename ContainerXType>
		void construct(ContainerXType&& X, const auto& Y, Type type = Type::Default) {
			if (X.size() != Y.size()) {
				std::cerr << "Error in function " << __FUNCTION__
						  << ": Vectors X (of size " << X.size()
						  << ") and Y (of size " << Y.size()
						  << ") are not the same size. Doing nothing..." << std::endl;
				return;
			}
			auto coeffs = compute_coeffs(X, Y, type);
			if (coeffs.empty() && !X.empty()) {
				std::cerr << "Error in function " << __FUNCTION__
						  << ": failed to construct coefficients. Not changing the object..." << std::endl;
				return;
			}
			_type = type;
			_X = std::forward<ContainerXType>(X);
			_coeffs = std::move(coeffs);
		}

		Number eval(Number x) const {
			return eval(x, std::distance(_X.begin(), util::misc::first_leq_or_begin(_X.begin(), _X.end(), x)));
		}
		Number eval(Number x, SizeT segment_index) const {
			if (_coeffs.empty()) {
				return NAN;
			} else if (_coeffs.size() == 1) {
				return _coeffs[0];
			}
			segment_index = std::max(std::min(segment_index, static_cast<SizeT>(_X.size() - 2)), static_cast<SizeT>(0));
			x = x - _X[segment_index];
			return ((_coeffs[4*segment_index] * x + _coeffs[4*segment_index+1]) * x + _coeffs[4*segment_index+2]) * x + _coeffs[4*segment_index+3];
		}

		Container<Number> eval(const Container<Number>& xx, bool sorted = true) const {
			Container result(xx.size());
			if (sorted) {
				for (SizeT i = 0, i_x = 0; i < xx.size(); ++i) {
					const Number evalx = xx[i];
					while (i_x + 1 < _X.size() && evalx >= _X[i_x+1])
						++i_x;
					result[i] = eval(evalx, i_x);
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

		Type type() const {
			return _type;
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

		static std::pair<SizeT, size_t> segment(size_t index) {
			return {4*index, 4*(index+1)};
		}

	private:
		Type _type = Type::Default;
		Container<Number> _X = {};
		Container<Number> _coeffs = {};
	};
}