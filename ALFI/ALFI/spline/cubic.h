#pragma once

#include <algorithm>
#include <iostream>
#include <cmath>
#include <variant>

#include "../config.h"
#include "../util/arrays.h"
#include "../util/misc.h"
#include "../util/spline.h"

namespace alfi::spline {
	template <typename Number = DefaultNumber, template <typename> class Container = DefaultContainer>
	class CubicSpline {
	public:
		struct Types final {
			/**
				Second derivatives at the end points are equal to zero.\n
				Has minimum curvature of all interpolating functions and is most stable.\n
				Equivalent to `FixedSecond{0, 0}`.
			 */
			struct Natural final {};
			/**
				The third derivative is continuous at the second and second-to-last points.\n
				In this way, the second and second-to-last points "cease to be knot points".
			 */
			struct NotAKnot final {};
			/**
				The first and the second derivatives at the end points are equal.
			 */
			struct Periodic final {};
			/**
				The first and the last segments are second degree curves (third derivative equals zero at the ends).\n
				Equivalent to `FixedThird{0, 0}`.
			 */
			struct ParabolicEnds final {};
			/**
				The first derivative equals `df_1` at the first point and `df_n` at the last point.
			 */
			struct Clamped final {
				Clamped(Number df_1, Number df_n) : df_1(std::move(df_1)), df_n(std::move(df_n)) {}
				Number df_1, df_n;
			};
			/**
				The second derivative equals `d2f_1` at the first point and `d2f_n` at the last point.
			 */
			struct FixedSecond final {
				FixedSecond(Number d2f_1, Number d2f_n) : d2f_1(std::move(d2f_1)), d2f_n(std::move(d2f_n)) {}
				Number d2f_1, d2f_n;
			};
			/**
				The third derivative equals `d3f_1` at the first point and `d3f_n` at the last point.\n
				If number of points equals two, the third derivative on the single segment equals (d3f_1+d3f_n)/2.
			 */
			struct FixedThird final {
				FixedThird(Number d3f_1, Number d3f_n) : d3f_1(std::move(d3f_1)), d3f_n(std::move(d3f_n)) {}
				Number d3f_1, d3f_n;
			};
			/**
				A cubic curve is constructed through the first four points. Then each subsequent segment is built sequentially.\n
				This is equivalent to the condition of continuity of the third derivative at the second and third points.\n
				In this way, the second and third points "cease to be knot points".
			 */
			struct NotAKnotStart final {};
			/**
				A cubic curve is constructed through the last four points. Then each previous segment is built sequentially.\n
				This is equivalent to the condition of continuity of the third derivative at the third-to-last and second-to-last points.\n
				In this way, the third-to-last and second-to-last points "cease to be knot points".
			 */
			struct NotAKnotEnd final {};
			/**
				The arithmetic mean_array between NotAKnotStart and NotAKnotEnd.
			 */
			struct SemiNotAKnot final {};
			using Default = Natural;
		};

		using Type = std::variant<typename Types::Natural,
								  typename Types::NotAKnot,
								  typename Types::Periodic,
								  typename Types::ParabolicEnds,
								  typename Types::Clamped,
								  typename Types::FixedSecond,
								  typename Types::FixedThird,
								  typename Types::NotAKnotStart,
								  typename Types::NotAKnotEnd,
								  typename Types::SemiNotAKnot>;

		static Container<Number> compute_coeffs(const auto& X, const auto& Y, Type type = typename Types::Default{}) {
			constexpr auto FUNCTION = __FUNCTION__;

			if (X.size() != Y.size()) {
				std::cerr << "Error in function " << FUNCTION
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

			std::visit(util::misc::overload{
				[&](const typename Types::Natural&) {
					coeffs = compute_coeffs(X, Y, typename Types::FixedSecond{0, 0});
				},
				[&](const typename Types::NotAKnot&) {
					if (n <= 4) {
						coeffs = util::spline::simple_spline<Number,Container>(X, Y, 3);
						return;
					}
					const auto dX = util::arrays::diff(X), dY = util::arrays::diff(Y);
					Container<Number> diag(n - 2), right(n - 2);
					for (SizeT i = 0; i < n - 2; ++i) {
						diag[i] = 2 * (dX[i] + dX[i+1]);
						right[i] = 3 * (dY[i+1]/dX[i+1] - dY[i]/dX[i]);
					}
					const auto C = util::linalg::tridiag_solve(
						dX[0] - dX[1], 2*dX[0] + dX[1], dX[0] / (dX[0]+dX[1]) * right[0],
						2*dX[n-2] + dX[n-3], dX[n-2] - dX[n-3], dX[n-2] / (dX[n-2]+dX[n-3]) * right[n-3],
						dX.begin(), diag.begin(), dX.begin() + 1, right.begin(), n
					);
					for (SizeT i = 0, index = 0; i < n - 1; ++i) {
						coeffs[index++] = (C[i+1] - C[i]) / (3*dX[i]);
						coeffs[index++] = C[i];
						coeffs[index++] = dY[i]/dX[i] - dX[i] * ((2*C[i] + C[i+1]) / 3);
						coeffs[index++] = Y[i];
					}
				},
				[&](const typename Types::Periodic&) {
					if (n <= 2) {
						coeffs = util::spline::simple_spline<Number,Container>(X, Y, 3);
						return;
					}

					const auto dX = util::arrays::diff(X), dY = util::arrays::diff(Y);

					Container<Number> C(n);
					if (n == 3) {
						C[0] = 3 * (dY[0]/dX[0] - dY[1]/dX[1]) / (dX[0] + dX[1]);
						C[1] = -C[0];
						C[2] = C[0];
					} else {
						Container<Number> diag(n - 1), right(n - 1);
						for (SizeT i = 0; i < n - 2; ++i) {
							diag[i] = 2 * (dX[i] + dX[i+1]);
							right[i] = 3 * (dY[i+1]/dX[i+1] - dY[i]/dX[i]);
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

						C[n-1] = right[n-2] / diag[n-2];
						C[n-2] = (right[n-3] - last_col[n-3] * C[n-1]) / diag[n-3];
						for (SizeT i = n - 3; i >= 1; --i) {
							C[i] = (right[i-1] - dX[i]*C[i+1] - last_col[i-1]*C[n-1]) / diag[i-1];
						}
						C[0] = C[n-1];
					}

					for (SizeT i = 0, index = 0; i < n - 1; ++i) {
						coeffs[index++] = (C[i+1] - C[i]) / (3*dX[i]);
						coeffs[index++] = C[i];
						coeffs[index++] = dY[i]/dX[i] - dX[i] * ((2*C[i] + C[i+1]) / 3);
						coeffs[index++] = Y[i];
					}
				},
				[&](const typename Types::ParabolicEnds&) {
					coeffs = compute_coeffs(X, Y, typename Types::FixedThird{0, 0});
				},
				[&](const typename Types::Clamped& c) {
					const auto dX = util::arrays::diff(X), dY = util::arrays::diff(Y);

					Container<Number> C(n);
					if (n == 2) {
						C[0] = (3*dY[0]/dX[0] - c.df_n - 2*c.df_1) / dX[0];
						C[1] = (c.df_1 + 2*c.df_n - 3*dY[0]/dX[0]) / dX[0];
					} else if (n == 3) {
						C[1] = (3*(dY[1]/dX[1] - dY[0]/dX[0]) + (c.df_1 - c.df_n)) / (dX[0] + dX[1]);
						C[0] = (3*(dY[0]/dX[0] - c.df_1)/dX[0] - C[1]) / 2;
						C[2] = (3*(c.df_n - dY[1]/dX[1])/dX[1] - C[1]) / 2;
					} else {
						Container<Number> diag(n - 2), right(n - 2);
						diag[0] = 1.5*dX[0] + 2*dX[1];
						right[0] = 3*dY[1]/dX[1] - 4.5*dY[0]/dX[0] + 1.5*c.df_1;
						for (SizeT i = 1; i < n - 3; ++i) {
							diag[i] = 2 * (dX[i] + dX[i+1]);
							right[i] = 3 * (dY[i+1]/dX[i+1] - dY[i]/dX[i]);
						}
						diag[n-3] = 2*dX[n-3] + 1.5*dX[n-2];
						right[n-3] = 4.5*dY[n-2]/dX[n-2] - 3*dY[n-3]/dX[n-3] - 1.5*c.df_n;

						for (SizeT i = 1; i < n - 2; ++i) {
							const auto m = dX[i] / diag[i-1];
							diag[i] -= m * dX[i];
							right[i] -= m * right[i-1];
						}

						C[n-2] = right[n-3] / diag[n-3];
						for (SizeT i = n - 3; i >= 1; --i) {
							C[i] = (right[i-1] - dX[i] * C[i+1]) / diag[i-1];
						}
						C[n-1] = 1.5*(c.df_n/dX[n-2] - dY[n-2]/dX[n-2]/dX[n-2]) - 0.5*C[n-2];
						C[0] = 1.5*(dY[0]/dX[0]/dX[0] - c.df_1/dX[0]) - 0.5*C[1];
					}

					for (SizeT i = 0, index = 0; i < n - 1; ++i) {
						coeffs[index++] = (C[i+1] - C[i]) / (3*dX[i]);
						coeffs[index++] = C[i];
						coeffs[index++] = dY[i]/dX[i] - dX[i] * ((2*C[i] + C[i+1]) / 3);
						coeffs[index++] = Y[i];
					}
				},
				[&](const typename Types::FixedSecond& s) {
					const auto dX = util::arrays::diff(X), dY = util::arrays::diff(Y);

					Container<Number> C(n);
					C[0] = s.d2f_1 / 2;
					C[n-1] = s.d2f_n / 2;
					if (n > 2) {
						Container<Number> diag(n - 2), right(n - 2);
						for (SizeT i = 0; i < n - 2; ++i) {
							diag[i] = 2 * (dX[i] + dX[i+1]);
							right[i] = 3 * (dY[i+1]/dX[i+1] - dY[i]/dX[i]);
						}
						right[0] -= C[0] * dX[0];
						right[n-3] -= C[n-1] * dX[n-2];
						for (SizeT i = 1; i < n - 2; ++i) {
							const auto m = dX[i] / diag[i-1];
							diag[i] -= m * dX[i];
							right[i] -= m * right[i-1];
						}
						C[n-2] = right[n-3] / diag[n-3];
						for (SizeT i = n - 3; i >= 1; --i) {
							C[i] = (right[i-1] - dX[i] * C[i+1]) / diag[i-1];
						}
					}

					for (SizeT i = 0, index = 0; i < n - 1; ++i) {
						coeffs[index++] = (C[i+1] - C[i]) / (3*dX[i]);
						coeffs[index++] = C[i];
						coeffs[index++] = dY[i]/dX[i] - dX[i] * ((2*C[i] + C[i+1]) / 3);
						coeffs[index++] = Y[i];
					}
				},
				[&](const typename Types::FixedThird& t) {
					const auto dX = util::arrays::diff(X), dY = util::arrays::diff(Y);

					Container<Number> C(n);
					if (n == 2) {
						C[0] = -dX[0] * (t.d3f_1 + t.d3f_n) / 8;
						C[1] = -C[0];
					} else {
						Container<Number> diag(n), right(n);
						diag[0] = 1;
						for (SizeT i = 1; i < n - 1; ++i) {
							diag[i] = 2 * (dX[i-1] + dX[i]);
						}
						diag[n-1] = 1;
						right[0] = -dX[0] * t.d3f_1 / 2;
						for (SizeT i = 1; i < n - 1; ++i) {
							right[i] = 3 * (dY[i]/dX[i] - dY[i-1]/dX[i-1]);
						}
						right[n-1] = dX[n-2] * t.d3f_n / 2;

						const auto m_1 = dX[0] / diag[0];
						diag[1] -= m_1 * (-1);
						right[1] -= m_1 * right[0];
						for (SizeT i = 2; i < n - 1; ++i) {
							const auto m = dX[i-1] / diag[i-1];
							diag[i] -= m * dX[i-1];
							right[i] -= m * right[i-1];
						}
						const auto m_n_1 = (-1) / diag[n-2];
						diag[n-1] -= m_n_1 * dX[n-2];
						right[n-1] -= m_n_1 * right[n-2];

						C[n-1] = right[n-1] / diag[n-1];
						for (SizeT i = n - 2; i >= 1; --i) {
							C[i] = (right[i] - dX[i] * C[i+1]) / diag[i];
						}
						C[0] = right[0] + C[1];
					}

					Container<Number> D(n - 1);
					for (SizeT i = 1; i < D.size() - 1; ++i) {
						D[i] = (C[i+1] - C[i]) / (3*dX[i]);
					}
					if (n == 2) {
						D[0] = (t.d3f_1 + t.d3f_n) / 12;
					} else {
						D[0] = t.d3f_1 / 6;
						D[n-2] = t.d3f_n / 6;
					}
					for (SizeT i = 0, index = 0; i < n - 1; ++i) {
						coeffs[index++] = D[i];
						coeffs[index++] = C[i];
						coeffs[index++] = dY[i]/dX[i] - dX[i] * ((2*C[i] + C[i+1]) / 3);
						coeffs[index++] = Y[i];
					}
				},
				[&](const typename Types::NotAKnotStart&) {
					if (n <= 4) {
						coeffs = util::spline::simple_spline<Number,Container>(X, Y, 3);
						return;
					}

					const auto dX = util::arrays::diff(X), dY = util::arrays::diff(Y);

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
							std::cerr << "Error in function " << FUNCTION << ": Could not compute. Returning an empty array..." << std::endl;
							coeffs = {};
							return;
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
				},
				[&](const typename Types::NotAKnotEnd&) {
					if (n <= 4) {
						coeffs = util::spline::simple_spline<Number,Container>(X, Y, 3);
						return;
					}

					const auto dX = util::arrays::diff(X), dY = util::arrays::diff(Y);

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
							std::cerr << "Error in function " << FUNCTION << ": Could not compute. Returning an empty array..." << std::endl;
							coeffs = {};
							return;
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
				},
				[&](const typename Types::SemiNotAKnot&) {
					coeffs = util::arrays::mean(compute_coeffs(X, Y, typename Types::NotAKnotStart{}), compute_coeffs(X, Y, typename Types::NotAKnotEnd{}));
				},
			}, type);

			return coeffs;
		}

		CubicSpline() = default;

		template <typename ContainerXType>
		CubicSpline(ContainerXType&& X, const auto& Y, Type type = typename Types::Default{}) {
			construct(std::forward<ContainerXType>(X), Y, type);
		}

		CubicSpline(const CubicSpline& other) = default;
		CubicSpline(CubicSpline&& other) noexcept = default;

		CubicSpline& operator=(const CubicSpline& other) = default;
		CubicSpline& operator=(CubicSpline&& other) noexcept = default;

		template <typename ContainerXType>
		void construct(ContainerXType&& X, const auto& Y, Type type = typename Types::Default{}) {
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
		Number eval(Number x, SizeT segment) const {
			if (_coeffs.empty()) {
				return NAN;
			} else if (_coeffs.size() == 1) {
				return _coeffs[0];
			}
			segment = std::clamp(segment, static_cast<SizeT>(0), static_cast<SizeT>(_X.size() - 2));
			x = x - _X[segment];
			return ((_coeffs[4*segment] * x + _coeffs[4*segment+1]) * x + _coeffs[4*segment+2]) * x + _coeffs[4*segment+3];
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
		Type _type = typename Types::Default{};
		Container<Number> _X = {};
		Container<Number> _coeffs = {};
	};
}