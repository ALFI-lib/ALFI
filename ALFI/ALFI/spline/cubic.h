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
	template <typename Number = DefaultNumber, template <typename, typename...> class Container = DefaultContainer>
	class CubicSpline {
	public:
		struct Conditions final {
			struct Clamped final {
				Clamped(SizeT point_idx, Number df) : point_idx(std::move(point_idx)), df(std::move(df)) {}
				SizeT point_idx;
				Number df;
			};
			struct FixedSecond final {
				FixedSecond(SizeT point_idx, Number d2f) : point_idx(std::move(point_idx)), d2f(std::move(d2f)) {}
				SizeT point_idx;
				Number d2f;
			};
			struct FixedThird final {
				FixedThird(SizeT segment_idx, Number d3f) : segment_idx(std::move(segment_idx)), d3f(std::move(d3f)) {}
				SizeT segment_idx;
				Number d3f;
			};
			struct NotAKnot final {
				explicit NotAKnot(SizeT point_idx) : point_idx(std::move(point_idx)) {}
				SizeT point_idx;
			};
		};

		using Condition = std::variant<typename Conditions::Clamped,
								  typename Conditions::FixedSecond,
								  typename Conditions::FixedThird,
								  typename Conditions::NotAKnot>;

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
				The arithmetic mean of NotAKnotStart and NotAKnotEnd.
			 */
			struct SemiNotAKnot final {};
			struct Custom final {
				Custom(Condition condition1, Condition condition2) : cond1(std::move(condition1)), cond2(std::move(condition2)) {}
				Condition cond1, cond2;
			};
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

		static Container<Number> compute_coeffs(const Container<Number>& X, const Container<Number>& Y, Type type = typename Types::Default{}) {
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

			if (std::holds_alternative<typename Types::Natural>(type)) {
				return compute_coeffs(X, Y, typename Types::Custom{typename Conditions::FixedSecond{0, 0}, typename Conditions::FixedSecond{n - 1, 0}});
			} else if (std::holds_alternative<typename Types::NotAKnot>(type)) {
				if (n == 2) {
					return compute_coeffs(X, Y, typename Types::Custom{typename Conditions::NotAKnot{0}, typename Conditions::NotAKnot{1}});
				}
				return compute_coeffs(X, Y, typename Types::Custom{typename Conditions::NotAKnot{1}, typename Conditions::NotAKnot{n - 2}});
			} else if (std::holds_alternative<typename Types::ParabolicEnds>(type)) {
				return compute_coeffs(X, Y, typename Types::Custom{typename Conditions::FixedThird{0, 0}, typename Conditions::FixedThird{n - 2, 0}});
			} else if (const auto* c = std::get_if<typename Types::Clamped>(&type)) {
				return compute_coeffs(X, Y, typename Types::Custom{typename Conditions::Clamped{0, c->df_1}, typename Conditions::Clamped{n - 1, c->df_n}});
			} else if (const auto* fs = std::get_if<typename Types::FixedSecond>(&type)) {
				return compute_coeffs(X, Y, typename Types::Custom{typename Conditions::FixedSecond{0, fs->d2f_1}, typename Conditions::FixedSecond{n - 1, fs->d2f_n}});
			} else if (const auto* ft = std::get_if<typename Types::FixedThird>(&type)) {
				return compute_coeffs(X, Y, typename Types::Custom{typename Conditions::FixedThird{0, ft->d3f_1}, typename Conditions::FixedThird{n - 2, ft->d3f_n}});
			} else if (std::holds_alternative<typename Types::NotAKnotStart>(type)) {
				if (n == 2) {
					return compute_coeffs(X, Y, typename Types::Custom{typename Conditions::NotAKnot{0}, typename Conditions::NotAKnot{1}});
				}
				return compute_coeffs(X, Y, typename Types::Custom{typename Conditions::NotAKnot{1}, typename Conditions::NotAKnot{2}});
			} else if (std::holds_alternative<typename Types::NotAKnotEnd>(type)) {
				if (n == 2) {
					return compute_coeffs(X, Y, typename Types::Custom{typename Conditions::NotAKnot{0}, typename Conditions::NotAKnot{1}});
				}
				return compute_coeffs(X, Y, typename Types::Custom{typename Conditions::NotAKnot{n - 3}, typename Conditions::NotAKnot{n - 2}});
			} else if (std::holds_alternative<typename Types::SemiNotAKnot>(type)) {
				return util::arrays::mean(compute_coeffs(X, Y, typename Types::NotAKnotStart{}), compute_coeffs(X, Y, typename Types::NotAKnotEnd{}));
			}

			Container<Number> coeffs(4 * (n - 1));

			SizeT i1 = 0, i2 = n - 1;

			const auto dX = util::arrays::diff(X), dY = util::arrays::diff(Y);

			if (std::holds_alternative<typename Types::Periodic>(type)) {
				if (n == 2) {
					return util::spline::simple_spline<Number,Container>(X, Y, 3);
				}

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
			} else if (auto* custom = std::get_if<typename Types::Custom>(&type)) {
				if (std::holds_alternative<typename Conditions::NotAKnot>(custom->cond1)
						&& std::holds_alternative<typename Conditions::NotAKnot>(custom->cond2)
						&& n <= 4) {
					return util::spline::simple_spline<Number,Container>(X, Y, 3);
						}
				if (const auto* c = std::get_if<typename Conditions::Clamped>(&custom->cond1)) {
					const auto i = c->point_idx;

				} else if (const auto* fs = std::get_if<typename Conditions::FixedSecond>(&custom->cond1)) {

				} else if (const auto* ft = std::get_if<typename Conditions::FixedThird>(&custom->cond1)) {

				} else if (const auto* nak = std::get_if<typename Conditions::NotAKnot>(&custom->cond1)) {

				} else {
					std::cerr << "Error in function " << __FUNCTION__ << ": Unknown first condition of 'Custom' type. This should not have happened."
							  << " Please report this to the developers. Returning an empty array..." << std::endl;
					return {};
				}
				if (const auto* c = std::get_if<typename Conditions::Clamped>(&custom->cond2)) {

				} else if (const auto* fs = std::get_if<typename Conditions::FixedSecond>(&custom->cond2)) {

				} else if (const auto* ft = std::get_if<typename Conditions::FixedThird>(&custom->cond2)) {

				} else if (const auto* nak = std::get_if<typename Conditions::NotAKnot>(&custom->cond2)) {

				} else {
					std::cerr << "Error in function " << __FUNCTION__ << ": Unknown second condition of 'Custom' type. This should not have happened."
							  << " Please report this to the developers. Returning an empty array..." << std::endl;
					return {};
				}
				if (auto c1 = std::get_if<typename Conditions::Clamped>(&custom->cond1),
						c2 = std::get_if<typename Conditions::Clamped>(&custom->cond2);
						c1 && c2) {
					if (c2->point_idx < c1->point_idx) {
						std::swap(custom->cond1, custom->cond2);
					}
				}
				if (auto fs1 = std::get_if<typename Conditions::FixedSecond>(&custom->cond1),
						fs2 = std::get_if<typename Conditions::FixedSecond>(&custom->cond2);
						fs1 && fs2) {
					if (fs2->point_idx < fs1->point_idx) {
						std::swap(custom->cond1, custom->cond2);
					}
				}
				if (auto ft1 = std::get_if<typename Conditions::FixedThird>(&custom->cond1),
						ft2 = std::get_if<typename Conditions::FixedThird>(&custom->cond2);
						ft1 && ft2) {
					if (ft2->segment_idx < ft1->segment_idx) {
						std::swap(custom->cond1, custom->cond2);
					}
				}
				if (auto n1 = std::get_if<typename Conditions::NotAKnot>(&custom->cond1),
						n2 = std::get_if<typename Conditions::NotAKnot>(&custom->cond2);
						n1 && n2) {
					if (n2->point_idx < n1->point_idx) {
						std::swap(custom->cond1, custom->cond2);
					}
				}

				if (std::holds_alternative<typename Conditions::FixedThird>(custom->cond1)
					&& std::holds_alternative<typename Conditions::FixedThird>(custom->cond2)
					&& n == 2) {
					// проверить индексы и приравнять третьи производные их среднему арифметическому
				}

				Container<Number> lower(n), diag(n), upper(n), right(n);
				for (SizeT i = 1; i < n - 1; ++i) {
					lower[i] = dX[i-1];
					diag[i] = 2 * (dX[i-1] + dX[i]);
					upper[i] = dX[i];
					right[i] = 3 * (dY[i]/dX[i] - dY[i-1]/dX[i-1]);
				}
				lower[0] = NAN;
				upper[n-1] = NAN;

				// уникальные условия

				const auto C = util::linalg::tridiag_solve<Number,Container>(
					std::move(lower), std::move(diag), std::move(upper), std::move(right));

				for (SizeT i = 0, index = 0; i < n - 1; ++i) {
					coeffs[index++] = (C[i+1] - C[i]) / (3*dX[i]);
					coeffs[index++] = C[i];
					coeffs[index++] = dY[i]/dX[i] - dX[i] * ((2*C[i] + C[i+1]) / 3);
					coeffs[index++] = Y[i];
				}

				// здесь делаю дополнительные правки

			} else {
				std::cerr << "Error in function " << __FUNCTION__ << ": Unknown type. This should not have happened."
						  << " Please report this to the developers. Returning an empty array..." << std::endl;
				return {};
			}

			for (SizeT iter = 0; iter < i1; ++iter) {
				const auto i = i1 - 1 - iter;
				coeffs[4*i] = (coeffs[4*(i+1)+1] - coeffs[4*(i+1)+2]/dX[i] + dY[i]/pow(dX[i], 2)) / dX[i];
				coeffs[4*i+1] = coeffs[4*(i+1)+1] - 3 * coeffs[4*i] * dX[i];
				coeffs[4*i+2] = dY[i]/dX[i] - coeffs[4*i]*pow(dX[i], 2) - coeffs[4*i+1]*dX[i];
				coeffs[4*i+3] = Y[i];
			}

			for (SizeT i = i2 + 1; i < n - 1; ++i) {
				coeffs[4*i+3] = Y[i];
				coeffs[4*i+2] = 3 * coeffs[4*(i-1)] * pow(dX[i-1], 2) + 2 * coeffs[4*(i-1)+1] * dX[i-1] + coeffs[4*(i-1)+2];
				coeffs[4*i+1] = 3 * coeffs[4*(i-1)] * dX[i-1] + coeffs[4*(i-1)+1];
				coeffs[4*i] = (dY[i]/dX[i] - coeffs[4*i+1] * dX[i] - coeffs[4*i+2]) / pow(dX[i], 2);
			}

			std::visit(util::misc::overload{
				[&](const typename Types::NotAKnot&) {
					// diag[0] = dX[0] - dX[1];
					// upper[0] = 2*dX[0] + dX[1];
					// right[0] = dX[0] / (dX[0]+dX[1]) * right[1];
					// lower[n-1] = 2*dX[n-2] + dX[n-3];
					// diag[n-1] = dX[n-2] - dX[n-3];
					// right[n-1] = dX[n-2] / (dX[n-2]+dX[n-3]) * right[n-2];
					// const auto C = util::linalg::tridiag_solve<Number,Container>(
					// 	std::move(lower), std::move(diag), std::move(upper), std::move(right));
					for (SizeT i = 0, index = 0; i < n - 1; ++i) {
						coeffs[index++] = (C[i+1] - C[i]) / (3*dX[i]);
						coeffs[index++] = C[i];
						coeffs[index++] = dY[i]/dX[i] - dX[i] * ((2*C[i] + C[i+1]) / 3);
						coeffs[index++] = Y[i];
					}
				},
				[&](const typename Types::Clamped& c) {
					// diag[0] = 2*dX[0];
					// upper[0] = dX[0];
					// right[0] = 3 * (dY[0]/dX[0] - c.df_1);
					// lower[n-1] = dX[n-2];
					// diag[n-1] = 2*dX[n-2];
					// right[n-1] = 3 * (c.df_n - dY[n-2]/dX[n-2]);
					// const auto C = util::linalg::tridiag_solve<Number,Container>(
					// 	std::move(lower), std::move(diag), std::move(upper), std::move(right));
					for (SizeT i = 0, index = 0; i < n - 1; ++i) {
						coeffs[index++] = (C[i+1] - C[i]) / (3*dX[i]);
						coeffs[index++] = C[i];
						coeffs[index++] = dY[i]/dX[i] - dX[i] * ((2*C[i] + C[i+1]) / 3);
						coeffs[index++] = Y[i];
					}
				},
				[&](const typename Types::FixedSecond& s) {
					// diag[0] = 1;
					// upper[0] = 0;
					// right[0] = s.d2f_1 / 2;
					// lower[n-1] = 0;
					// diag[n-1] = 1;
					// right[n-1] = s.d2f_n / 2;
					// const auto C = util::linalg::tridiag_solve<Number,Container>(
					// 	std::move(lower), std::move(diag), std::move(upper), std::move(right));
					for (SizeT i = 0, index = 0; i < n - 1; ++i) {
						coeffs[index++] = (C[i+1] - C[i]) / (3*dX[i]);
						coeffs[index++] = C[i];
						coeffs[index++] = dY[i]/dX[i] - dX[i] * ((2*C[i] + C[i+1]) / 3);
						coeffs[index++] = Y[i];
					}
				},
				[&](const typename Types::FixedThird& t) {
					Container<Number> C(n);
					if (n == 2) {
						C[0] = -dX[0] * (t.d3f_1 + t.d3f_n) / 8;
						C[1] = -C[0];
					} else {
						// diag[0] = 1;
						// upper[0] = -1;
						// right[0] = -dX[0] * t.d3f_1 / 2;
						// lower[n-1] = -1;
						// diag[n-1] = 1;
						// right[n-1] = dX[n-2] * t.d3f_n / 2;
						// C = util::linalg::tridiag_solve<Number,Container>(
						// 	std::move(lower), std::move(diag), std::move(upper), std::move(right));
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
				}
			}, type);

			return coeffs;
		}

		CubicSpline() = default;

		template <typename ContainerXType>
		CubicSpline(ContainerXType&& X, const Container<Number>& Y, Type type = typename Types::Default{}) {
			construct(std::forward<ContainerXType>(X), Y, type);
		}

		CubicSpline(const CubicSpline& other) = default;
		CubicSpline(CubicSpline&& other) noexcept = default;

		CubicSpline& operator=(const CubicSpline& other) = default;
		CubicSpline& operator=(CubicSpline&& other) noexcept = default;

		template <typename ContainerXType>
		void construct(ContainerXType&& X, const Container<Number>& Y, Type type = typename Types::Default{}) {
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
			Container<Number> result(xx.size());
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