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
	class QuadraticSpline {
	public:
		struct Types final {
			/**
				A parabola is constructed through the first three points. Then each subsequent segment is built sequentially.\n
				This is equivalent to the condition of continuity of the second derivative at the second point.\n
				In this way, the second point "ceases to be a knot point".
			 */
			struct NotAKnotStart final {};
			/**
				A parabola is constructed through the last three points. Then each previous segment is built sequentially.\n
				This is equivalent to the condition of continuity of the second derivative at the second-to-last point.\n
				In this way, the second-to-last point "ceases to be a knot point".
			 */
			struct NotAKnotEnd final {};
			/**
				The arithmetic mean_array between NotAKnotStart and NotAKnotEnd.
			 */
			struct SemiNotAKnot final {};
			/**
				The second derivative at the first point is equal to zero.
			 */
			struct NaturalStart final {};
			/**
				The second derivative at the last point is equal to zero.
			 */
			struct NaturalEnd final {};
			/**
				The arithmetic mean_array between NaturalStart and NaturalEnd.
			 */
			struct SemiNatural final {};
			/**
				The arithmetic mean_array between SemiNotAKnot and SemiNatural.
			 */
			struct SemiSemi final {};
			using Default = SemiNotAKnot;
		};

		using Type = std::variant<typename Types::NotAKnotStart,
								  typename Types::NotAKnotEnd,
								  typename Types::SemiNotAKnot,
								  typename Types::NaturalStart,
								  typename Types::NaturalEnd,
								  typename Types::SemiNatural,
								  typename Types::SemiSemi>;

		static Container<Number> compute_coeffs(const auto& X, const auto& Y, Type type = typename Types::Default{}) {
			if (X.size() != Y.size()) {
				std::cerr << "Error in function " << __FUNCTION__
						  << ": Vectors X (of size " << X.size()
						  << ") and Y (of size " << Y.size()
						  << ") are not the same size. Returning an empty array..." << std::endl;
				return {};
			}

			const SizeT n = X.size();

			if (n <= 1) {
				return util::spline::simple_spline<Number,Container>(X, Y, 2);
			}

			Container<Number> coeffs(3 * (n - 1));

			std::visit(util::misc::overload{
				[&](const typename Types::NotAKnotStart&) {
					if (n <= 3) {
						coeffs = util::spline::simple_spline<Number,Container>(X, Y, 2);
						return;
					}
					{ // first three points
						const Number h1 = X[1] - X[0], h2 = X[2] - X[1];
						const Number d1 = (Y[1] - Y[0]) / h1, d2 = (Y[2] - Y[1]) / h2;
						coeffs[0] = (d2 - d1) / (h1 + h2); // a1
						coeffs[1] = d1 - h1 * coeffs[0]; // b1
						coeffs[2] = Y[0]; // c1
						coeffs[3] = coeffs[0]; // a2
						coeffs[4] = d1 + h1 * coeffs[0]; // b2
						coeffs[5] = Y[1]; // c2
					}
					for (SizeT i = 2; i < n - 1; ++i) {
						const Number hi = X[i+1] - X[i], hi_1 = X[i] - X[i-1];
						const Number d = 2 * coeffs[3*(i-1)] * hi_1 + coeffs[3*(i-1)+1];
						coeffs[3*i] = ((Y[i+1] - Y[i])/hi - d) / hi;
						coeffs[3*i+1] = d;
						coeffs[3*i+2] = Y[i];
					}
				},
				[&](const typename Types::NotAKnotEnd&) {
					if (n <= 3) {
						coeffs = util::spline::simple_spline<Number,Container>(X, Y, 2);
						return;
					}
					{ // last three points
						const Number hn_2 = X[n-2] - X[n-3], hn_1 = X[n-1] - X[n-2];
						const Number dn_2 = (Y[n-2] - Y[n-3]) / hn_2, dn_1 = (Y[n-1] - Y[n-2]) / hn_1;
						coeffs[3*(n-3)] = (dn_1 - dn_2) / (hn_2 + hn_1); // an-2
						coeffs[3*(n-3)+1] = dn_2 - hn_2 * coeffs[3*(n-3)]; // bn-2
						coeffs[3*(n-3)+2] = Y[n-3]; // cn-2
						coeffs[3*(n-2)] = coeffs[3*(n-3)]; // an-1
						coeffs[3*(n-2)+1] = dn_2 + hn_2 * coeffs[3*(n-3)]; // bn-1
						coeffs[3*(n-2)+2] = Y[n-2]; // cn-1
					}
					for (SizeT iter = 0; iter <= n - 4; ++iter) {
						const SizeT i = n - 4 - iter;
						const Number hi = X[i+1] - X[i];
						coeffs[3*i] = (coeffs[3*(i+1)+1] - (Y[i+1] - Y[i])/hi) / hi;
						coeffs[3*i+1] = coeffs[3*(i+1)+1] - 2*coeffs[3*i]*hi;
						// another way: coeffs[3*i+1] = (Y[i+1] - Y[i])/hi - coeffs[3*i]*hi;
						coeffs[3*i+2] = Y[i];
					}
				},
				[&](const typename Types::SemiNotAKnot&) {
					coeffs = util::arrays::mean(compute_coeffs(X, Y, typename Types::NotAKnotStart{}), compute_coeffs(X, Y, typename Types::NotAKnotEnd{}));
				},
				[&](const typename Types::NaturalStart&) {
					if (n <= 2) {
						coeffs = util::spline::simple_spline<Number,Container>(X, Y, 2);
						return;
					}
					coeffs[0] = 0; // a1 = 0
					coeffs[1] = (Y[1] - Y[0]) / (X[1] - X[0]); // b1
					coeffs[2] = Y[0]; // c1
					for (SizeT i = 1; i < n - 1; ++i) {
						coeffs[3*i+1] = 2 * (Y[i] - Y[i-1]) / (X[i] - X[i-1]) - coeffs[3*(i-1)+1];
						coeffs[3*i] = ((Y[i+1] - Y[i]) / (X[i+1] - X[i]) - coeffs[3*i+1]) / (X[i+1] - X[i]);
						// another way: coeffs[3*i] = (coeffs[3*(i+1)+1] - coeffs[3*i+1]) / 2 / (X[i+1] - X[i]);
						coeffs[3*i+2] = Y[i];
					}
				},
				[&](const typename Types::NaturalEnd&) {
					if (n <= 2) {
						coeffs = util::spline::simple_spline<Number,Container>(X, Y, 2);
						return;
					}
					coeffs[3*(n-2)] = 0; // an-1 = 0
					coeffs[3*(n-2)+1] = (Y[n-1] - Y[n-2]) / (X[n-1] - X[n-2]); // bn-1
					coeffs[3*(n-2)+2] = Y[n-2]; // cn-1
					for (SizeT iter = 0; iter <= n - 3; ++iter) {
						const SizeT i = n - 3 - iter;
						coeffs[3*i+1] = 2 * (Y[i+1] - Y[i]) / (X[i+1] - X[i]) - coeffs[3*(i+1)+1];
						// another way: coeffs[3*i] = ((Y[i+1] - Y[i]) / (X[i+1] - X[i]) - coeffs[3*i+1]) / (X[i+1] - X[i]);
						coeffs[3*i] = (coeffs[3*(i+1)+1] - coeffs[3*i+1]) / 2 / (X[i+1] - X[i]);
						coeffs[3*i+2] = Y[i];
					}
				},
				[&](const typename Types::SemiNatural&) {
					coeffs = util::arrays::mean(compute_coeffs(X, Y, typename Types::NaturalStart{}), compute_coeffs(X, Y, typename Types::NaturalEnd{}));
				},
				[&](const typename Types::SemiSemi&) {
					coeffs = util::arrays::mean(compute_coeffs(X, Y, typename Types::SemiNotAKnot{}), compute_coeffs(X, Y, typename Types::SemiNatural{}));
				},
			}, type);

			return coeffs;
		}

		QuadraticSpline() = default;

		template <typename ContainerXType>
		QuadraticSpline(ContainerXType&& X, const auto& Y, Type type = typename Types::Default{}) {
			construct(std::forward<ContainerXType>(X), Y, type);
		}

		QuadraticSpline(const QuadraticSpline& other) = default;
		QuadraticSpline(QuadraticSpline&& other) noexcept = default;

		QuadraticSpline& operator=(const QuadraticSpline& other) = default;
		QuadraticSpline& operator=(QuadraticSpline&& other) noexcept = default;

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
			return (_coeffs[3*segment] * x + _coeffs[3*segment+1]) * x + _coeffs[3*segment+2];
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

		static std::pair<SizeT, SizeT> segment(SizeT index) {
			return {3*index, 3*(index+1)};
		}

	private:
		Type _type = typename Types::Default{};
		Container<Number> _X = {};
		Container<Number> _coeffs = {};
	};
}