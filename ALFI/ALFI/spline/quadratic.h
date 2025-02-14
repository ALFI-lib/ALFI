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
				The arithmetic mean of NotAKnotStart and NotAKnotEnd.
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
				The arithmetic mean of NaturalStart and NaturalEnd.
			 */
			struct SemiNatural final {};
			/**
				The arithmetic mean of SemiNotAKnot and SemiNatural.
			 */
			struct SemiSemi final {};
			/**
				The first derivative equals `df` at the first point.
			 */
			struct ClampedStart final {
				explicit ClampedStart(Number df) : df(std::move(df)) {}
				Number df;
			};
			/**
				The first derivative equals `df` at the last point.
			 */
			struct ClampedEnd final {
				explicit ClampedEnd(Number df) : df(std::move(df)) {}
				Number df;
			};
			/**
				The arithmetic mean of ClampedStart and ClampedEnd.
			 */
			struct SemiClamped final {
				SemiClamped(Number df_1, Number df_n) : df_1(std::move(df_1)), df_n(std::move(df_n)) {}
				Number df_1, df_n;
			};
			/**
				The second derivative equals `d2f` on the first segment.
			 */
			struct FixedSecondStart final {
				explicit FixedSecondStart(Number d2f) : d2f(std::move(d2f)) {}
				Number d2f;
			};
			/**
				The second derivative equals `d2f` on the last segment.
			 */
			struct FixedSecondEnd final {
				explicit FixedSecondEnd(Number d2f) : d2f(std::move(d2f)) {}
				Number d2f;
			};
			/**
				The arithmetic mean of FixedSecondStart and FixedSecondEnd.
			 */
			struct SemiFixedSecond final {
				SemiFixedSecond(Number d2f_1, Number d2f_n) : d2f_1(std::move(d2f_1)), d2f_n(std::move(d2f_n)) {}
				Number d2f_1, d2f_n;
			};
			/**
				The first derivative equals `df` at the point with index `point_idx`.
			 */
			struct Clamped final {
				Clamped(SizeT point_idx, Number df) : point_idx(std::move(point_idx)), df(std::move(df)) {}
				SizeT point_idx;
				Number df;
			};
			/**
				The second derivative equals `d2f` on the segment with index `segment_idx`.
			 */
			struct FixedSecond final {
				FixedSecond(SizeT segment_idx, Number d2f) : segment_idx(std::move(segment_idx)), d2f(std::move(d2f)) {}
				SizeT segment_idx;
				Number d2f;
			};
			/**
				The second derivative is continuous at the point with index `point_idx`.\n
				This also means, that the second derivative is the same on the segments with indices `point_idx-1` and `point_idx`.\n
				If number of points is less than or equal to 2, then `point_idx` is ignored.
			 */
			struct NotAKnot final {
				explicit NotAKnot(SizeT point_idx) : point_idx(std::move(point_idx)) {}
				SizeT point_idx;
			};
			using Default = SemiNotAKnot;
		};

		using Type = std::variant<typename Types::NotAKnotStart,
								  typename Types::NotAKnotEnd,
								  typename Types::SemiNotAKnot,
								  typename Types::NaturalStart,
								  typename Types::NaturalEnd,
								  typename Types::SemiNatural,
								  typename Types::SemiSemi,
								  typename Types::ClampedStart,
								  typename Types::ClampedEnd,
								  typename Types::SemiClamped,
								  typename Types::FixedSecondStart,
								  typename Types::FixedSecondEnd,
								  typename Types::SemiFixedSecond,
								  typename Types::Clamped,
								  typename Types::FixedSecond,
								  typename Types::NotAKnot>;

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

			if (std::holds_alternative<typename Types::NotAKnotStart>(type)) {
				return compute_coeffs(X, Y, typename Types::NotAKnot(1));
			} else if (std::holds_alternative<typename Types::NotAKnotEnd>(type)) {
				return compute_coeffs(X, Y, typename Types::NotAKnot(n - 2));
			} else if (std::holds_alternative<typename Types::SemiNotAKnot>(type)) {
				return util::arrays::mean(compute_coeffs(X, Y, typename Types::NotAKnot{1}), compute_coeffs(X, Y, typename Types::NotAKnot{n - 2}));
			} else if (std::holds_alternative<typename Types::NaturalStart>(type)) {
				return compute_coeffs(X, Y, typename Types::FixedSecond(0, 0));
			} else if (std::holds_alternative<typename Types::NaturalEnd>(type)) {
				return compute_coeffs(X, Y, typename Types::FixedSecond(n - 2, 0));
			} else if (std::holds_alternative<typename Types::SemiNatural>(type)) {
				return util::arrays::mean(compute_coeffs(X, Y, typename Types::NaturalStart{}), compute_coeffs(X, Y, typename Types::NaturalEnd{}));
			} else if (std::holds_alternative<typename Types::SemiSemi>(type)) {
				return util::arrays::mean(compute_coeffs(X, Y, typename Types::SemiNotAKnot{}), compute_coeffs(X, Y, typename Types::SemiNatural{}));
			} else if (const auto* cs = std::get_if<typename Types::ClampedStart>(&type)) {
				return compute_coeffs(X, Y, typename Types::Clamped(0, cs->df));
			} else if (const auto* ce = std::get_if<typename Types::ClampedEnd>(&type)) {
				return compute_coeffs(X, Y, typename Types::Clamped(n - 1, ce->df));
			} else if (const auto* sc = std::get_if<typename Types::SemiClamped>(&type)) {
				return util::arrays::mean(compute_coeffs(X, Y, typename Types::Clamped(0, sc->df_1)), compute_coeffs(X, Y, typename Types::Clamped(n - 1, sc->df_n)));
			} else if (const auto* fss = std::get_if<typename Types::FixedSecondStart>(&type)) {
				return compute_coeffs(X, Y, typename Types::FixedSecond(0, fss->d2f));
			} else if (const auto* fse = std::get_if<typename Types::FixedSecondEnd>(&type)) {
				return compute_coeffs(X, Y, typename Types::FixedSecond(n - 2, fse->d2f));
			} else if (const auto* sfs = std::get_if<typename Types::SemiFixedSecond>(&type)) {
				return util::arrays::mean(compute_coeffs(X, Y, typename Types::FixedSecond(0, sfs->d2f_1)), compute_coeffs(X, Y, typename Types::FixedSecond(n - 2, sfs->d2f_n)));
			}

			Container<Number> coeffs(3 * (n - 1));

			SizeT i1 = 0, i2 = n - 1;

			if (const auto* clamped = std::get_if<typename Types::Clamped>(&type)) {
				const auto i = clamped->point_idx;
				if (i < 0 || i >= n) {
					std::cerr << "Error in function " << __FUNCTION__
							  << ": point index for type 'Clamped' is out of bounds:"
							  << "expected to be non-negative and less than " << n << ", but got " << i
							  << ". Returning an empty array..." << std::endl;
					return {};
				}
				if (i > 0) {
					const auto dx1 = X[i] - X[i-1], dy1 = Y[i] - Y[i-1];
					coeffs[3*(i-1)+2] = Y[i-1];
					coeffs[3*(i-1)+1] = 2*dy1/dx1 - clamped->df;
					coeffs[3*(i-1)] = dy1/dx1/dx1 - coeffs[3*(i-1)+1]/dx1;
					i1 = i - 1;
				}
				if (i < n - 1) {
					const auto dx = X[i+1] - X[i], dy = Y[i+1] - Y[i];
					coeffs[3*i+2] = Y[i];
					coeffs[3*i+1] = clamped->df;
					coeffs[3*i] = dy/dx/dx - coeffs[3*i+1]/dx;
					i2 = i;
				}
			} else if (const auto* fixed_second = std::get_if<typename Types::FixedSecond>(&type)) {
				const auto i = fixed_second->segment_idx;
				if (i < 0 || i >= n - 1) {
					std::cerr << "Error in function " << __FUNCTION__
							  << ": point index for type 'FixedSecond' is out of bounds:"
							  << "expected to be non-negative and less than " << n - 1 << ", but got " << i
							  << ". Returning an empty array..." << std::endl;
					return {};
				}
				coeffs[3*i] = fixed_second->d2f / 2;
				coeffs[3*i+1] = (Y[i+1]-Y[i])/(X[i+1]-X[i]) - coeffs[3*i]*(X[i+1]-X[i]);
				coeffs[3*i+2] = Y[i];
				i1 = i2 = i;
			} else if (const auto* not_a_knot = std::get_if<typename Types::NotAKnot>(&type)) {
				if (n <= 2) {
					return util::spline::simple_spline<Number,Container>(X, Y, 2);
				}
				const auto i = not_a_knot->point_idx;
				if (i < 1 || i >= n - 1) {
					std::cerr << "Error in function " << __FUNCTION__
							  << ": point index for type 'NotAKnot' is out of bounds:"
							  << "expected to be positive and less than " << n - 1 << ", but got " << i
							  << ". Returning an empty array..." << std::endl;
					return {};
				}
				const auto dx1 = X[i] - X[i-1], dy1 = Y[i] - Y[i-1];
				const auto dx = X[i+1] - X[i], dy = Y[i+1] - Y[i];
				coeffs[3*(i-1)+2] = Y[i-1];
				coeffs[3*i+2] = Y[i];
				coeffs[3*(i-1)] = coeffs[3*i] = (dy/dx - dy1/dx1) / (dx1 + dx);
				coeffs[3*(i-1)+1] = dy1/dx1 - coeffs[3*(i-1)]*dx1;
				coeffs[3*i+1] = dy/dx - coeffs[3*i]*dx;
				i1 = i - 1;
				i2 = i;
			} else {
				std::cerr << "Error in function " << __FUNCTION__ << ": Unknown type. This should not have happened."
						  << " Please report this to the developers. Returning an empty array..." << std::endl;
				return {};
			}

			for (SizeT iter = 0; iter < i1; ++iter) {
				const auto i = i1 - 1 - iter;
				const auto dx = X[i+1] - X[i], dy = Y[i+1] - Y[i];
				coeffs[3*i] = coeffs[3*(i+1)+1]/dx - dy/dx/dx;
				coeffs[3*i+1] = 2*dy/dx - coeffs[3*(i+1)+1];
				coeffs[3*i+2] = Y[i];
			}

			for (SizeT i = i2 + 1; i < n - 1; ++i) {
				const auto dx = X[i+1] - X[i], dy = Y[i+1] - Y[i];
				coeffs[3*i+2] = Y[i];
				coeffs[3*i+1] = coeffs[3*(i-1)+1] + 2*coeffs[3*(i-1)]*(X[i]-X[i-1]);
				coeffs[3*i] = dy/dx/dx - coeffs[3*i+1]/dx;
			}

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

		static std::pair<SizeT, SizeT> segment(SizeT index) {
			return {3*index, 3*(index+1)};
		}

	private:
		Type _type = typename Types::Default{};
		Container<Number> _X = {};
		Container<Number> _coeffs = {};
	};
}