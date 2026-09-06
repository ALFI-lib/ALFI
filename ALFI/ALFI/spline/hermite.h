#pragma once

#include <algorithm>
#include <iostream>
#include <cmath>
#include <variant>

#include "../config.h"
#include "../util/misc.h"
#include "../util/spline.h"

namespace alfi::spline {
	template <typename Number = DefaultNumber, template <typename, typename...> class Container = DefaultContainer>
	class HermiteSpline {
	public:
		struct Types final {
			/**
				Known as "Finite difference".\n
				One-sided secant derivatives at the end points and arithmetic mean of adjacent secants at interior points.
			 */
			struct Classic final {};
			/**
				Cardinal spline, sometimes called a canonical spline, with tension parameter `c`.\n
				The derivative at each point is `(1-c)` times the slope of the chord between adjacent points.
			 */
			struct Cardinal final {
				explicit Cardinal(Number c) : c(std::move(c)) {}
				Number c;
			};
			/**
				Cardinal spline with tension parameter `c = 0`.\n
				This implementation is suitable for non-uniformly distributed points.
			 */
			struct CatmullRom final {};
			/**
				Akima interpolation, using weighted averages of neighboring secant derivatives.
			 */
			struct Akima final {};
			/**
				Modified Akima interpolation with additional weights that reduce overshoot on flat data.
			 */
			struct ModifiedAkima final {};
			/**
				Monotonicity-preserving cubic interpolation using Steffen's derivative estimates.
			 */
			struct Steffen final {};
			/**
				Zero first derivatives at all points.\n
				Each segment is monotonic between its end points.
			 */
			struct Zero final {};
			/**
				Explicitly specified first derivatives at all points.
			 */
			struct Explicit final {
				explicit Explicit(Container<Number> derivatives) : derivatives(std::move(derivatives)) {}
				Container<Number> derivatives;
			};
			using Default = Classic;
		};

		using Type = std::variant<typename Types::Classic,
								  typename Types::Cardinal,
								  typename Types::CatmullRom,
								  typename Types::Akima,
								  typename Types::ModifiedAkima,
								  typename Types::Steffen,
								  typename Types::Zero,
								  typename Types::Explicit>;

		struct Boundaries final {
			/**
				Use the boundary behavior defined by the spline type.
			 */
			struct Inherit final {};
			/**
				The first derivative at each end point equals the slope of the adjacent segment.
			 */
			struct Linear final {};
			/**
				The first derivative at each end point is obtained from a quadratic interpolating polynomial.
			 */
			struct Quadratic final {};
			/**
				The first derivative at each end point is obtained from a cubic interpolating polynomial.
			 */
			struct Cubic final {};
			/**
				The first derivative at each end point is obtained from an interpolating polynomial of the specified degree.
			 */
			struct Polynomial final {
				explicit Polynomial(SizeT degree) : degree(std::move(degree)) {}
				SizeT degree;
			};
			/**
				The first derivative equals `left` at the first point and `right` at the last point.
			 */
			struct Clamped final {
				explicit Clamped(Number left, Number right) : left(std::move(left)), right(std::move(right)) {}
				Number left, right;
			};
			/**
				The first derivative is periodic: the derivatives at the first and last points are equal,
				and periodic neighbors are used when estimating boundary derivatives.

				Has no effect when used with `Types::Explicit`.
			 */
			struct Periodic final {};
			using Default = Inherit;
		};

		using BoundariesType = std::variant<typename Boundaries::Inherit,
											typename Boundaries::Linear,
											typename Boundaries::Quadratic,
											typename Boundaries::Cubic,
											typename Boundaries::Polynomial,
											typename Boundaries::Clamped,
											typename Boundaries::Periodic>;

		static Container<Number> compute_coeffs(
				const Container<Number>& X,
				const Container<Number>& Y,
				Type type = typename Types::Default{},
				BoundariesType boundaries_type = typename Boundaries::Default{}
		) {
			if (X.size() != Y.size()) {
				std::cerr << "Error in function " << __FUNCTION__
						  << ": Vectors X (of size " << X.size()
						  << ") and Y (of size " << Y.size()
						  << ") are not the same size. Returning an empty array..." << std::endl;
				return {};
			}

			const auto n = X.size();

			if (n <= 1) {
				return util::spline::simple_spline<Number,Container>(X, Y, 3);
			}

			Container<Number> delta(n - 1);

			for (SizeT i = 0; i < n - 1; ++i) {
				delta[i] = (Y[i+1] - Y[i]) / (X[i+1] - X[i]);
			}

			Container<Number> derivatives(n);

			if (std::holds_alternative<typename Types::Classic>(type)) {
				// Classic
				for (SizeT i = 1; i < n - 1; ++i) {
					derivatives[i] = (delta[i-1] + delta[i]) / 2;
				}
				if (std::holds_alternative<typename Boundaries::Inherit>(boundaries_type)) {
					derivatives[0] = delta[0];
					derivatives[n-1] = delta[n-2];
				} else if (std::holds_alternative<typename Boundaries::Periodic>(boundaries_type)) {
					derivatives[0] = derivatives[n-1] = (delta[n-2] + delta[0]) / 2;
				}
			} else if (const auto* cardinal = std::get_if<typename Types::Cardinal>(&type)) {
				// Cardinal
				for (SizeT i = 1; i < n - 1; ++i) {
					derivatives[i] = (1 - cardinal->c) * (Y[i+1] - Y[i-1]) / (X[i+1] - X[i-1]);
				}
				if (std::holds_alternative<typename Boundaries::Inherit>(boundaries_type)) {
					derivatives[0] = (1 - cardinal->c) * delta[0];
					derivatives[n-1] = (1 - cardinal->c) * delta[n-2];
				} else if (std::holds_alternative<typename Boundaries::Periodic>(boundaries_type)) {
					derivatives[0] = derivatives[n-1] = (1 - cardinal->c) * (Y[1] - Y[n-1]) / (X[1] - X[n-1]);
				}
			} else if (std::holds_alternative<typename Types::CatmullRom>(type)) {
				// CatmullRom
				for (SizeT i = 1; i < n - 1; ++i) {
					derivatives[i] = (Y[i+1] - Y[i-1]) / (X[i+1] - X[i-1]);
				}
				if (std::holds_alternative<typename Boundaries::Inherit>(boundaries_type)) {
					derivatives[0] = delta[0];
					derivatives[n-1] = delta[n-2];
				} else if (std::holds_alternative<typename Boundaries::Periodic>(boundaries_type)) {
					derivatives[0] = derivatives[n-1] = (Y[1] - Y[n-1]) / (X[1] - X[n-1]);
				}
			} else if (std::holds_alternative<typename Types::Akima>(type) || std::holds_alternative<typename Types::ModifiedAkima>(type)) {
				// Akima / ModifiedAkima
				const auto akima = [&](Number d_im2, Number d_im1, Number d_i, Number d_ip1) {
					auto w1 = std::abs(d_ip1 - d_i);
					auto w2 = std::abs(d_im1 - d_im2);
					if (std::holds_alternative<typename Types::ModifiedAkima>(type)) {
						w1 += std::abs(d_ip1 + d_i) / 2;
						w2 += std::abs(d_im1 + d_im2) / 2;
					}
					return w1 + w2 == 0 ? (d_im1 + d_i) / 2 : (w1 * d_im1 + w2 * d_i) / (w1 + w2);
				};

				for (SizeT i = 2; i < n - 2; ++i) {
					derivatives[i] = akima(delta[i-2], delta[i-1], delta[i], delta[i+1]);
				}

				if (std::holds_alternative<typename Boundaries::Inherit>(boundaries_type)) {
					if (n == 2) {
						derivatives[0] = derivatives[1] = delta[0];
					} else {
						const auto delta_n_2 = 3 * delta[0] - 2 * delta[1];
						const auto delta_n_1 = 2 * delta[0] - delta[1];
						const auto delta_0 = 2 * delta[n-2] - delta[n-3];
						const auto delta_1 = 3 * delta[n-2] - 2 * delta[n-3];
						derivatives[0] = akima(delta_n_2, delta_n_1, delta[0], delta[1]);
						if (n == 3) {
							derivatives[1] = akima(delta_n_1, delta[0], delta[1], delta_0);
						} else {
							derivatives[1] = akima(delta_n_1, delta[0], delta[1], delta[2]);
							derivatives[n-2] = akima(delta[n-4], delta[n-3], delta[n-2], delta_0);
						}
						derivatives[n-1] = akima(delta[n-3], delta[n-2], delta_0, delta_1);
					}
				} else if (std::holds_alternative<typename Boundaries::Periodic>(boundaries_type)) {
					if (n == 2) {
						derivatives[0] = derivatives[1] = delta[0];
					} else if (n == 3) {
						derivatives[0] = derivatives[1] = derivatives[2] = (delta[0] + delta[1]) / 2;
					} else {
						derivatives[0] = akima(delta[n-3], delta[n-2], delta[0], delta[1]);
						derivatives[1] = akima(delta[n-2], delta[0], delta[1], delta[2]);
						derivatives[n-2] = akima(delta[n-4], delta[n-3], delta[n-2], delta[0]);
						derivatives[n-1] = derivatives[0];
					}
				}
			} else if (std::holds_alternative<typename Types::Steffen>(type)) {
				// Steffen
				for (SizeT i = 1; i < n - 1; ++i) {
					const auto s1 = (delta[i-1] > 0) - (delta[i-1] < 0), s2 = (delta[i] > 0) - (delta[i] < 0);
					const auto hi_1 = X[i] - X[i-1], hi = X[i+1] - X[i];
					derivatives[i] = (s1 + s2) * std::min({std::abs(delta[i-1]), std::abs(delta[i]), std::abs((delta[i-1] * hi_1 + delta[i] * hi) / (hi_1 + hi) / 2)});
				}
				if (std::holds_alternative<typename Boundaries::Inherit>(boundaries_type)) {
					if (n == 2) {
						derivatives[0] = derivatives[1] = delta[0];
					} else {
						const auto p1 = delta[0] * (1 + (X[1] - X[0]) / (X[2] - X[0])) - delta[1] * ((X[1] - X[0]) / (X[2] - X[0]));
						const auto s1 = (delta[0] > 0) - (delta[0] < 0), s2 = (p1 > 0) - (p1 < 0);
						derivatives[0] = (s1 + s2) * std::min(std::abs(delta[0]), std::abs(p1 / 2));
						const auto pn = delta[n-2] * (1 + (X[n-1] - X[n-2]) / (X[n-1] - X[n-3])) - delta[n-3] * ((X[n-1] - X[n-2]) / (X[n-1] - X[n-3]));
						const auto s3 = (delta[n-2] > 0) - (delta[n-2] < 0), s4 = (pn > 0) - (pn < 0);
						derivatives[n-1] = (s3 + s4) * std::min(std::abs(delta[n-2]), std::abs(pn / 2));
					}
				} else if (std::holds_alternative<typename Boundaries::Periodic>(boundaries_type)) {
					const auto sn_2 = (delta[n-2] > 0) - (delta[n-2] < 0), s0 = (delta[0] > 0) - (delta[0] < 0);
					const auto hn_2 = X[n-1] - X[n-2], h0 = X[1] - X[0];
					derivatives[0] = derivatives[n-1] = (sn_2 + s0) * std::min({std::abs(delta[n-2]), std::abs(delta[0]), std::abs((delta[n-2] * hn_2 + delta[0] * h0) / (hn_2 + h0) / 2)});
				}
			} else if (std::holds_alternative<typename Types::Zero>(type)) {
				// Zero
				std::fill(derivatives.begin(), derivatives.end(), 0);
			} else if (const auto* expl = std::get_if<typename Types::Explicit>(&type)) {
				// Explicit
				if (expl->derivatives.size() != n) {
					std::cerr << "Error in function " << __FUNCTION__
							<< ": Explicit derivatives (of size " << expl->derivatives.size()
							<< ") and points (of size " << n
							<< ") are not the same size. Returning an empty array..."
							<< std::endl;
					return {};
				}
				derivatives = expl->derivatives;
			}

			const auto polynomial_endpoint_derivative =
				[&](bool left, SizeT degree) -> Number {
					if (degree == 0) {
						return 0;
					}
					const auto count = std::min(degree + 1, n);
					if (count < 2) {
						return 0;
					}
					const auto first = left ? 0 : n - count;
					const auto last = first + count - 1;
					const auto r = left ? first : last;
					Number result = 0;
					for (auto j = first; j <= last; ++j) {
						if (j == r) {
							continue;
						}
						auto basis_derivative = 1 / (X[j] - X[r]);
						for (auto k = first; k <= last; ++k) {
							if (k == j || k == r) {
								continue;
							}
							basis_derivative *= (X[r] - X[k]) / (X[j] - X[k]);
						}
						result += (Y[j] - Y[r]) * basis_derivative;
					}
					return result;
				};

			std::visit(util::misc::overload{
				[&](const typename Boundaries::Inherit&) {
					// Do nothing. This case is covered above
				},
				[&](const typename Boundaries::Linear&) {
					derivatives[0] = delta[0];
					derivatives[n-1] = delta[n-2];
				},
				[&](const typename Boundaries::Quadratic&) {
					derivatives[0] = polynomial_endpoint_derivative(true, 2);
					derivatives[n-1] = polynomial_endpoint_derivative(false, 2);
				},
				[&](const typename Boundaries::Cubic&) {
					derivatives[0] = polynomial_endpoint_derivative(true, 3);
					derivatives[n-1] = polynomial_endpoint_derivative(false, 3);
				},
				[&](const typename Boundaries::Polynomial& polynomial) {
					derivatives[0] = polynomial_endpoint_derivative(true, polynomial.degree);
					derivatives[n-1] = polynomial_endpoint_derivative(false, polynomial.degree);
				},
				[&](const typename Boundaries::Clamped& clamped) {
					derivatives[0] = clamped.left;
					derivatives[n-1] = clamped.right;
				},
				[&](const typename Boundaries::Periodic&) {
					// Do nothing. This case is covered above
				},
			}, boundaries_type);

			Container<Number> coeffs;
			coeffs.resize(4 * (n - 1));

			for (SizeT i = 0; i < n - 1; ++i) {
				const auto h = X[i+1] - X[i];
				const auto h2 = h * h;

				const auto m0 = derivatives[i];
				const auto m1 = derivatives[i+1];

				const auto y0 = Y[i];
				const auto y1 = Y[i+1];

				const auto s = (y1 - y0) / h;

				coeffs[4*i+0] = (m0 + m1 - 2 * s) / h2;
				coeffs[4*i+1] = (3 * s - 2 * m0 - m1) / h;
				coeffs[4*i+2] = m0;
				coeffs[4*i+3] = y0;
			}

			return coeffs;
		}

		HermiteSpline() = default;

		template <typename ContainerXType>
		HermiteSpline(
				ContainerXType&& X,
				const Container<Number>& Y,
				Type type = typename Types::Default{},
				BoundariesType boundaries_type = typename Boundaries::Default{}
		) {
			construct(std::forward<ContainerXType>(X), Y, type, boundaries_type);
		}

		HermiteSpline(const HermiteSpline& other) = default;
		HermiteSpline(HermiteSpline&& other) noexcept = default;

		HermiteSpline& operator=(const HermiteSpline& other) = default;
		HermiteSpline& operator=(HermiteSpline&& other) noexcept = default;

		template <typename ContainerXType>
		void construct(
				ContainerXType&& X,
				const Container<Number>& Y,
				Type type = typename Types::Default{},
				BoundariesType boundaries_type = typename Boundaries::Default{}
		) {
			if (X.size() != Y.size()) {
				std::cerr << "Error in function " << __FUNCTION__
						  << ": Vectors X (of size " << X.size()
						  << ") and Y (of size " << Y.size()
						  << ") are not the same size. Doing nothing..." << std::endl;
				return;
			}
			auto coeffs = compute_coeffs(X, Y, type, boundaries_type);
			if (coeffs.empty() && !X.empty()) {
				std::cerr << "Error in function " << __FUNCTION__
						  << ": failed to construct coefficients. Not changing the object..." << std::endl;
				return;
			}
			_type = type;
			_boundaries_type = boundaries_type;
			_X = std::forward<ContainerXType>(X);
			_coeffs = std::move(coeffs);
		}

		Number eval(const Number& x) const {
			return eval(x, std::distance(_X.begin(), util::misc::first_leq_or_begin(_X.begin(), _X.end(), x)));
		}
		Number eval(const Number& x, SizeT segment) const {
			if (_coeffs.empty()) {
				return NAN;
			} else if (_coeffs.size() == 1) {
				return _coeffs[0];
			}
			segment = std::clamp(segment, static_cast<SizeT>(0), static_cast<SizeT>(_X.size() - 2));
			const auto x_seg = x - _X[segment];
			return ((_coeffs[4*segment] * x_seg + _coeffs[4*segment+1]) * x_seg + _coeffs[4*segment+2]) * x_seg + _coeffs[4*segment+3];
		}

		Container<Number> eval(const Container<Number>& xx, bool sorted = true) const {
			Container<Number> result(xx.size());
			if (sorted) {
				for (SizeT i = 0, i_x = 0; i < xx.size(); ++i) {
					const Number& x = xx[i];
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

		Number operator()(const Number& x) const {
			return eval(x);
		}
		Container<Number> operator()(const Container<Number>& xx) const {
			return eval(xx);
		}

		Type type() const {
			return _type;
		}
		Type boundaries_type() const {
			return _boundaries_type;
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
			return {4*index, 4*(index+1)};
		}

	private:
		Type _type = typename Types::Default{};
		BoundariesType _boundaries_type = typename Boundaries::Default{};
		Container<Number> _X = {};
		Container<Number> _coeffs = {};
	};
}