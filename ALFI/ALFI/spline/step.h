#pragma once

#include <algorithm>
#include <iostream>
#include <cmath>
#include <variant>

#include "../config.h"
#include "../util/misc.h"

namespace alfi::spline {
	template <typename Number = DefaultNumber, template <typename, typename...> class Container = DefaultContainer>
	class StepSpline {
	public:
		struct Types final {
			struct Left final {};
			struct Middle final {};
			struct Right final {};
			struct Fraction final {
				explicit Fraction(Number fraction) : f(std::move(fraction)) {}
				Number f;
			};
			struct Proportion final {
				Proportion(Number left, Number right) : l(std::move(left)), r(std::move(right)) {}
				Number l, r;
			};
			using Default = Left;
		};

		using Type = std::variant<typename Types::Left,
								  typename Types::Middle,
								  typename Types::Right,
								  typename Types::Fraction,
								  typename Types::Proportion>;

		StepSpline() = default;

		template <typename ContainerXType, typename ContainerYType>
		StepSpline(ContainerXType&& X, ContainerYType&& Y, Type type = typename Types::Default{}) {
			construct(std::forward<ContainerXType>(X), std::forward<ContainerYType>(Y), type);
		}

		StepSpline(const StepSpline& other) = default;
		StepSpline(StepSpline&& other) noexcept = default;

		StepSpline& operator=(const StepSpline& other) = default;
		StepSpline& operator=(StepSpline&& other) noexcept = default;

		template <typename ContainerXType, typename ContainerYType>
		void construct(ContainerXType&& X, ContainerYType&& Y, Type type = typename Types::Default{}) {
			if (X.size() != Y.size()) {
				std::cerr << "Error in function " << __FUNCTION__
						  << ": Vectors X (of size " << X.size()
						  << ") and Y (of size " << Y.size()
						  << ") are not the same size. Doing nothing..." << std::endl;
				return;
			}
			_type = type;
			_X = std::forward<ContainerXType>(X);
			_Y = std::forward<ContainerYType>(Y);
		}

		Number eval(const Number& x) const {
			return eval(x, std::distance(_X.begin(), util::misc::first_leq_or_begin(_X.begin(), _X.end(), x)));
		}
		Number eval(const Number& x, SizeT segment) const {
			if (_Y.empty()) {
				return NAN;
			} else if (_Y.size() == 1) {
				return _Y[0];
			}
			segment = std::clamp(segment, static_cast<SizeT>(0), static_cast<SizeT>(_X.size() - 2));
			if (x <= _X[segment]) {
				return _Y[segment];
			}
			if (x >= _X[segment+1]) {
				return _Y[segment+1];
			}
			return std::visit(util::misc::overload{
				[&](const typename Types::Left&) { return _Y[segment]; },
				[&](const typename Types::Middle&) { return (_Y[segment] + _Y[segment+1]) / 2; },
				[&](const typename Types::Right&) { return _Y[segment+1]; },
				[&](const typename Types::Fraction& f) { return _Y[segment] + f.f*(_Y[segment+1] - _Y[segment]); },
				[&](const typename Types::Proportion& p) { return (p.r*_Y[segment] + p.l*_Y[segment+1]) / (p.l + p.r); },
			}, _type);
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

		const Container<Number>& X() const & {
			return _X;
		}
		Container<Number>&& X() && {
			return std::move(_X);
		}

		const Container<Number>& Y() const & {
			return _Y;
		}
		Container<Number>&& Y() && {
			return std::move(_Y);
		}

	private:
		Type _type = typename Types::Default{};
		Container<Number> _X = {};
		Container<Number> _Y = {};
	};
}