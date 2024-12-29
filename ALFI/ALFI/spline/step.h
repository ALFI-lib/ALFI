#pragma once

#include <iostream>
#include <cmath>

#include "../config.h"
#include "../util/misc.h"

namespace alfi::spline {
	template <typename Number = DefaultNumber, template <typename> class Container = DefaultContainer>
	class StepSpline {
	public:
		enum class Type {
			LEFT,
			RIGHT,
			MIDDLE,

			Default = LEFT,
		};

		StepSpline() = default;

		template <typename ContainerXType, typename ContainerYType>
		StepSpline(ContainerXType&& X, ContainerYType&& Y, Type type = Type::Default) {
			construct(std::forward<ContainerXType>(X), std::forward<ContainerYType>(Y), type);
		}

		StepSpline(const StepSpline& other) = default;
		StepSpline(StepSpline&& other) noexcept = default;

		StepSpline& operator=(const StepSpline& other) = default;
		StepSpline& operator=(StepSpline&& other) noexcept = default;

		template <typename ContainerXType, typename ContainerYType>
		void construct(ContainerXType&& X, ContainerYType&& Y, Type type = Type::Default) {
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

		Number eval(Number x) const {
			return eval(x, std::distance(_X.begin(), util::misc::first_leq_or_begin(_X.begin(), _X.end(), x)));
		}
		Number eval(Number x, SizeT segment_index) const {
			if (_Y.empty()) {
				return NAN;
			} else if (_Y.size() == 1) {
				return _Y[0];
			}
			segment_index = std::max(std::min(segment_index, static_cast<SizeT>(_X.size() - 2)), static_cast<SizeT>(0));
			if (x <= _X[segment_index]) {
				return _Y[segment_index];
			}
			if (x >= _X[segment_index + 1]) {
				return _Y[segment_index + 1];
			}
			switch (_type) {
				case Type::LEFT: return _Y[segment_index];
				case Type::RIGHT: return _Y[segment_index + 1];
				case Type::MIDDLE: return (_Y[segment_index] + _Y[segment_index + 1]) / 2;
				default: return NAN;
			}
		}

		Container<Number> eval(const Container<Number>& xx, bool sorted = true) const {
			Container result(xx.size());
			if (sorted) {
				for (SizeT i = 0, i_x = 0; i < xx.size(); ++i) {
					const Number x = xx[i];
					while (i_x < _X.size() && x >= _X[i_x+1])
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
		Type _type = Type::Default;
		Container<Number> _X = {};
		Container<Number> _Y = {};
	};
}