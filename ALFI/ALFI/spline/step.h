#pragma once

#include <iostream>
#include <cmath>

#include "../config.h"
#include "../util.h"

namespace alfi::spline {
	template <typename Number = DefaultNumber, typename Container = DefaultContainer<Number>>
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
			return eval(x, std::distance(_X.begin(), util::first_leq_or_begin(_X.begin(), _X.end(), x)));
		}
		Number eval(Number x, size_t segment_index) const {
			if (_X.empty())
				return NAN;
			segment_index = std::max(std::min(segment_index, _X.size() - 2), static_cast<size_t>(0));
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

		Type type() const {
			return _type;
		}

		const Container& X() const & {
			return _X;
		}
		Container&& X() && {
			return std::move(_X);
		}

		const Container& Y() const & {
			return _Y;
		}
		Container&& Y() && {
			return std::move(_Y);
		}

	private:
		Type _type = Type::Default;
		Container _X = {};
		Container _Y = {};
	};
}