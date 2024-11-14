#pragma once

#include <cmath>
#include <functional>
#include <memory>

namespace alfi::dist {
	enum class Type {
		GENERAL,
		UNIFORM,
		CHEBYSHEV,
		CHEBYSHEV_STRETCHED,
		CHEBYSHEV_ELLIPSE,
		CHEBYSHEV_ELLIPSE_STRETCHED,
		CIRCLE_PROJECTION,
		ELLIPSE_PROJECTION,
		SIGMOID,
		SIGMOID_STRETCHED,
		ERF,
		ERF_STRETCHED,
	};

	template <typename Number = double>
	class Base {
	public:
		virtual ~Base() = default;

		Type type() const { return _impl->type; }
		size_t n() const { return _impl->n; }
		const Number& a() const { return _impl->a; }
		const Number& b() const { return _impl->b; }

		Number operator[](size_t index) const {
			return _impl->get(index);
		}
		Number get(size_t index) const {
			return _impl->get(index);
		}

		Number first() const {
			return _impl->get(0);
		}
		Number last() const {
			return _impl->get(_impl->n - 1);
		}

	protected:
		class BaseImpl {
		public:
			BaseImpl(Type type, size_t n, Number a, Number b) : type(type), n(n), a(std::move(a)), b(std::move(b)) {}
			virtual ~BaseImpl() = default;
			virtual Number get(size_t index) const = 0;
			const Type type;
			const size_t n;
			const Number a, b;
		};

	public:
		class Iterator final : public std::iterator<std::random_access_iterator_tag, Number> {
		public:
			Iterator(std::shared_ptr<BaseImpl> impl, size_t index)
				: _impl(impl), _index(index) {}

			Iterator& operator++() {
				++_index;
				return *this;
			}
			Iterator operator++(int) {
				Iterator copy = *this;
				++_index;
				return copy;
			}

			Iterator& operator--() {
				--_index;
				return *this;
			}
			Iterator operator--(int) {
				Iterator copy = *this;
				--_index;
				return copy;
			}

			Iterator& operator+=(ptrdiff_t offset) {
				_index += offset;
				return *this;
			}
			Iterator& operator-=(ptrdiff_t offset) {
				_index -= offset;
				return *this;
			}

			bool operator==(const Iterator& other) const {
				return _index == other._index && _impl == other._impl;
			}
			bool operator!=(const Iterator& other) const {
				return !(*this == other);
			}

			Number operator*() const {
				return _impl->get(_index);
			}
			Number operator[](ptrdiff_t offset) const {
				return _impl->get(_index + offset);
			}

		private:
			std::shared_ptr<BaseImpl> _impl;
			size_t _index;
		};

		using ReverseIterator = std::reverse_iterator<Iterator>;

		Iterator begin() const {
			return Iterator(_impl, 0);
		}
		Iterator end() const {
			return Iterator(_impl, _impl->n);
		}

		ReverseIterator rbegin() const {
			return ReverseIterator(end());
		}
		ReverseIterator rend() const {
			return ReverseIterator(begin());
		}

		template <typename Container>
		Container collect() const {
			Container result;
			result.reserve(n());
			for (auto it = begin(); it != end(); ++it) {
				result.push_back(*it);
			}
			return result;
		}

	protected:
		explicit Base(std::shared_ptr<BaseImpl> impl) : _impl(impl) {}

		std::shared_ptr<BaseImpl> _impl;
	};

	template <typename Number = double>
	class General final : public Base<Number> {
	public:
		General(size_t n, Number a, Number b, std::function<Number(size_t)> lambda)
			: Base<Number>(std::make_shared<GeneralImpl>(n, std::move(a), std::move(b), lambda)) {}

	private:
		class GeneralImpl final : public Base<Number>::BaseImpl {
		public:
			GeneralImpl(size_t n, Number a, Number b, std::function<Number(size_t)> lambda)
				: Base<Number>::BaseImpl(Type::GENERAL, n, std::move(a), std::move(b)), _lambda(lambda) {}

			virtual Number get(size_t index) const override {
				return _lambda(index);
			}

		private:
			const std::function<Number(size_t)> _lambda;
		};
	};

	template <typename Number = double>
	class Uniform final : public Base<Number> {
	public:
		Uniform(size_t n, Number a, Number b) : Base<Number>(std::make_shared<UniformImpl>(n, std::move(a), std::move(b))) {}
	private:

		class UniformImpl final : public Base<Number>::BaseImpl {
		public:
			UniformImpl(size_t n, Number a, Number b) : Base<Number>::BaseImpl(Type::UNIFORM, n, std::move(a), std::move(b)) {}

			virtual Number get(size_t index) const override {
				if (this->n == 0)
					return static_cast<Number>(NAN);
				if (index == 0 && this->n == 1)
					return (this->a + this->b) / 2;
				return this->a + (this->b - this->a) * static_cast<Number>(index) / static_cast<Number>(this->n - 1);
			}
		};
	};

	template <typename Number = double>
	class Chebyshev final : public Base<Number> {
	public:
		Chebyshev(size_t n, Number a, Number b) : Base<Number>(std::make_shared<ChebyshevImpl>(n, std::move(a), std::move(b))) {}

	private:
		class ChebyshevImpl final : public Base<Number>::BaseImpl {
		public:
			ChebyshevImpl(size_t n, Number a, Number b) : Base<Number>::BaseImpl(Type::CHEBYSHEV, n, std::move(a), std::move(b)) {}

			virtual Number get(size_t index) const override {
				if (this->n == 0)
					return static_cast<Number>(NAN);
				if (index * 2 == this->n - 1)
					return (this->a + this->b) / 2;
				Number x = std::cos(M_PI * (2 * (static_cast<Number>(this->n) - static_cast<Number>(index)) - 1) / (2 * static_cast<Number>(this->n)));
				x = (x + 1) * (this->b - this->a) / 2 + this->a;
				return x;
			}
		};
	};

	template <typename Number = double>
	class ChebyshevStretched final : public Base<Number> {
	public:
		ChebyshevStretched(size_t n, Number a, Number b) 
			: Base<Number>(std::make_shared<ChebyshevStretchedImpl>(n, std::move(a), std::move(b))) {}
	
		const Number& stretched_a() const { return static_cast<ChebyshevStretchedImpl&>(*this->_impl).stretched_a; }
		const Number& stretched_b() const { return static_cast<ChebyshevStretchedImpl&>(*this->_impl).stretched_b; }
	
	private:
		class ChebyshevStretchedImpl final : public Base<Number>::BaseImpl {
		public:
			ChebyshevStretchedImpl(size_t n, Number a, Number b) 
				: Base<Number>::BaseImpl(Type::CHEBYSHEV_STRETCHED, n, std::move(a), std::move(b)) {
				const Number mid = this->a + (this->b - this->a) / 2;
				if (this->n > 1) {
					const Number min = this->a + (this->b - this->a) / 2 * (std::cos(M_PI * (2 * static_cast<Number>(this->n) - 1) / (2 * static_cast<Number>(this->n))) + 1);
					const Number max = this->a + (this->b - this->a) / 2 * (std::cos(M_PI / (2 * static_cast<Number>(this->n))) + 1);
					stretched_a = mid - (mid - this->a) * (mid - this->a) / (mid - min);
					stretched_b = mid + (this->b - mid) * (this->b - mid) / (max - mid);
				} else {
					stretched_a = this->a;
					stretched_b = this->b;
				}
			}
	
			virtual Number get(size_t index) const override {
				if (this->n == 0)
					return static_cast<Number>(NAN);
				if (index * 2 == this->n - 1)
					return (this->a + this->b) / 2;
				if (index == 0)
					return this->a;
				if (index == this->n - 1)
					return this->b;
				Number x = std::cos(M_PI * (2 * (static_cast<Number>(this->n) - static_cast<Number>(index)) - 1) / (2 * static_cast<Number>(this->n)));
				x = (x + 1) * (stretched_b - stretched_a) / 2 + stretched_a;
				return x;
			}
	
			Number stretched_a;
			Number stretched_b;
		};
	};

	template <typename Number = double>
	class ChebyshevEllipse final : public Base<Number> {
	public:
		ChebyshevEllipse(size_t n, Number a, Number b, Number B)
			: Base<Number>(std::make_shared<ChebyshevEllipseImpl>(n, std::move(a), std::move(b), std::move(B))) {}

		const Number& ratio() const { return static_cast<ChebyshevEllipseImpl&>(*this->_impl).ratio; }

	private:
		class ChebyshevEllipseImpl final : public Base<Number>::BaseImpl {
		public:
			ChebyshevEllipseImpl(size_t n, Number a, Number b, Number B)
				: Base<Number>::BaseImpl(Type::CHEBYSHEV_ELLIPSE, n, std::move(a), std::move(b)), ratio(std::move(B)) {}

			virtual Number get(size_t index) const override {
				if (this->n == 0)
					return static_cast<Number>(NAN);
				if (index * 2 == this->n - 1)
					return (this->a + this->b) / 2;
				const Number theta = M_PI * (2 * static_cast<Number>(index) + 1) / (2 * static_cast<Number>(this->n));
				Number x = 1 / std::sqrt(1 + std::pow(std::tan(theta) / ratio, 2));
				if (index < this->n / 2) {
					x = 1 - x;
				} else {
					x = 1 + x;
				}
				return this->a + (this->b - this->a) / 2 * x;
			}

			const Number ratio;
		};
	};

	template <typename Number = double>
	class ChebyshevEllipseStretched final : public Base<Number> {
	public:
		ChebyshevEllipseStretched(size_t n, Number a, Number b, Number B)
			: Base<Number>(std::make_shared<ChebyshevEllipseStretchedImpl>(n, std::move(a), std::move(b), std::move(B))) {}

		const Number& ratio() const { return static_cast<ChebyshevEllipseStretchedImpl&>(*this->_impl).ratio; }
		const Number& stretched_a() const { return static_cast<ChebyshevEllipseStretchedImpl&>(*this->_impl).stretched_a; }
		const Number& stretched_b() const { return static_cast<ChebyshevEllipseStretchedImpl&>(*this->_impl).stretched_b; }

	private:
		class ChebyshevEllipseStretchedImpl final : public Base<Number>::BaseImpl {
		public:
			ChebyshevEllipseStretchedImpl(size_t n, Number a, Number b, Number B)
				: Base<Number>::BaseImpl(Type::CHEBYSHEV_ELLIPSE_STRETCHED, n, std::move(a), std::move(b)), ratio(std::move(B)) {
				const Number mid = this->a + (this->b - this->a) / 2;
				if (this->n > 1) {
					const Number theta = M_PI / (2 * static_cast<Number>(this->n));
					const Number x = 1 / std::sqrt(1 + std::pow(std::tan(theta) / B, 2));
					const Number min = this->a + (this->b - this->a) / 2 * (1 - x);
					const Number max = this->a + (this->b - this->a) / 2 * (1 + x);
					stretched_a = mid - (mid - this->a) * (mid - this->a) / (mid - min);
					stretched_b = mid + (this->b - mid) * (this->b - mid) / (max - mid);
				} else {
					stretched_a = this->a;
					stretched_b = this->b;
				}
			}

			virtual Number get(size_t index) const override {
				if (this->n == 0)
					return static_cast<Number>(NAN);
				if (index * 2 == this->n - 1)
					return (this->a + this->b) / 2;
				if (index == 0)
					return this->a;
				if (index == this->n - 1)
					return this->b;
				const Number theta = M_PI * (2 * static_cast<Number>(index) + 1) / (2 * static_cast<Number>(this->n));
				Number x = 1 / std::sqrt(1 + std::pow(std::tan(theta) / ratio, 2));
				if (index < this->n / 2) {
					x = 1 - x;
				} else {
					x = 1 + x;
				}
				return stretched_a + (stretched_b - stretched_a) / 2 * x;
			}

			const Number ratio;
			Number stretched_a;
			Number stretched_b;
		};
	};

	template <typename Number = double>
	class CircleProj final : public Base<Number> {
	public:
		CircleProj(size_t n, Number a, Number b)
			: Base<Number>(std::make_shared<CircleProjImpl>(n, std::move(a), std::move(b))) {}

	private:
		class CircleProjImpl final : public Base<Number>::BaseImpl {
		public:
			CircleProjImpl(size_t n, Number a, Number b)
				: Base<Number>::BaseImpl(Type::CIRCLE_PROJECTION, n, std::move(a), std::move(b)) {}

			virtual Number get(size_t index) const override {
				if (this->n == 0)
					return static_cast<Number>(NAN);
				if (index == 0 && this->n == 1)
					return (this->a + this->b) / 2;
				if (index == 0)
					return this->a;
				if (index == this->n - 1)
					return this->b;
				const Number angle = M_PI * static_cast<Number>(index) / (static_cast<Number>(this->n) - 1);
				const Number x = 1 - std::cos(angle);
				return this->a + (this->b - this->a) * x / 2;
				// Alternative method using sinus.
				// Since 1 + sin(pi * (x - 0.5)) = 1 - cos(pi * x).
				// const Number x = static_cast<Number>(index) / (static_cast<Number>(this->n) - 1);
				// return this->a + (this->b - this->a) / 2 * (1 + sin(M_PI * (x - 0.5)));
			}
		};
	};

	template <typename Number = double>
	class EllipseProj final : public Base<Number> {
	public:
		EllipseProj(size_t n, Number a, Number b, Number B)
			: Base<Number>(std::make_shared<EllipseProjImpl>(n, std::move(a), std::move(b), std::move(B))) {}
	
		const Number& ratio() const { return static_cast<EllipseProjImpl&>(*this->_impl).ratio; }

	private:
		class EllipseProjImpl final : public Base<Number>::BaseImpl {
		public:
			EllipseProjImpl(size_t n, Number a, Number b, Number B)
				: Base<Number>::BaseImpl(Type::ELLIPSE_PROJECTION, n, std::move(a), std::move(b)), ratio(std::move(B)) {}
			
			virtual Number get(size_t index) const override {
				if (this->n == 0)
					return static_cast<Number>(NAN);
				if (index * 2 == this->n - 1)
					return (this->a + this->b) / 2;
				if (index == 0)
					return this->a;
				if (index == this->n - 1)
					return this->b;
				const Number theta = M_PI * static_cast<Number>(index) / (static_cast<Number>(this->n) - 1);
				Number x = 1 / std::sqrt(1 + std::pow(std::tan(theta) / ratio, 2));
				if (index < this->n / 2) {
					x = 1 - x;
				} else {
					x = 1 + x;
				}
				return this->a + (this->b - this->a) / 2 * x;
			}

			const Number ratio;
		};
	};

	template <typename Number = double>
	class Sigmoid final : public Base<Number> {
	public:
		Sigmoid(size_t n, Number a, Number b, Number steepness)
			: Base<Number>(std::make_shared<SigmoidImpl>(n, std::move(a), std::move(b), std::move(steepness))) {}

		const Number& steepness() const { return static_cast<SigmoidImpl&>(*this->_impl).steepness; }

	private:
		class SigmoidImpl final : public Base<Number>::BaseImpl {
		public:
			SigmoidImpl(size_t n, Number a, Number b, Number steepness)
				: Base<Number>::BaseImpl(Type::SIGMOID, n, std::move(a), std::move(b)), steepness(std::move(steepness)) {}

			virtual Number get(size_t index) const override {
				if (this->n == 0)
					return static_cast<Number>(NAN);
				if (index * 2 == this->n - 1)
					return (this->a + this->b) / 2;
				const Number x = static_cast<Number>(index) / (static_cast<Number>(this->n) - 1);
				const Number sigmoidValue = 1.0 / (1.0 + exp(-steepness * (x - 0.5)));
				return this->a + (this->b - this->a) * sigmoidValue;
			}

			const Number steepness;
		};
	};

	template <typename Number = double>
	class SigmoidStretched final : public Base<Number> {
	public:
		SigmoidStretched(size_t n, Number a, Number b, Number steepness)
			: Base<Number>(std::make_shared<SigmoidStretchedImpl>(n, std::move(a), std::move(b), std::move(steepness))) {}

		const Number& steepness() const { return static_cast<SigmoidStretchedImpl&>(*this->_impl).steepness; }
		const Number& stretched_a() const { return static_cast<SigmoidStretchedImpl&>(*this->_impl).stretched_a; }
		const Number& stretched_b() const { return static_cast<SigmoidStretchedImpl&>(*this->_impl).stretched_b; }
		
	private:
		class SigmoidStretchedImpl final : public Base<Number>::BaseImpl {
		public:
			SigmoidStretchedImpl(size_t n, Number a, Number b, Number steepness)
				: Base<Number>::BaseImpl(Type::SIGMOID_STRETCHED, n, std::move(a), std::move(b)), steepness(std::move(steepness)) {
				const Number mid = this->a + (this->b - this->a) / 2;
				if (this->n > 1) {
					const Number min = this->a + (this->b - this->a) / (1 + exp(steepness / 2));
					const Number max = this->a + (this->b - this->a) / (1 + exp(-steepness / 2));
					stretched_a = mid - (mid - this->a) * (mid - this->a) / (mid - min);
					stretched_b = mid + (this->b - mid) * (this->b - mid) / (max - mid);
				} else {
					stretched_a = this->a;
					stretched_b = this->b;
				}
			}

			virtual Number get(size_t index) const override {
				if (this->n == 0)
					return static_cast<Number>(NAN);
				if (index * 2 == this->n - 1)
					return (this->a + this->b) / 2;
				if (index == 0)
					return this->a;
				if (index == this->n - 1)
					return this->b;
				const Number x = static_cast<Number>(index) / (static_cast<Number>(this->n) - 1);
				const Number sigmoidValue = 1.0 / (1.0 + exp(-steepness * (x - 0.5)));
				return stretched_a + (stretched_b - stretched_a) * sigmoidValue;
			}

			const Number steepness;
			Number stretched_a;
			Number stretched_b;
		};
	};

	template <typename Number = double>
	class Erf final : public Base<Number> {
	public:
		Erf(size_t n, Number a, Number b, Number steepness)
			: Base<Number>(std::make_shared<ErfImpl>(n, std::move(a), std::move(b), std::move(steepness))) {}

		const Number& steepness() const { return static_cast<ErfImpl&>(*this->_impl).steepness; }

	private:
		class ErfImpl final : public Base<Number>::BaseImpl {
		public:
			ErfImpl(size_t n, Number a, Number b, Number steepness)
				: Base<Number>::BaseImpl(Type::ERF, n, std::move(a), std::move(b)), steepness(std::move(steepness)) {}

			virtual Number get(size_t index) const override {
				if (this->n == 0)
					return static_cast<Number>(NAN);
				if (index * 2 == this->n - 1)
					return (this->a + this->b) / 2;
				const Number x = static_cast<Number>(index) / (static_cast<Number>(this->n) - 1);
				const Number erf_value = erf(steepness * (x - 0.5));
				return this->a + (this->b - this->a) * (1 + erf_value) / 2;
			}

			const Number steepness;
		};
	};

	template <typename Number = double>
	class ErfStretched final : public Base<Number> {
	public:
		ErfStretched(size_t n, Number a, Number b, Number steepness)
			: Base<Number>(std::make_shared<ErfStretchedImpl>(n, std::move(a), std::move(b), std::move(steepness))) {}

		const Number& steepness() const { return static_cast<ErfStretchedImpl&>(*this->_impl).steepness; }
		const Number& stretched_a() const { return static_cast<ErfStretchedImpl&>(*this->_impl).stretched_a; }
		const Number& stretched_b() const { return static_cast<ErfStretchedImpl&>(*this->_impl).stretched_b; }

	private:
		class ErfStretchedImpl final : public Base<Number>::BaseImpl {
		public:
			ErfStretchedImpl(size_t n, Number a, Number b, Number steepness)
				: Base<Number>::BaseImpl(Type::ERF_STRETCHED, n, std::move(a), std::move(b)), steepness(std::move(steepness)) {
				const Number mid = this->a + (this->b - this->a) / 2;
				if (this->n > 1) {
					const Number min = this->a + (this->b - this->a) * (1 + erf(-steepness / 2)) / 2;
					const Number max = this->a + (this->b - this->a) * (1 + erf(steepness / 2)) / 2;
					stretched_a = mid - (mid - this->a) * (mid - this->a) / (mid - min);
					stretched_b = mid + (this->b - mid) * (this->b - mid) / (max - mid);
				} else {
					stretched_a = this->a;
					stretched_b = this->b;
				}
			}

			virtual Number get(size_t index) const override {
				if (this->n == 0)
					return static_cast<Number>(NAN);
				if (index * 2 == this->n - 1)
					return (this->a + this->b) / 2;
				if (index == 0)
					return this->a;
				if (index == this->n - 1)
					return this->b;
				const Number x = static_cast<Number>(index) / (static_cast<Number>(this->n) - 1);
				const Number erf_value = erf(steepness * (x - 0.5));
				return stretched_a + (stretched_b - stretched_a) * (1 + erf_value) / 2;
			}

			const Number steepness;
			Number stretched_a;
			Number stretched_b;
		};
	};
}