#pragma once

#include <vector>

namespace alfi::poly {

	template <typename Number = double>
	class Polynomial final {
	public:
		Polynomial() = default;

		Polynomial(const Polynomial& other) = default;

		Polynomial(Polynomial&& other) noexcept
			: _coeffs(std::move(other._coeffs)) {}

		Polynomial(std::initializer_list<Number> coeffs)
			: _coeffs(std::move(coeffs)) {}

		// ReSharper disable once CppNonExplicitConvertingConstructor
		Polynomial(std::vector<Number> coeffs) // NOLINT(*-explicit-constructor)
			: _coeffs(std::move(coeffs)) {}

		template <typename Container,
			std::enable_if_t<
				!std::is_same_v<std::remove_cvref_t<Container>, Polynomial>
				&& !std::is_same_v<std::remove_cvref_t<Container>, std::vector<Number>>, bool> = true>
		// ReSharper disable once CppNonExplicitConvertingConstructor
		Polynomial(const Container& coeffs) // NOLINT(*-explicit-constructor)
			: _coeffs(coeffs.begin(), coeffs.end()) {}

		template <typename Container,
			std::enable_if_t<
				!std::is_same_v<std::remove_cvref_t<Container>, Polynomial>
				&& !std::is_same_v<std::remove_cvref_t<Container>, std::vector<Number>>, bool> = true>
		// ReSharper disable once CppNonExplicitConvertingConstructor
		Polynomial(Container&& coeffs) // NOLINT(*-explicit-constructor)
			: _coeffs(std::make_move_iterator(coeffs.begin()), std::make_move_iterator(coeffs.end())) {}

		template <typename InputIterator>
		Polynomial(InputIterator begin, InputIterator end)
			: _coeffs(begin, end) {}

		Polynomial& operator=(Polynomial other) noexcept {
			_coeffs = std::move(other._coeffs);
			return *this;
		}

		size_t size() const { return _coeffs.size(); }

		size_t degree() const {
			return size() == 0 ? 0 : size() - 1;
		}

		size_t max_degree() const {
			size_t max_deg = degree();
			while (max_deg > 0 && at_degree(max_deg) == 0) {
				--max_deg;
			}
			return max_deg;
		}

		void set_size(size_t s) {
			if (s > size()) {
				pad_left(s - size());
			} else if (s < size()) {
				crop_left(size() - s);
			}
		}
		void set_degree(size_t deg) {
			set_size(deg + 1);
		}

		void reserve_size(size_t s) {
			if (s > size()) {
				pad_left(s - size());
			}
		}
		void reserve_degree(size_t deg) {
			reserve_size(deg + 1);
		}

		[[nodiscard]] bool empty() const { return size() == 0; }

		void clear() { _coeffs.clear(); }

		Number& at_degree(size_t deg) {
			reserve_degree(deg);
			return _coeffs[degree() - deg];
		}
		const Number& at_degree(size_t deg) const {
			if (empty() || deg > degree()) {
				static const Number zero = 0;
				return zero;
			}
			return _coeffs[degree() - deg];
		}

		Number& at_index(size_t index) {
			return _coeffs[index];
		}
		const Number& at_index(size_t index) const {
			return _coeffs[index];
		}

		std::vector<Number>& coeffs() & {
			return _coeffs;
		}
		std::vector<Number> coeffs() && {
			return std::move(_coeffs);
		}
		const std::vector<Number>& coeffs() const & {
			return _coeffs;
		}
		std::vector<Number> coeffs() const && {
			return std::move(_coeffs);
		}
		std::vector<Number> move_coeffs() {
			std::vector<Number> coeffs = std::move(_coeffs);
			_coeffs = {};
			return coeffs;
		}

		bool operator==(const Polynomial& other) const {
			if (size() < other.size()) {
				return other == *this;
			}
			// size() >= other.size()
			const size_t size_diff = size() - other.size();
			for (size_t i = 0; i < size_diff; ++i) {
				if (at_index(i) != 0) {
					return false;
				}
			}
			for (size_t i = size_diff; i < size(); ++i) {
				if (at_index(i) != other.at_index(i - size_diff)) {
					return false;
				}
			}
			return true;
		}
		bool operator!=(const Polynomial& other) const {
			return !(*this == other);
		}

		bool is_exactly_equal(const Polynomial& other) const {
			return _coeffs == other._coeffs;
		}
		bool is_not_exactly_equal(const Polynomial& other) const {
			return _coeffs != other._coeffs;
		}

		Number operator()(Number x) const {
			return eval(x);
		}
		Number eval(const Number& x) const {
			if (empty()) {
				return 0;
			}
			Number result = _coeffs[0];
			for (size_t i = 1; i < size(); ++i) {
				result = result * x + _coeffs[i];
			}
			return result;
		}

		template <template <typename> class Container = std::vector>
		Container<Number> operator()(const Container<Number>& X) const {
			return eval(X);
		}
		template <template <typename> class Container = std::vector>
		Container<Number> eval(const Container<Number> X) const {
			Container<Number> result;
			result.reserve(X.size());
			for (const Number& x : X) {
				result.push_back(eval(x));
			}
			return result;
		}

		Polynomial derivative() const {
			if (empty()) {
				return Polynomial({});
			}
			if (size() == 1) {
				return Polynomial({0});
			}
			Polynomial result (std::vector<Number>(size() - 1));
			for (size_t i = 0; i < result.size(); ++i) {
				result.at_degree(i) = at_degree(i + 1) * (i + 1);
			}
			return result;
		}

		Polynomial derivative(size_t order) const {
			if (order == 0 || empty()) {
				return *this;
			}
			if (size() == 1) {
				return Polynomial({0});
			}
			Polynomial result (std::vector<Number>(size() - order));
			for (size_t i = 0; i < result.size(); ++i) {
				Number cur_coeff = at_degree(i + order);
				for (size_t o = 1; o <= order; ++o) {
					cur_coeff *= i + o;
				}
				result.at_degree(i) = cur_coeff;
			}
			return result;
		}

		Polynomial integral(const Number& C = 0) const {
			if (empty()) {
				return Polynomial({C});
			}
			Polynomial result (std::vector<Number>(size() + 1));
			for (size_t i = 0; i < size(); ++i) {
				result.at_degree(i + 1) = at_degree(i) / (i + 1);
			}
			result.at_degree(0) = C;
			return result;
		}

		Polynomial integral(const Number& C, size_t order) const {
			if (order == 0) {
				return *this;
			}
			Polynomial result (std::vector<Number>(size() + order));
			for (size_t i = 0; i < size(); ++i) {
				Number cur_coeff = at_degree(i);
				for (size_t o = 1; o <= order; ++o) {
					cur_coeff /= i + o;
				}
				result.at_degree(i + order) = cur_coeff;
			}
			result.at_degree(0) = C;
			for (size_t d = 1; d < order; ++d) {
				result.at_degree(d) = result.at_degree(d - 1) / d;
			}
			return result;
		}

		template <template <typename> class Container = std::vector>
		Polynomial integral(const Container<Number>& constants) const {
			if (constants.size() == 0) {
				return *this;
			}
			Polynomial result (std::vector<Number>(size() + constants.size()));
			for (size_t i = 0; i < size(); ++i) {
				Number cur_coeff = at_degree(i);
				for (size_t o = 1; o <= constants.size(); ++o) {
					cur_coeff /= i + o;
				}
				result.at_degree(i + constants.size()) = cur_coeff;
			}
			result.at_degree(0) = constants[0];
			Number divisor = 1;
			for (size_t d = 1; d < constants.size(); ++d) {
				result.at_degree(d) = constants[d] / divisor;
				divisor *= static_cast<Number>(d);
			}
			return result;
		}

		Polynomial copy() const {
			return Polynomial(*this);
		}

		Polynomial& normalize() {
			if (empty()) {
				_coeffs.push_back(0);
				return *this;
			}

			size_t non_zero_index = 0;
			while (non_zero_index < _coeffs.size() && _coeffs[non_zero_index] == 0) {
				++non_zero_index;
			}

			if (non_zero_index == _coeffs.size()) {
				_coeffs.clear();
				_coeffs.push_back(0);
			} else {
				_coeffs.erase(_coeffs.begin(), _coeffs.begin() + non_zero_index);
			}
			return *this;
		}

		Polynomial& pad_left(size_t pad_size, const Number& pad_value = 0) {
			_coeffs.insert(_coeffs.begin(), pad_size, pad_value);
			return *this;
		}
		Polynomial& pad_right(size_t pad_size, const Number& pad_value = 0) {
			_coeffs.insert(_coeffs.end(), pad_size, pad_value);
			return *this;
		}

		Polynomial& crop_left(size_t crop_size) {
			if (crop_size < size()) {
				_coeffs.erase(_coeffs.begin(), _coeffs.begin() + crop_size);
			} else {
				_coeffs.clear();
			}
			return *this;
		}
		Polynomial& crop_right(size_t crop_size) {
			if (crop_size < size()) {
				_coeffs.resize(size() - crop_size);
			} else {
				_coeffs.clear();
			}
			return *this;
		}

		Polynomial& replace_nan(const Number& value = 0) {
			for (Number& c : _coeffs) {
				if (std::isnan(c)) {
					c = value;
				}
			}
			return *this;
		}
		Polynomial& replace_inf(const Number& value = 0) {
			for (Number& c : _coeffs) {
				if (std::isinf(c)) {
					c = value;
				}
			}
			return *this;
		}

		Polynomial operator+() const {
			return Polynomial(*this);
		}
		Polynomial operator-() const {
			Polynomial result(*this);
			for (Number& c : result._coeffs) {
				c = -c;
			}
			return result;
		}

		Polynomial& operator+=(const Number& value) {
			if (empty()) {
				_coeffs.push_back(value);
			} else {
				at_degree(0) += value;
			}
			return *this;
		}
		Polynomial& operator-=(const Number& value) {
			if (empty()) {
				_coeffs.push_back(-value);
			} else {
				at_degree(0) -= value;
			}
			return *this;
		}
		Polynomial& operator*=(const Number& value) {
			for (auto& c : _coeffs)
				c *= value;
			return *this;
		}
		Polynomial& operator/=(const Number& value) {
			for (auto& c : _coeffs)
				c /= value;
			return *this;
		}

		Polynomial operator+(const Number& value) const { return Polynomial(*this) += value; }
		Polynomial operator-(const Number& value) const { return Polynomial(*this) -= value; }
		Polynomial operator*(const Number& value) const { return Polynomial(*this) *= value; }
		Polynomial operator/(const Number& value) const { return Polynomial(*this) /= value; }

		Polynomial& operator+=(const Polynomial& other) {
			if (size() >= other.size()) {
				const size_t size_diff = size() - other.size();
				for (size_t i = size_diff; i < size(); ++i) {
					_coeffs[i] += other._coeffs[i - size_diff];
				}
			} else {
				const size_t size_diff = other.size() - size();
				_coeffs.insert(_coeffs.begin(), other._coeffs.begin(), other._coeffs.begin() + size_diff);
				for (size_t i = size_diff; i < size(); ++i) {
					_coeffs[i] += other._coeffs[i];
				}
			}
			return *this;
		}

		Polynomial& operator-=(const Polynomial& other) {
			if (size() >= other.size()) {
				const size_t size_diff = size() - other.size();
				for (size_t i = size_diff; i < size(); ++i) {
					_coeffs[i] -= other._coeffs[i - size_diff];
				}
			} else {
				const size_t size_diff = other.size() - size();
				_coeffs.insert(_coeffs.begin(), size_diff, 0);
				for (size_t i = 0; i < size(); ++i) {
					_coeffs[i] -= other._coeffs[i];
				}
			}
			return *this;
		}

		Polynomial& operator*=(const Polynomial& other) {
			if (empty()) {
				return *this;
			}
			if (other.empty()) {
				_coeffs = {};
				return *this;
			}
			std::vector<Number> new_coeffs (size() + other.size() - 1, 0);
			for (size_t i = 0; i < size(); ++i) {
				for (size_t j = 0; j < other.size(); ++j) {
					new_coeffs[i + j] += _coeffs[i] * other._coeffs[j];
				}
			}
			_coeffs = std::move(new_coeffs);
			return *this;
		}

		Polynomial operator+(const Polynomial& other) const { return Polynomial(*this) += other; }
		Polynomial operator-(const Polynomial& other) const { return Polynomial(*this) -= other; }
		Polynomial operator*(const Polynomial& other) const { return Polynomial(*this) *= other; }

		class EachView final {
		public:
			// ReSharper disable once CppNonExplicitConvertingConstructor
			EachView(Polynomial& poly) : _poly(poly) {} // NOLINT(*-explicit-constructor)

			Polynomial& operator+=(const Number& value) {
				for (Number& c : _poly._coeffs) {
					c += value;
				}
				return _poly;
			}
			Polynomial& operator-=(const Number& value) {
				for (Number& c : _poly._coeffs) {
					c -= value;
				}
				return _poly;
			}
			Polynomial& operator*=(const Number& value) { return _poly *= value; }
			Polynomial& operator/=(const Number& value) { return _poly /= value; }

			Polynomial operator+(const Number& value) const { return _poly.copy().each() += value; }
			Polynomial operator-(const Number& value) const { return _poly.copy().each() -= value; }
			Polynomial operator*(const Number& value) const { return _poly * value; }
			Polynomial operator/(const Number& value) const { return _poly / value; }

			Polynomial& operator+=(const Polynomial& other) { return _poly += other; }
			Polynomial& operator-=(const Polynomial& other) { return _poly -= other; }

			Polynomial operator+(const Polynomial& other) const { return _poly + other; }
			Polynomial operator-(const Polynomial& other) const { return _poly - other; }

		private:
			Polynomial& _poly;
		};

		class ConstEachView final {
		public:
			// ReSharper disable once CppNonExplicitConvertingConstructor
			ConstEachView(const Polynomial& poly) : _poly(poly) {} // NOLINT(*-explicit-constructor)

			Polynomial operator+(const Number& value) const { return _poly.copy().each() += value; }
			Polynomial operator-(const Number& value) const { return _poly.copy().each() -= value; }
			Polynomial operator*(const Number& value) const { return _poly * value; }
			Polynomial operator/(const Number& value) const { return _poly / value; }

			Polynomial operator+(const Polynomial& other) const { return _poly + other; }
			Polynomial operator-(const Polynomial& other) const { return _poly - other; }

		private:
			const Polynomial& _poly;
		};

		EachView each() {
			return EachView(*this);
		}
		ConstEachView each() const {
			return ConstEachView(*this);
		}

		using Iterator = typename std::vector<Number>::iterator;
		using ConstIterator = typename std::vector<Number>::const_iterator;
		using ReverseIterator = typename std::vector<Number>::reverse_iterator;
		using ConstReverseIterator = typename std::vector<Number>::const_reverse_iterator;

		Iterator begin() {
			return _coeffs.begin();
		}
		Iterator end() {
			return _coeffs.end();
		}
		ConstIterator begin() const {
			return _coeffs.begin();
		}
		ConstIterator end() const {
			return _coeffs.end();
		}

		ReverseIterator rbegin() {
			return _coeffs.rbegin();
		}
		ReverseIterator rend() {
			return _coeffs.rend();
		}
		ConstReverseIterator rbegin() const {
			return _coeffs.rbegin();
		}
		ConstReverseIterator rend() const {
			return _coeffs.rend();
		}

	private:
		std::vector<Number> _coeffs;
	};
}