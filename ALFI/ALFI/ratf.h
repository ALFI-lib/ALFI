#pragma once

#include <cmath>

#include "config.h"
#include "util/numeric.h"
#include "util/poly.h"

/**
	@namespace alfi::ratf
	@brief Namespace providing support for rational functions.

	This namespace provides types and functions for representing, computing and evaluating rational functions of the form
	\f[
		f(x) = \frac{A(x)}{B(x)}
	\f]
	where \f(A(x)\f) and \f(B(x)\f) are polynomials.
*/
namespace alfi::ratf {
	/**
		@brief Represents a rational function \f(\displaystyle f(x) = \frac{A(x)}{B(x)}\f), where \f(A(x)\f) and \f(B(x)\f) are polynomials.

		A pair (`std::pair`) of polynomials `{.first as numerator, .second as denominator}` stored as containers of coefficients in descending degree order.
	*/
	template <typename Number = DefaultNumber, template <typename, typename...> class Container = DefaultContainer>
	using RationalFunction = std::pair<Container<Number>, Container<Number>>;

	/**
		@brief Evaluates the rational function at a scalar point using Horner's method.

		Computes \f(f(x) = \frac{A(x)}{B(x)}\f) by evaluating the values of the numerator \f(A\f) and the denominator \f(B\f)
		using Horner's method, and then dividing the results.

		@ref val utilizes this function for \f(|x| \leq 1\f).
	*/
	template <typename Number = DefaultNumber, template <typename, typename...> class Container = DefaultContainer>
	std::enable_if_t<!traits::has_size<Number>::value, Number>
	val_mul(const RationalFunction<Number, Container>& rf, const Number& x) {
		Number n = 0;
		for (const auto& c : rf.first) {
			n = n * x + c;
		}
		Number d = 0;
		for (const auto& c : rf.second) {
			d = d * x + c;
		}
		return n / d;
	}

	/**
		@brief Evaluates the rational function at each point in the container using @ref val_mul for scalar values.
	*/
	template <typename Number = DefaultNumber, template <typename, typename...> class Container = DefaultContainer>
	std::enable_if_t<traits::has_size<Container<Number>>::value, Container<Number>>
	val_mul(const RationalFunction<Number, Container>& rf, const Container<Number>& xx) {
		Container<Number> result(xx.size());
#if defined(_OPENMP) && !defined(ALFI_DISABLE_OPENMP)
#pragma omp parallel for
#endif
		for (SizeT i = 0; i < xx.size(); ++i) {
			result[i] = val_mul(rf, xx[i]);
		}
		return result;
	}

	/**
		@brief Evaluates the rational function at a scalar point by factoring out powers of x.

		Computes \f(f(x) = \frac{A(x)}{B(x)}\f) by evaluating both numerator and denominator in reverse order,
		effectively factoring out the dominant power of \f(x\f) to improve numerical stability for large \f(|x|\f).

		@ref val utilizes this function for \f(|x| > 1\f).
	*/
	template <typename Number = DefaultNumber, template <typename, typename...> class Container = DefaultContainer>
	std::enable_if_t<!traits::has_size<Number>::value, Number>
	val_div(const RationalFunction<Number, Container>& rf, const Number& x) {
		const auto& numerator = rf.first;
		const auto& denominator = rf.second;
		Number n = 0;
		for (auto i = numerator.rbegin(); i != numerator.rend(); ++i) {
			n = n / x + *i;
		}
		Number d = 0;
		for (auto i = denominator.rbegin(); i != denominator.rend(); ++i) {
			d = d / x + *i;
		}
		const auto numerator_degree = numerator.empty() ? 0 : numerator.size() - 1;
		const auto denominator_degree = denominator.empty() ? 0 : denominator.size() - 1;
		if (numerator_degree >= denominator_degree) {
			return n / d * util::numeric::binpow(x, numerator_degree - denominator_degree);
		} else {
			return n / d / util::numeric::binpow(x, denominator_degree - numerator_degree);
		}
	}

	/**
		@brief Evaluates the rational function at each point in the container using @ref val_div for scalar values.
	*/
	template <typename Number = DefaultNumber, template <typename, typename...> class Container = DefaultContainer>
	std::enable_if_t<traits::has_size<Container<Number>>::value, Container<Number>>
	val_div(const RationalFunction<Number, Container>& rf, const Container<Number>& xx) {
		Container<Number> result(xx.size());
#if defined(_OPENMP) && !defined(ALFI_DISABLE_OPENMP)
#pragma omp parallel for
#endif
		for (SizeT i = 0; i < xx.size(); ++i) {
			result[i] = val_div(rf, xx[i]);
		}
		return result;
	}

	/**
		@brief Evaluates the rational function at a scalar point.

		Calls @ref val_mul for \f(|x| \leq 1\f) and @ref val_div otherwise for the sake of numerical stability.
	*/
	template <typename Number = DefaultNumber, template <typename, typename...> class Container = DefaultContainer>
	std::enable_if_t<!traits::has_size<Number>::value, Number>
	val(const RationalFunction<Number, Container>& rf, const Number& x) {
		if (std::abs(x) <= 1) {
			return val_mul(rf, x);
		} else {
			return val_div(rf, x);
		}
	}

	/**
		@brief Evaluates the rational function at each point in the container using @ref val for scalar values.
	*/
	template <typename Number = DefaultNumber, template <typename, typename...> class Container = DefaultContainer>
	std::enable_if_t<traits::has_size<Container<Number>>::value, Container<Number>>
	val(const RationalFunction<Number, Container>& rf, const Container<Number>& xx) {
		Container<Number> result(xx.size());
#if defined(_OPENMP) && !defined(ALFI_DISABLE_OPENMP)
#pragma omp parallel for
#endif
		for (SizeT i = 0; i < xx.size(); ++i) {
			result[i] = val(rf, xx[i]);
		}
		return result;
	}

	/**
		@brief Computes the [@p n / @p m] Pade approximant of the polynomial @p P about the point \f(x = 0\f).

		The Pade approximant is given by the formula
		\f[
			P(x) = \frac{A(x)}{B(x)} = \frac{\sum_{i=0}^{n}{a_ix^i}}{b_0+\sum_{j=0}^{m}{b_jx^j}}
		\f]
		where \f(A(x)\f) and \f(B(x)\f) are the numerator and denominator polynomials, respectively; \f(b_0 \neq 0\f).

		This function computes the Pade approximant of a polynomial by applying a modified extended Euclidean algorithm for polynomials.

		The modification consists in that:
		- the algorithm may terminate early if the numerator's degree already meets the requirements,
		- or perform an extra iteration involving a division by a zero polynomial in special cases.

		The latter is necessary to avoid false negatives, for example, when computing the `[2/2]` approximant of the function \f(x^5\f).

		Without the additional check, this also may lead to false positives, as in the case of computing the `[2/2]` approximant of \f(x^4\f).@n
		This is prevented by verifying that the constant term of the denominator is non-zero after the algorithm completes.

		@param P the polynomial to approximate (a container of coefficients in descending degree order)
		@param n the maximum degree for the numerator
		@param m the maximum degree for the denominator
		@param epsilon the tolerance used to determine whether a coefficient is considered zero (default is machine epsilon)
		@return a pair `{numerator, denominator}` representing the Pade approximant; if an approximant does not exist, an empty pair is returned
	 */
	template <typename Number = DefaultNumber, template <typename, typename...> class Container = DefaultContainer>
	RationalFunction<Number,Container> pade(Container<Number> P, SizeT n, SizeT m, const Number& epsilon = std::numeric_limits<Number>::epsilon()) {
		if constexpr (std::is_signed_v<SizeT>) {
			if (n < 0 || m < 0) {
				return {{}, {}};
			}
		}

		util::poly::normalize(P);

		Container<Number> Xnm1(n + m + 2, 0);
		Xnm1[0] = 1;

		// Modified extended Euclidean algorithm `gcd(a,b)=as+bt` without s variable
		// a = Xnm1
		// b = P
		Container<Number> old_r, r, old_t, t;
		if (Xnm1.size() >= P.size()) {
			old_r = std::move(Xnm1), r = std::move(P);
			old_t = {0}, t = {1};
		} else {
			old_r = std::move(P), r = std::move(Xnm1);
			old_t = {1}, t = {0};
		}

		// `old_r.size()` strictly decreases, except maybe the first iteration
		// ReSharper disable once CppDFALoopConditionNotUpdated
		while (old_r.size() > n + 1) {
			auto [q, new_r] = util::poly::div(old_r, r, epsilon);

			const auto qt = util::poly::mul(q, t);
			const auto new_t_size = std::max(old_t.size(), qt.size());
			Container<Number> new_t(new_t_size, 0);
			for (SizeT i = 0, offset = new_t_size - old_t.size(); i < old_t.size(); ++i) {
				new_t[offset+i] = old_t[i];
			}
			for (SizeT i = 0, offset = new_t_size - qt.size(); i < qt.size(); ++i) {
				new_t[offset+i] -= qt[i];
			}

			old_r = std::move(r);
			r = std::move(new_r);
			old_t = std::move(t);
			t = std::move(new_t);

			util::poly::normalize(old_r, epsilon);
		}

		util::poly::normalize(old_t, epsilon);

		if (old_t.size() > m + 1 || std::abs(old_t[old_t.size()-1]) <= epsilon) {
			return {{}, {}};
		}

		return std::make_pair(old_r, old_t);
	}
}