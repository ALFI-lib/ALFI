#include <gtest/gtest.h>

#include <ALFI.h>

using alfi::poly::Polynomial;

TEST(PolynomialsTest, Creation) {
	const Polynomial<double> p1;
	EXPECT_EQ(p1, Polynomial<double>({}));
	const std::vector<double> v {0, 1, 2};
	const Polynomial<double> p2 (v);
	EXPECT_EQ(p2, Polynomial<double>({0, 1, 2}));
	const Polynomial<double> p3 (std::vector<double>({0, 1, 2}));
	EXPECT_EQ(p3, Polynomial<double>({0, 1, 2}));
	const std::array<double, 3> a {0, 1, 2};
	const Polynomial<double> p4 (a);
	EXPECT_EQ(p4, Polynomial<double>({0, 1, 2}));
	const Polynomial<double> p5 (std::array<double, 3>({0, 1, 2}));
	EXPECT_EQ(p5, Polynomial<double>({0, 1, 2}));
	const Polynomial<double> p6 (v.begin(), v.end());
	EXPECT_EQ(p6, Polynomial<double>({0, 1, 2}));
}

TEST(PolynomialsTest, Empty) {
	Polynomial<double> p = {};
	EXPECT_EQ(p.size(), 0);
	EXPECT_EQ(p.degree(), 0);
	p = +p;
	EXPECT_EQ(p.size(), 0);
	EXPECT_EQ(p.degree(), 0);
	p = -p;
	EXPECT_EQ(p.size(), 0);
	EXPECT_EQ(p.degree(), 0);
}

TEST(PolynomialsTest, General) {
	// x^2 + 2x + 3
	Polynomial<double> p = {1, 2, 3};
	EXPECT_EQ(p.size(), 3);
	EXPECT_EQ(p.degree(), 2);
	EXPECT_EQ(p.max_degree(), 2);
	EXPECT_EQ(p.at_degree(0), 3);
	EXPECT_EQ(p.at_degree(1), 2);
	EXPECT_EQ(p.at_degree(2), 1);
	EXPECT_EQ(p.at_index(0), 1);
	EXPECT_EQ(p.at_index(1), 2);
	EXPECT_EQ(p.at_index(2), 3);

	// (x^2 + 2x + 3) + (7x + 3) = x^2 + 9x + 6
	p += {7, 3};
	EXPECT_EQ(p, Polynomial<double>({1, 9, 6}));

	// (x^2 + 9x + 6) * 3 = 3x^2 + 27x + 18
	p.each() *= 3;
	EXPECT_EQ(p, Polynomial<double>({3, 27, 18}));

	// (3x^2 + 27x + 18) - (-x^3 + 5x^2 - 8) = x^3 - 2x^2 + 27x + 26
	p -= {-1, 5, 0, -8};
	EXPECT_EQ(p, Polynomial<double>({1, -2, 27, 26}));

	// (x^3 - 2x^2 + 27x + 26) * (x + 2) = x^4 + 23x^2 + 80x + 52
	p *= {1.0, 2.0};
	EXPECT_EQ(p, Polynomial<double>({1, 0, 23, 80, 52}));

	// x^4 + 23x^2 + 80x + 52 = 0x^5 + x^4 + 23x^2 + 80x + 52
	p.set_degree(5);
	EXPECT_EQ(p.degree(), 5);
	EXPECT_EQ(p.max_degree(), 4);
	EXPECT_EQ(p, Polynomial<double>({1, 0, 23, 80, 52}));
	EXPECT_EQ(p, Polynomial<double>({0, 1, 0, 23, 80, 52}));

	// 0x^5 + x^4 + 23x^2 + 80x + 52 => 80x + 52
	p.set_degree(1);
	EXPECT_EQ(p.degree(), 1);
	EXPECT_EQ(p.max_degree(), 1);
	EXPECT_EQ(p, Polynomial<double>({80, 52}));
	EXPECT_EQ(p, Polynomial<double>({0, 80, 52}));

	p = {15};
	p.set_degree(1);
	EXPECT_EQ(p.size(), 2);
	EXPECT_EQ(p.degree(), 1);
	EXPECT_EQ(p, Polynomial<double>({0, 15}));

	p = {};
	p.set_degree(1);
	EXPECT_EQ(p.size(), 2);
	EXPECT_EQ(p.degree(), 1);
	EXPECT_EQ(p, Polynomial<double>({0, 0}));

	p = {};
	p.set_degree(2);
	EXPECT_EQ(p.size(), 3);
	EXPECT_EQ(p.degree(), 2);
	EXPECT_EQ(p, Polynomial<double>({0, 0, 0}));

	p = {};
	p.set_size(2);
	EXPECT_EQ(p.size(), 2);
	EXPECT_EQ(p.degree(), 1);
	EXPECT_EQ(p, Polynomial<double>({0, 0}));
}

TEST(PolynomialsTest, SingleValueArithmetic) {
	Polynomial<double> p = {3, 4, 5};
	p += 2;
	EXPECT_EQ(p, Polynomial<double>({3, 4, 7}));
	p -= 3;
	EXPECT_EQ(p, Polynomial<double>({3, 4, 4}));
	p *= 2;
	EXPECT_EQ(p, Polynomial<double>({6, 8, 8}));
	p /= 10;
	EXPECT_EQ(p, Polynomial<double>({0.6, 0.8, 0.8}));
	p = {};
	p += 2;
	EXPECT_EQ(p, Polynomial<double>({2}));
}

TEST(PolynomialsTest, ElementWiseSingleValueArithmetic) {
	Polynomial<double> p = {3, 4, 5};
	p.each() += 2;
	EXPECT_EQ(p, Polynomial<double>({5, 6, 7}));
	p.each() -= 3;
	EXPECT_EQ(p, Polynomial<double>({2, 3, 4}));
	p.each() *= 2;
	EXPECT_EQ(p, Polynomial<double>({4, 6, 8}));
	p.each() /= 10;
	EXPECT_EQ(p, Polynomial<double>({0.4, 0.6, 0.8}));
	p = {};
	p.each() += 2;
	EXPECT_EQ(p, Polynomial<double>({}));
}

TEST(PolynomialsTest, SelfOperators) {
	Polynomial<double> p = {1, 2, 3};
	p.each() += p;
	EXPECT_EQ(p, Polynomial<double>({2, 4, 6}));
	p.each() -= p;
	EXPECT_EQ(p, Polynomial<double>({0, 0, 0}));
}

TEST(PolynomialsTest, UnaryOperators) {
	// ReSharper disable once CppLocalVariableMayBeConst
	Polynomial<double> p1 = {1, 2, 3};
	EXPECT_EQ(+p1, p1);
	EXPECT_EQ(-+-p1, p1);
	EXPECT_EQ(-p1, Polynomial<double>({-1, -2, -3}));

	const Polynomial<double> p1_1 = +-+-p1 += {1, 0, 2};
	EXPECT_EQ(p1_1, Polynomial<double>({2, 2, 5}));

	const Polynomial<double> p2 = {3, 2, 1};
	EXPECT_EQ(+p2, p2);
	EXPECT_EQ(-+-p2, p2);
	EXPECT_EQ(-p2, Polynomial<double>({-3, -2, -1}));

	const Polynomial<double> p2_1 = +-+-p2 += {1, 0, 2};
	EXPECT_EQ(p2_1, Polynomial<double>({4, 2, 3}));
}

TEST(PolynomialsTest, Convolution) {
	Polynomial<double> p1 = {1, -2, 3};
	Polynomial<double> p2 = {4, -1};
	EXPECT_EQ(p1 * p2, Polynomial<double>({4, -9, 14, -3}));
	p1 = {};
	p2 = {1};
	EXPECT_EQ(p1 * p2, Polynomial<double>({}));
	p1 = {3, 4};
	p2 = {-1};
	EXPECT_EQ(p1 * p2, Polynomial<double>({-3, -4}));
}

TEST(PolynomialsTest, Differentiation) {
	EXPECT_EQ(Polynomial().derivative(), Polynomial());
	EXPECT_EQ(Polynomial({0}).derivative(), Polynomial({0}));
	EXPECT_EQ(Polynomial({15}).derivative(), Polynomial({0}));

	const Polynomial<double> p = {4, -3, 2, 5, 1, -1}; // 4x^5 - 3x^4 + 2x^3 + 5x^2 + x - 1
	EXPECT_EQ(p.derivative(), Polynomial<double>({20, -12, 6, 10, 1}));
	EXPECT_EQ(p.derivative(2), Polynomial<double>({80, -36, 12, 10}));
}

TEST(PolynomialsTest, Integration) {
	EXPECT_EQ(Polynomial().integral(5), Polynomial{5.0});
	EXPECT_EQ(Polynomial().integral(5, 2), Polynomial({5.0, 5.0}));

	EXPECT_EQ(Polynomial({0.0}).integral(5), Polynomial{5.0});
	EXPECT_EQ(Polynomial({-2.0}).integral(5), Polynomial({-2.0, 5.0}));
	EXPECT_EQ(Polynomial({0.0}).integral(5, 2), Polynomial({5.0, 5.0}));
	EXPECT_EQ(Polynomial({-2.0}).integral(5, 2), Polynomial({-1.0, 5.0, 5.0}));

	const Polynomial<double> p = {4, -3, 2, 5, 1, -1}; // 4x^5 - 3x^4 + 2x^3 + 5x^2 + x - 1
	EXPECT_EQ(p.integral(200), Polynomial({4.0/6.0, -3.0/5.0, 0.5, 5.0/3.0, 0.5, -1, 200}));
	EXPECT_EQ(p.integral(200, 2), Polynomial({4.0/6.0/7.0, -3.0/5.0/6.0, 0.1, 5.0/3.0/4.0, 0.5/3.0, -1/2.0, 200, 200}));
	EXPECT_EQ(p.integral({1, 2}), Polynomial({4.0/6.0/7.0, -3.0/5.0/6.0, 0.1, 5.0/3.0/4.0, 0.5/3.0, -1/2.0, 2, 1}));
}

TEST(PolynomialsTest, Normalization) {
	EXPECT_EQ(Polynomial().normalize(), Polynomial({0.0}));
	EXPECT_EQ(Polynomial({0.0}).normalize(), Polynomial({0.0}));
	EXPECT_EQ(Polynomial({0.0, 1}).normalize(), Polynomial({1.0}));
}

TEST(PolynomialsTest, PaddingAndCropping) {
	Polynomial<double> p = {0, 1, 2};
	p.pad_left(3, 2);
	EXPECT_EQ(p, Polynomial<double>({2, 2, 2, 0, 1, 2}));
	p.pad_right(2);
	EXPECT_EQ(p, Polynomial<double>({2, 2, 2, 0, 1, 2, 0, 0}));
	p.crop_left(1);
	EXPECT_EQ(p, Polynomial<double>({2, 2, 0, 1, 2, 0, 0}));
	p.crop_right(5);
	EXPECT_EQ(p, Polynomial<double>({2, 2}));
	p.crop_left(100);
	EXPECT_EQ(p, Polynomial<double>({}));

	p = {};
	p.pad_left(1);
	EXPECT_EQ(p.size(), 1);
	p = {};
	p.pad_left(2);
	EXPECT_EQ(p.size(), 2);
}

TEST(PolynomialsTest, Replacing) {
	Polynomial<double> p = {1, NAN, 0, -5, NAN, NAN, 2, NAN};
	p.replace_nan();
	EXPECT_EQ(p, Polynomial<double>({1, 0, 0, -5, 0, 0, 2, 0}));
	p = {1, NAN, 0, -5, NAN, NAN, 2, NAN};
	p.replace_nan(100);
	EXPECT_EQ(p, Polynomial<double>({1, 100, 0, -5, 100, 100, 2, 100}));
	p = {INFINITY, -1.0 * INFINITY, 0, 1, 2};
	p.replace_inf(8);
	EXPECT_EQ(p, Polynomial<double>({8, 8, 0, 1, 2}));
}

TEST(PolynomialsTest, Evaluation) {
	Polynomial<double> p = {};
	EXPECT_EQ(p.eval(0), 0);
	EXPECT_EQ(p.eval(1), 0);
	EXPECT_EQ(p.eval(-100), 0);
	EXPECT_EQ(p.eval(NAN), 0);
	EXPECT_EQ(p.eval(INFINITY), 0);

	p = {0};
	EXPECT_EQ(p.eval(0), 0);
	EXPECT_EQ(p.eval(1), 0);
	EXPECT_EQ(p.eval(-100), 0);
	EXPECT_EQ(p.eval(NAN), 0);
	EXPECT_EQ(p.eval(INFINITY), 0);

	p = {5};
	EXPECT_EQ(p.eval(0), 5);
	EXPECT_EQ(p.eval(1), 5);
	EXPECT_EQ(p.eval(-100), 5);
	EXPECT_EQ(p.eval(NAN), 5);
	EXPECT_EQ(p.eval(INFINITY), 5);

	p = {-1, 0};
	EXPECT_EQ(p.eval(0), 0);
	EXPECT_EQ(p.eval(1), -1);
	EXPECT_EQ(p.eval(-100), 100);
	EXPECT_TRUE(std::isnan(p.eval(NAN)));
	EXPECT_EQ(p.eval(INFINITY), -INFINITY);

	p = {2, -3};
	EXPECT_EQ(p.eval(0), -3);
	EXPECT_EQ(p.eval(1), -1);
	EXPECT_EQ(p.eval(-100), -203);
	EXPECT_TRUE(std::isnan(p.eval(NAN)));
	EXPECT_EQ(p.eval(INFINITY), INFINITY);

	p = {1, 2, 1};
	EXPECT_EQ(p.eval(0), 1);
	EXPECT_EQ(p.eval(1), 4);
	EXPECT_EQ(p.eval(-100), 10000 - 200 + 1);
	EXPECT_TRUE(std::isnan(p.eval(NAN)));
	EXPECT_EQ(p.eval(INFINITY), INFINITY);
	EXPECT_EQ(p.eval({-2, -1, 0, 1, 2}), std::vector<double>({1, 0, 1, 4, 9}));

	EXPECT_EQ(p(50), p.eval(50));
	EXPECT_EQ(p(-1000), p.eval(-1000));
	EXPECT_EQ(p({-50, 0, 1, 2, 10}), p.eval({-50, 0, 1, 2, 10}));
}