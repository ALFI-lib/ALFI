#include <benchmark/benchmark.h>

#include <ALFI/spline/cubic.h>

#include "../bench_utils.h"
#include "../bench_data.h"

const std::vector<alfi::spline::CubicSpline<>::Type> cubic_spline_types = {
	alfi::spline::CubicSpline<>::Types::Natural{},
	alfi::spline::CubicSpline<>::Types::NotAKnot{},
	alfi::spline::CubicSpline<>::Types::ParabolicEnds{},
};

static void BM_cubic_spline(benchmark::State& state) {
	const auto& [func_name, lambda] = funcs_and_ints[state.range(0)].first;
	const auto& interval = funcs_and_ints[state.range(0)].second;
	const auto& dist_type = dists[state.range(1)];
	const auto& n = state.range(2);
	const auto& type = cubic_spline_types[state.range(3)];
	const std::vector<double> X = of_type<double>(dist_type, n, interval.first, interval.second);
	const std::vector<double> Y = apply_func(X, lambda);
	alfi::spline::CubicSpline<> spline;
	for (auto _ : state) {
		spline = alfi::spline::CubicSpline<>(X, Y, type);
		benchmark::DoNotOptimize(spline);
		benchmark::ClobberMemory();
	}
}
BENCHMARK(BM_cubic_spline)
	->ArgsProduct({
		benchmark::CreateDenseRange(0, static_cast<int64_t>(funcs_and_ints.size()) - 1, 1),
		benchmark::CreateDenseRange(0, static_cast<int64_t>(dists.size()) - 1, 1),
		benchmark::CreateRange(1 << 6, 1 << 21, 1 << 3),
		benchmark::CreateDenseRange(0, static_cast<int64_t>(cubic_spline_types.size()) - 1, 1),
	})
	->ArgNames({"func", "dist", "n", "type"});

static void BM_cubic_spline_values(benchmark::State& state) {
	const auto& [func_name, lambda] = funcs_and_ints[state.range(0)].first;
	const auto& interval = funcs_and_ints[state.range(0)].second;
	const auto& dist_type = dists[state.range(1)];
	const auto& n = state.range(2);
	const auto& nn = state.range(3);
	const auto& type = cubic_spline_types[state.range(4)];
	const std::vector<double> X = of_type<double>(dist_type, n, interval.first, interval.second);
	const std::vector<double> Y = apply_func(X, lambda);
	const std::vector<double> xx = of_type<double>(alfi::dist::Type::UNIFORM, nn, interval.first, interval.second);
	const alfi::spline::CubicSpline<> spline(X, Y, type);
	std::vector<double> result;
	for (auto _ : state) {
		result = spline.eval(xx);
		benchmark::DoNotOptimize(result);
		benchmark::ClobberMemory();
	}
	const std::vector<double> expected = apply_func(xx, lambda);
	add_error_metrics(state, result, expected);
}
BENCHMARK(BM_cubic_spline_values)
	->ArgsProduct({
		benchmark::CreateDenseRange(0, static_cast<int64_t>(funcs_and_ints.size()) - 1, 1),
		benchmark::CreateDenseRange(0, static_cast<int64_t>(dists.size()) - 1, 1),
		concat_vectors<int64_t>(benchmark::CreateRange(8, 32, 2), benchmark::CreateRange(1 << 6, 1 << 21, 1 << 3)),
		{nn},
		benchmark::CreateDenseRange(0, static_cast<int64_t>(cubic_spline_types.size()) - 1, 1),
	})
	->ArgNames({"func", "dist", "n", "nn", "type"});

BENCHMARK_MAIN();