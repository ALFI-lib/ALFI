#include <benchmark/benchmark.h>

#include <ALFI/spline/linear.h>

#include "../bench_utils.h"
#include "../bench_data.h"

static void BM_linear_spline(benchmark::State& state) {
	const auto& [func_name, lambda] = funcs_and_ints[state.range(0)].first;
	const auto& interval = funcs_and_ints[state.range(0)].second;
	const auto& dist_type = dists[state.range(1)];
	const auto& n = state.range(2);
	const std::vector<double> X = of_type<double>(dist_type, n, interval.first, interval.second);
	const std::vector<double> Y = apply_func(X, lambda);
	alfi::spline::LinearSpline<> spline;
	for (auto _ : state) {
		spline = alfi::spline::LinearSpline<>(X, Y);
		benchmark::DoNotOptimize(spline);
		benchmark::ClobberMemory();
	}
}
BENCHMARK(BM_linear_spline)
	->ArgsProduct({
		benchmark::CreateDenseRange(0, static_cast<int64_t>(funcs_and_ints.size()) - 1, 1),
		benchmark::CreateDenseRange(0, static_cast<int64_t>(dists.size()) - 1, 1),
		benchmark::CreateRange(1 << 6, 1 << 21, 1 << 3),
	})
	->ArgNames({"func", "dist", "n"});

static void BM_linear_spline_values(benchmark::State& state) {
	const auto& [func_name, lambda] = funcs_and_ints[state.range(0)].first;
	const auto& interval = funcs_and_ints[state.range(0)].second;
	const auto& dist_type = dists[state.range(1)];
	const auto& n = state.range(2);
	const auto& nn = state.range(3);
	const std::vector<double> X = of_type<double>(dist_type, n, interval.first, interval.second);
	const std::vector<double> Y = apply_func(X, lambda);
	const std::vector<double> xx = of_type<double>(alfi::dist::Type::UNIFORM, nn, interval.first, interval.second);
	const alfi::spline::LinearSpline<> spline(X, Y);
	std::vector<double> result;
	for (auto _ : state) {
		result = spline.eval(xx);
		benchmark::DoNotOptimize(result);
		benchmark::ClobberMemory();
	}
	const std::vector<double> expected = apply_func(xx, lambda);
	add_error_metrics(state, result, expected);
}
BENCHMARK(BM_linear_spline_values)
	->ArgsProduct({
		benchmark::CreateDenseRange(0, static_cast<int64_t>(funcs_and_ints.size()) - 1, 1),
		benchmark::CreateDenseRange(0, static_cast<int64_t>(dists.size()) - 1, 1),
		concat_vectors<int64_t>(benchmark::CreateRange(8, 32, 2), benchmark::CreateRange(1 << 6, 1 << 21, 1 << 3)),
		{nn},
	})
	->ArgNames({"func", "dist", "n", "nn"});

BENCHMARK_MAIN();