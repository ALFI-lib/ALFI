#pragma once

#include <random>
#include <vector>

#include <benchmark/benchmark.h>

#include <ALFI/util/arrays.h>
#include <ALFI/util/stat.h>

template <typename T>
T random_value(T min, T max) {
	static std::mt19937_64 generator(std::random_device{}());

	if constexpr (std::is_integral_v<T>) {
		static std::uniform_int_distribution<T> dist(min, max);
		return dist(generator);
	} else if constexpr (std::is_floating_point_v<T>) {
		static std::uniform_real_distribution<T> dist(min, max);
		return dist(generator);
	} else {
		static_assert(false, "Unsupported type: must be integral or floating-point");
		return {};
	}
}

template <typename T, typename F>
std::vector<T> apply_func(const std::vector<T>& input, F func) {
	std::vector<T> result(input.size());
	std::transform(input.begin(), input.end(), result.begin(), func);
	return result;
}

template <typename T, typename V1, typename V2>
std::vector<T> concat_vectors(V1&& a, V2&& b) {
	std::vector<T> result = std::forward<V1>(a);
	result.reserve(std::size(a) + std::size(b));
	result.insert(result.end(), std::make_move_iterator(std::begin(b)), std::make_move_iterator(std::end(b)));
	return result;
}

template <typename T>
void add_error_metrics(benchmark::State& state, const std::vector<T>& result, const std::vector<T>& expected) {
	const std::vector<double> error = alfi::util::arrays::sub(result, expected);
	state.counters["1.ME"] = alfi::util::stat::mean(error);
	state.counters["2.MAE"] = alfi::util::stat::mean_abs(error);
	state.counters["3.RMSE"] = alfi::util::stat::rms(error);
	state.counters["4.Variance"] = alfi::util::stat::var(error);
}