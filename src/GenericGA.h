#ifndef _GENERIC_GA_H
#define _GENERIC_GA_H

#include <vector>
#include <functional>
#include <math.h>
#include "random.h"
#include "Utils.h"

template <typename T>
struct Evaluated {
	double score;
	T value;
};

template <typename T>
struct GAState {
	int generation;
	std::vector<T> population;
};

template <typename T>
Evaluated<T> __attribute__((pure)) minimumEvaluated(std::vector<Evaluated<T>> v) {
	return fold1([](Evaluated<T> best, Evaluated<T> current) {
		if (current.score <= best.score) {
			return current;
		}
		return best;
	}, v);
}

template <typename T>
Evaluated<T> __attribute__((pure)) makeEvaluated(T value, double score) {
	Evaluated<T> res;
	res.value = value;
	res.score = score;
	return res;
}

template <typename T, typename F, typename G>
std::function<RNGFUNC(std::vector<Evaluated<T>>)(std::vector<Evaluated<T>>)>
	lambdaPlusN
			( F fitness
			, G mutation
			, int n
			) {

	static_assert(std::is_convertible<F, std::function<double(T)>> ::value,
			"lamdaPlusN's fitness function must be of type T -> double");
	static_assert(std::is_convertible<G, std::function<RNGFUNC(T)(T)>> ::value,
			"lambdaPlusN's mutation function must be of type T -> RNGFUNC(T)");

	return [=](std::vector<Evaluated<T>> population) {
		return bind(mapM(mutation, replicate(n, population[0].value)),
				[=](std::vector<T> newIndividuals) {
		    newIndividuals.insert(newIndividuals.begin(), population[0].value);
			auto fitnesses = map(fitness, newIndividuals);
			auto evaluated = zipWith(makeEvaluated<T>, newIndividuals, fitnesses);
            return pure(std::vector<Evaluated<T>>{ minimumEvaluated(evaluated) });
		});
	};
}

template <typename T, typename F>
std::function<RNGFUNC(GAState<T>)(GAState<T>)> __attribute__((pure)) makeGAFunction
    (F strategy) {

	static_assert(std::is_convertible<F, std::function<RNGFUNC(std::vector<T>)(std::vector<T>)>> ::value,
			"makeGAFunction's strategy function must be of type std::vector<T> -> RNGFUNC(std::vector<T>)");

	return [strategy](GAState<T> state) {
		return bind
				( strategy(state.population)
				, [=](std::vector<T> newPopulation) {
			GAState<T> newState;
			newState.population = newPopulation;
			newState.generation = state.generation + 1;
			return pure(newState);
		});
	};
}

#endif
