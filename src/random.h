/*
 * random.h
 *
 *  Created on: 01/11/2016
 *      Author: Vitor Coimbra
 */

#ifndef RANDOM_H_
#define RANDOM_H_

#include <functional>
#include <vector>
#include "Utils.h"

typedef unsigned int random_type;

template <typename T>
struct RNGProduct {
	random_type newRng;
	T product;
};

#define RNGFUNC(T) std::function<RNGProduct<T>(random_type)>

// Periodo maximal de 32 bits, ou seja,
// gera todos os outros numeros antes de voltar para o original.
random_type  rand_(random_type current) {
	auto bit = ((current >> 0) ^ (current >> 10) ^ (current >> 30) ^ (current >> 31)) & 1;
	return (current >> 1) | (bit << 31);
}

template <typename T>
RNGProduct<T>  makeRNGProduct(T value, random_type newRng) {
	RNGProduct<T> res;
	res.product = value;
	res.newRng = newRng;
	return res;
}

template <typename T>
T  evalState(RNGFUNC(T) f, random_type rng) {
	return f(rng).product;
}

template <typename T>
RNGProduct<T>  runRngFunc(RNGFUNC(T) f, random_type rng) {
	return f(rng);
}

template <typename T>
RNGFUNC(T)  pure(T val) {
	return [val](random_type rng) {
		return makeRNGProduct<T>(val, rng);
	};
}

RNGFUNC(random_type)  getRandom() {
	return [](random_type rng) {
		auto next = rand_(rng);
		return makeRNGProduct<random_type>(rng, next);
	};
}

template <typename T, typename F>
auto  bind(RNGFUNC(T) m, F f) -> decltype(f(m(0).product)) {
	using U = decltype(f(m(0).product)(0).product);
	static_assert(std::is_convertible<F, std::function<RNGFUNC(U)(T)>> ::value,
			"bind requires a function type RNGFUNC(U)(T).");

	return [m, f](random_type rng) mutable {
		auto res = m(rng);
		auto b = f(res.product);
		return b(res.newRng);
	};
}

template <typename T, typename F>
auto  operator>>(RNGFUNC(T) m, F f) -> decltype(f(m(0).product)) {
	using U = decltype(f(m(0).product)(0).product);
	static_assert(std::is_convertible<F, std::function<RNGFUNC(U)(T)>> ::value,
			"bind requires a function type RNGFUNC(U)(T).");

	return [m, f](random_type rng) mutable {
		auto res = m(rng);
		auto b = f(res.product);
		return b(res.newRng);
	};
}

template <typename T, typename U>
RNGFUNC(U)  rmap(std::function<U(T)> f, RNGFUNC(T) a) {
	return bind(a, [f](T val) {
		return pure(f(val));
	});
}

// This imperative implementation pains me, but it'll do.
template<typename T>
RNGFUNC(std::vector<T>) sequence(std::vector<RNGFUNC(T)> seq) {
	return [=](random_type rng) {
		std::vector<T> result;
		random_type currentRng = rng;
		for (auto f : seq) {
			auto seqRes = f(currentRng);
			result.emplace_back(seqRes.product);
			currentRng = seqRes.newRng;
		}

		return makeRNGProduct(result, currentRng);
	};
}

template <typename T, typename F>
auto mapM
	( F f
    , std::vector<T> v
	) -> RNGFUNC(std::vector<decltype(f(v[0])(0).product)>) {

	using U = decltype(f(v[0])(0).product);
	static_assert(std::is_convertible<F, std::function<RNGFUNC(U)(T)>> ::value,
			"mapM's function must be of type T -> RNGFUNC(U)");

	return sequence(map(f, v));
}

template <typename T, typename U, typename F>
auto foldM
	( F f
    , U init
	, std::vector<T> ins
	) -> RNGFUNC(U) {

	static_assert(std::is_convertible<F, std::function<RNGFUNC(U)(U, T)>> ::value,
			"foldM's function must be of type U -> T -> RNGFUNC(U)");

	return [=](random_type rand) mutable {
		auto curRand = rand;
		auto curVal = init;

		for (T in : ins) {
			auto res = runRngFunc(f(curVal, in), curRand);
			curRand = res.newRng;
			curVal = res.product;
		}

		return makeRNGProduct(curVal, curRand);
	};
}

template <typename T, typename F, typename G>
auto  iterateWhileM
	( F cond
    , G iter
	, T val
	) -> RNGFUNC(T) {

	static_assert(std::is_convertible<F, std::function<bool(T)>> ::value,
			"iterateWhileM's conditional function must be of type T -> bool");
	static_assert(std::is_convertible<G, std::function<RNGFUNC(T)(T)>> ::value,
			"iterateWhileM's iteration function must be of type T -> RNGFUNC(T)");

	return [=](random_type rand) mutable {
		while (cond(val)) {
			auto res = runRngFunc(iter(val), rand);
			val = res.product;
			rand = res.newRng;
		}
		return makeRNGProduct(val, rand);
	};
}

template <typename C, typename T, typename F>
auto whileFoldM
	( C cond
    , F f
	, std::vector<T> ins
	) -> RNGFUNC(decltype(f(ins[0])(0).product)) {

	using U = decltype(f(ins[0])(0).product);
	static_assert(std::is_convertible<C, std::function<bool(T, U)>> ::value,
			"whileFoldM's conditional function must be of type T -> U -> bool");
	static_assert(std::is_convertible<F, std::function<RNGFUNC(U)(T)>> ::value,
			"whileFoldM's function must be of type T -> RNGFUNC(U)");

	return [=](random_type rand) mutable {
		auto curRand = rand;
		U curVal;

		for (T in : ins) {
		    auto wow = runRngFunc(f(in), curRand);
		    curRand = wow.newRng;
		    curVal = wow.product;
		    if (!cond(in, curVal)) {
		        break;
		    }
		}

		return makeRNGProduct(curVal, curRand);
	};
}

#endif /* RANDOM_H_ */
