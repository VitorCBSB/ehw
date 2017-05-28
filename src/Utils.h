/*
 * Utils.h
 *
 *  Created on: 12/11/2016
 *      Author: Vitor Coimbra
 */

#ifndef SRC_UTILS_H_
#define SRC_UTILS_H_

#include <functional>
#include <algorithm>
#include <tuple>
#include <vector>

template <typename T, typename U>
struct Pair {
	T first;
	U second;
};

template <typename T, typename U>
Pair<T, U>  makePair(T t, U u) {
	Pair<T, U> p;
	p.first = t;
	p.second = u;
	return p;
}

template <typename T>
std::vector<std::vector<T>> chunksOf(std::vector<T> v, unsigned int length) {
    std::vector<std::vector<T>> result;

    for (unsigned int i = 0; i < v.size(); i = i + length) {
        std::vector<T> tempRes;
        for (unsigned int j = 0; j < length && (i + j) < v.size(); j++) {
            tempRes.push_back(v[i + j]);
        }
        result.emplace_back(tempRes);
    }

    return result;
}

template <typename T>
std::vector<T>  replicate(int times, T element) {
	std::vector<T> result;
	for (int i = 0; i < times; i++) {
		result.emplace_back(element);
	}
	return result;
}

template <typename T>
std::vector<T> concat(std::vector<std::vector<T>> vs) {
	std::vector<T> result;

	for (auto v : vs) {
		result.insert(result.end(), v.begin(), v.end());
	}

	return result;
}

template <typename T, typename F>
auto  map
	( F f
    , std::vector<T> va
	) -> std::vector<decltype(f(va[0]))> {

	using U = decltype(f(va[0]));
	static_assert(std::is_convertible<F, std::function<U(T)>> ::value,
			"map's function must be of type T -> U");

	std::vector<U> result;
	for (auto a : va) {
		result.push_back(f(a));
	}
	return result;
}

template <typename T>
std::vector<std::vector<T>>  transpose(std::vector<std::vector<T>> v) {
	std::vector<std::vector<T>> res(v[0].size(), std::vector<T>(v.size()));

	for (unsigned int i = 0; i < v[0].size(); i++) {
		for (unsigned int j = 0; j < v.size(); j++) {
			res[i][j] = v[j][i];
		}
	}

	return res;
}

template <typename T, typename U>
std::vector<std::tuple<T, U>> zip(std::vector<T> a, std::vector<U> b) {
	std::vector<std::tuple<T, U>> result;

	auto resultSize = std::min(a.size(), b.size());
	for (unsigned int i = 0; i < resultSize; i++) {
		result.emplace_back(std::make_tuple(a[i], b[i]));
	}

	return result;
}

// Non-inclusive on to.
std::vector<unsigned int>  vectorFromTo(unsigned int from, unsigned int to) {
	std::vector<unsigned int> res;
	for (unsigned int i = from; i < to; i++) {
		res.emplace_back(i);
	}
	return res;
}

template <typename T>
std::vector<T>  tail(std::vector<T> v) {
	v.erase(v.begin());
	return v;
}

template <typename F, typename T>
std::vector<T>  filter
	( F f
    , std::vector<T> vs
	) {

	static_assert(std::is_convertible<F, std::function<bool(T)>> ::value,
			"filter's function must be of type T -> bool");

	std::vector<T> result;

	for (auto v : vs) {
		if (f(v)) {
			result.emplace_back(v);
		}
	}

	return result;
}

template <typename T, typename U>
U  foldWhile
	( std::function<bool(U)> cond
	, std::function<U(U, T)> f
	, U init
	, std::vector<T> in
	) {
    auto res = init;
    unsigned cur = 0;

    while (cur < in.size() && cond(res)) {
        res = f(res, in[cur]);
        cur++;
    }
    return res;
}

template <typename T, typename U, typename F>
U  fold
	( F f
	, U init
	, std::vector<T> in
	) {

	static_assert(std::is_convertible<F, std::function<U(U, T)>> ::value,
			"fold's function must be of type U -> T -> U");

	unsigned int cur = 0;
	auto res = init;

	while (cur < in.size()) {
		res = f(res, in[cur]);
		cur++;
	}
	return res;
}

// 'in' vector MUST be non-empty.
template <typename T, typename F>
auto  fold1
	( F f
	, std::vector<T> in
	) -> decltype(f(in[0], in[0])) {

	static_assert(std::is_convertible<F, std::function<T(T,T)>> ::value,
			"fold1's function must be of type T -> T -> T");

	return fold<T, T>(f, in[0], tail(in));
}

template <typename T, typename F, typename G>
T  iterateWhile
	( F cond
	, G iter
	, T init
	) {

	static_assert(std::is_convertible<F, std::function<bool(T)>> ::value,
			"iterateWhile's conditional function must be of type T -> bool");
	static_assert(std::is_convertible<G, std::function<T(T)>> ::value,
			"iterateWhile's iteration function must be of type T -> T");

	auto res = init;

	while (cond(res)) {
		res = iter(res);
	}
	return res;
}

template <typename T, typename U, typename F>
auto
 zipWith
	( F f
    , std::vector<T> tv
	, std::vector<U> uv
	) -> std::vector<decltype(f(tv[0], uv[0]))>{

	using V = decltype(f(tv[0], uv[0]));
	static_assert(std::is_convertible<F, std::function<V(T, U)>> ::value,
			"zipWith's function must be of type T -> U -> V");

	std::vector<V> res;
	for (int i = 0; i < std::min(tv.size(), uv.size()); i++) {
		res.push_back(f(tv[i], uv[i]));
	}
	return res;
}

template <typename T, typename F>
bool  any
	( F f
	, std::vector<T> vs
    ) {

	static_assert(std::is_convertible<F, std::function<bool(T)>> ::value,
			"any's function must be of type T -> bool");

	return fold(std::logical_or<bool>(), false, map(f, vs));
}

#endif /* SRC_UTILS_H_ */
