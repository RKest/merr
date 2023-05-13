#include <iostream>

#include "boost/math/differentiation/autodiff.hpp"
#include "boost/multiprecision/cpp_bin_float.hpp"

namespace d = boost::math::differentiation;

constexpr static auto cw = 4196;
constexpr static auto ck = 896;
constexpr static auto t0 = 273.16;

template<typename T1, typename T2, typename MK, typename MKW, typename MKWL>
d::promote<T1, T2, MK, MKW, MKWL>qt(const T1& t1, const T2 t2, const MK& mk, const MKW& mkw, const MKWL& mkwl)
{
	return ((cw * (mkw - mk) + ck * mk) * (t1 - t2) - cw * (mkwl - mkw) * (t2 - t0)) / mkwl - mkw;
}

int main()
{
	using float50 = boost::multiprecision::cpp_bin_float_50;
	constexpr unsigned max_order = 2;
	auto const vars = d::make_ftuple<float50, max_order, max_order, max_order, max_order, max_order>(296.05, 285.65, 180.280 / 1000, 372.754 / 1000, 401.800 / 1000);
	auto const t1 = std::get<0>(vars);
	auto const t2 = std::get<1>(vars);
	auto const mk = std::get<2>(vars);
	auto const mkw = std::get<3>(vars);
	auto const mkwl = std::get<4>(vars);
	auto const res = qt(t1,t2,mk,mkw,mkwl);

	std::cout << "Result = " << res.derivative(max_order, max_order, max_order, max_order, max_order) << '\n';
}

