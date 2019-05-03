#ifndef VECTOR_SPACE
#define VECTOR_SPACE

#include <iostream>
#include <tuple>
#include <iostream>
#include <math.h>

namespace VectorSpaceUtils{
	template<typename... VectorTypes> class ProductSpace;

	template<typename... VectorTypes>
	std::ostream& operator<<(std::ostream& os, const ProductSpace<VectorTypes...>& thing);

	template<typename first, typename... VectorTypes>
	struct printHelper{
		static inline void formatPrint(std::ostream& os, const first& item, const VectorTypes&... others){
			os << item << ", ";
			printHelper<VectorTypes...>::formatPrint(os, others...);
		}
	};

	template<typename last>
	struct printHelper<last>{
		static inline void formatPrint(std::ostream& os, const last& item){
			os << item;
		}
	};

	template<typename ...VectorTypes>
	class ProductSpace{
		typedef std::tuple<VectorTypes...> data;
		typedef ProductSpace<VectorTypes...> pType;
		data values;
		static inline auto seq(){
			return std::index_sequence_for<VectorTypes...>();
		}

		static inline void init(VectorTypes&... valsList){
			(void(valsList = VectorTypes()), ...);
		}
		template<size_t... Is>
		static inline void initTuple(data& input,
			std::index_sequence<Is...> seq)
		{
			init(std::get<Is>(input)...);
		}

		static inline void init(VectorTypes&... valsList, const VectorTypes&... initValues){
			(void(valsList = VectorTypes(initValues)), ...);
		}
		template<size_t... Is>
		static inline void initTuple(data& input, const VectorTypes&... initValues,
			std::index_sequence<Is...> seq)
		{
			init(std::get<Is>(input)..., initValues...);
		}

		template<size_t... Is>
		static inline void toStreamTuple(std::ostream& os, const data& thing,
			std::index_sequence<Is...> seq)
		{
			os << "<";
			printHelper<VectorTypes...>::formatPrint(os, std::get<Is>(thing)...);
			os << ">";
		}

		static inline void add(VectorTypes&... left, const VectorTypes&... right){
			(void(left += right),...);
		}
		template<size_t... Is>
		static inline void addTuple(data& left, const data& right,
			std::index_sequence<Is...> seq)
		{
			add(std::get<Is>(left)..., std::get<Is>(right)...);
		}

		static void sub(VectorTypes&... left, const VectorTypes&... right){
			(void(left -= right),...);
		}
		template<size_t... Is>
		static inline void subTuple(data& left, const data& right,
			std::index_sequence<Is...> seq)
		{
			sub(std::get<Is>(left)..., std::get<Is>(right)...);
		}

		static void mult(VectorTypes&... left, double right){
			(void(left *= right),...);
		}
		template<size_t... Is>
		static inline void multTuple(data& left, double right,
			std::index_sequence<Is...> seq)
		{
			mult(std::get<Is>(left)..., right);
		}

		static void div(VectorTypes&... left, double right){
			(void(left /= right),...);
		}
		template<size_t... Is>
		static inline void divTuple(data& left, double right,
			std::index_sequence<Is...> seq)
		{
			div(std::get<Is>(left)..., right);
		}

		static double dot(const VectorTypes&... left, const VectorTypes&... right){
			return ((left * right) + ...);
		}

		template<size_t... Is>
		static inline double dotTuple(const data& left, const data& right,
			std::index_sequence<Is...> seq)
		{
			return dot(std::get<Is>(left)..., std::get<Is>(right)...);
		}
	public:
		ProductSpace(){
			initTuple(values, seq());
		}

		ProductSpace(const VectorTypes&... initValues){
			initTuple(values, initValues..., seq());
		}

		template<size_t first, size_t... Is>
		friend struct Projection;

		friend std::ostream& operator<<<VectorTypes...>(std::ostream& os, const pType& thing);
		
		pType& operator+=(const pType& other){
			addTuple(values, other.values, seq());
			return *this;
		}

		pType& operator-=(const pType& other){
			subTuple(values, other.values, seq());
			return *this;
		}

		pType& operator*=(double other){
			multTuple(values, other, seq());
			return *this;
		}

		pType& operator*=(float other){
			multTuple(values, other, seq());
			return *this;
		}

		pType& operator*=(int other){
			multTuple(values, other, seq());
			return *this;
		}

		template<typename T>
		pType& operator/=(T other){
			divTuple(values, other, seq());
			return *this;
		}

		double defaultDot(const pType& other) const{
			return dotTuple(values, other.values, seq());
		}

		double defaultSquareMagnitude() const{
			return defaultDot(*this);
		}

		double defaultMagnitude() const
		{
			return std::sqrt(defaultSquareMagnitude());
		}

		data getData() const
		{
			return values;
		}
	};

	template<>
	class ProductSpace<>{
		ProductSpace<> operator+=(const ProductSpace<>& other){
			return *this;
		}
		ProductSpace<> operator-=(const ProductSpace<>& other){
			return *this;
		}
		ProductSpace<> operator*=(double other){
			return *this;
		}
		ProductSpace<> operator/=(double other){
			return *this;
		}
	};

	template<typename... VectorTypes>
	std::ostream& operator<<(std::ostream& os, const ProductSpace<VectorTypes...>& thing){
		ProductSpace<VectorTypes...>::toStreamTuple(os, thing.values, ProductSpace<VectorTypes...>::seq());
		return os;
	}

	template<typename... VectorTypes>
	ProductSpace<VectorTypes...> operator+(ProductSpace<VectorTypes...> left, const ProductSpace<VectorTypes...>& right){
		return left += right;
	}

	template<typename... VectorTypes>
	ProductSpace<VectorTypes...> operator-(ProductSpace<VectorTypes...> left, const ProductSpace<VectorTypes...>& right){
		return left -= right;
	}

	template<typename T, typename... VectorTypes>
	ProductSpace<VectorTypes...> operator*(ProductSpace<VectorTypes...> left, T right){
		return left *= right;
	}

	template<typename T, typename... VectorTypes>
	ProductSpace<VectorTypes...> operator*(T left, ProductSpace<VectorTypes...> right){
		return right *= left;
	}

	template<typename... VectorTypes>
	double operator*(ProductSpace<VectorTypes...> left, ProductSpace<VectorTypes...> right){
		return left.defualtDot(right);
	}

	template<typename T, typename... VectorTypes>
	ProductSpace<VectorTypes...> operator/(ProductSpace<VectorTypes...> left, T right){
		return left /= right;
	}

	template<size_t first, size_t... Is>
	struct Projection{
		template<typename... VectorTypes>
		static auto get(ProductSpace<VectorTypes...> input){
			return Projection<Is...>::get(std::get<first>(input.values));
		}
	};

	template<size_t last>
	struct Projection<last>{
		template<typename... VectorTypes>
		static auto get(ProductSpace<VectorTypes...> input){
			return std::get<last>(input.values);
		}
	};

	template<size_t... Is, typename... VectorTypes>
	auto project(ProductSpace<VectorTypes...> input){
		return Projection<Is...>::get(input);
	}
}

namespace vsu = VectorSpaceUtils;

#endif
