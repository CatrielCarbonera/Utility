#pragma once
#include <type_traits>

template<typename T, typename S, bool isNumeric>
struct NumericAPI;

template<typename T, typename S>
struct NumericAPI<T, S, false>
{
    using type = T;
    using value_type = S;
};

template<typename T, typename S>
struct NumericAPI<T, S, true>
{
    using type = T;
    using value_type = S;

    friend type & operator+=(type & a, type const & b)
    {
        a.get() += b.get();
        return a;
    }

    friend type & operator-=(type & a, type const & b)
    {
        a.get() -= b.get();
        return a;
    }

    friend type & operator *= (type & a, value_type const & b)
    {
        a.get() *= b;
        return a;
    }

    friend type & operator /= (type & a, value_type const & b)
    {
        a.get() /= b;
        return a;
    }

    friend type operator + (type a, type const & b) { return(a += b); }
    friend type operator - (type a, type const & b) { return(a -= b); }
    friend type operator*(type a, value_type const & b) { return(a *= b); }
    friend type operator*(value_type a, type const & b) { return(b *= a); }
    friend type operator / (type a, value_type const & b) { return(a /= b); }
};

template  <typename T, typename V>
struct TypeSafe : public NumericAPI<TypeSafe<T, V>, V, std::is_arithmetic_v<V>>
{
    using value_type = V;

    TypeSafe() = default;
    TypeSafe(TypeSafe const &) = default;
    TypeSafe(TypeSafe &&) = default;

    explicit TypeSafe(V const & v) : m_v(v) { }

    ~TypeSafe() = default;
    value_type & get() { return m_v; }
    value_type const & get() const { return m_v; }

protected:
    value_type m_v;
};

#define DECLARE_SAFE_TYPE(TP, VL) using TP = TypeSafe<struct TP##_s, VL>;
