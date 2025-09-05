#ifndef VEC_HPP
#define VEC_HPP

#include <cmath>

template <unsigned int N>
struct vec;

using vec2 = vec<2>;
using vec3 = vec<3>;
using vec4 = vec<4>;

template <>
struct vec<2>
{
    float x;
    float y;

    [[nodiscard]] constexpr bool
    operator==(const vec2 &) const noexcept = default;
    [[nodiscard]] constexpr bool
    operator!=(const vec2 &) const noexcept = default;
};

template <>
struct vec<3>
{
    float x;
    float y;
    float z;

    [[nodiscard]] constexpr bool
    operator==(const vec3 &) const noexcept = default;
    [[nodiscard]] constexpr bool
    operator!=(const vec3 &) const noexcept = default;
};

template <>
struct vec<4>
{
    float x;
    float y;
    float z;
    float w;

    [[nodiscard]] constexpr bool
    operator==(const vec4 &) const noexcept = default;
    [[nodiscard]] constexpr bool
    operator!=(const vec4 &) const noexcept = default;
};

[[nodiscard]] constexpr vec2 operator+(const vec2 &v) noexcept
{
    return v;
}

[[nodiscard]] constexpr vec2 operator-(const vec2 &v) noexcept
{
    return {-v.x, -v.y};
}

[[nodiscard]] constexpr vec2 operator+(const vec2 &u, const vec2 &v) noexcept
{
    return {u.x + v.x, u.y + v.y};
}

[[nodiscard]] constexpr vec2 operator-(const vec2 &u, const vec2 &v) noexcept
{
    return {u.x - v.x, u.y - v.y};
}

[[nodiscard]] constexpr vec2 operator*(const vec2 &u, const vec2 &v) noexcept
{
    return {u.x * v.x, u.y * v.y};
}

[[nodiscard]] constexpr vec2 operator/(const vec2 &u, const vec2 &v) noexcept
{
    return {u.x / v.x, u.y / v.y};
}

[[nodiscard]] constexpr vec2 operator+(const vec2 &v, float f) noexcept
{
    return {v.x + f, v.y + f};
}

[[nodiscard]] constexpr vec2 operator-(const vec2 &v, float f) noexcept
{
    return {v.x - f, v.y - f};
}

[[nodiscard]] constexpr vec2 operator*(const vec2 &v, float f) noexcept
{
    return {v.x * f, v.y * f};
}

[[nodiscard]] constexpr vec2 operator/(const vec2 &v, float f) noexcept
{
    return {v.x / f, v.y / f};
}

[[nodiscard]] constexpr vec2 operator+(float f, const vec2 &v) noexcept
{
    return {f + v.x, f + v.y};
}

[[nodiscard]] constexpr vec2 operator-(float f, const vec2 &v) noexcept
{
    return {f - v.x, f - v.y};
}

[[nodiscard]] constexpr vec2 operator*(float f, const vec2 &v) noexcept
{
    return {f * v.x, f * v.y};
}

[[nodiscard]] constexpr vec2 operator/(float f, const vec2 &v) noexcept
{
    return {f / v.x, f / v.y};
}

[[nodiscard]] constexpr vec3 operator+(const vec3 &v) noexcept
{
    return v;
}

[[nodiscard]] constexpr vec3 operator-(const vec3 &v) noexcept
{
    return {-v.x, -v.y, -v.z};
}

[[nodiscard]] constexpr vec3 operator+(const vec3 &u, const vec3 &v) noexcept
{
    return {u.x + v.x, u.y + v.y, u.z + v.z};
}

[[nodiscard]] constexpr vec3 operator-(const vec3 &u, const vec3 &v) noexcept
{
    return {u.x - v.x, u.y - v.y, u.z - v.z};
}

[[nodiscard]] constexpr vec3 operator*(const vec3 &u, const vec3 &v) noexcept
{
    return {u.x * v.x, u.y * v.y, u.z * v.z};
}

[[nodiscard]] constexpr vec3 operator/(const vec3 &u, const vec3 &v) noexcept
{
    return {u.x / v.x, u.y / v.y, u.z / v.z};
}

[[nodiscard]] constexpr vec3 operator+(const vec3 &v, float f) noexcept
{
    return {v.x + f, v.y + f, v.z + f};
}

[[nodiscard]] constexpr vec3 operator-(const vec3 &v, float f) noexcept
{
    return {v.x - f, v.y - f, v.z - f};
}

[[nodiscard]] constexpr vec3 operator*(const vec3 &v, float f) noexcept
{
    return {v.x * f, v.y * f, v.z * f};
}

[[nodiscard]] constexpr vec3 operator/(const vec3 &v, float f) noexcept
{
    return {v.x / f, v.y / f, v.z / f};
}

[[nodiscard]] constexpr vec3 operator+(float f, const vec3 &v) noexcept
{
    return {f + v.x, f + v.y, f + v.z};
}

[[nodiscard]] constexpr vec3 operator-(float f, const vec3 &v) noexcept
{
    return {f - v.x, f - v.y, f - v.z};
}

[[nodiscard]] constexpr vec3 operator*(float f, const vec3 &v) noexcept
{
    return {f * v.x, f * v.y, f * v.z};
}

[[nodiscard]] constexpr vec3 operator/(float f, const vec3 &v) noexcept
{
    return {f / v.x, f / v.y, f / v.z};
}

template <unsigned int N>
constexpr vec<N> &operator+=(vec<N> &u, const vec<N> &v) noexcept
{
    u = u + v;
    return u;
}

template <unsigned int N>
constexpr vec<N> &operator-=(vec<N> &u, const vec<N> &v) noexcept
{
    u = u - v;
    return u;
}

template <unsigned int N>
constexpr vec<N> &operator*=(vec<N> &u, const vec<N> &v) noexcept
{
    u = u * v;
    return u;
}

template <unsigned int N>
constexpr vec<N> &operator/=(vec<N> &u, const vec<N> &v) noexcept
{
    u = u / v;
    return u;
}

template <unsigned int N>
constexpr vec<N> &operator+=(vec<N> &v, float f) noexcept
{
    v = v + f;
    return v;
}

template <unsigned int N>
constexpr vec<N> &operator-=(vec<N> &v, float f) noexcept
{
    v = v - f;
    return v;
}

template <unsigned int N>
constexpr vec<N> &operator*=(vec<N> &v, float f) noexcept
{
    v = v * f;
    return v;
}

template <unsigned int N>
constexpr vec<N> &operator/=(vec<N> &v, float f) noexcept
{
    v = v / f;
    return v;
}

[[nodiscard]] constexpr float dot(const vec2 &u, const vec2 &v) noexcept
{
    return u.x * v.x + u.y * v.y;
}

[[nodiscard]] constexpr float dot(const vec3 &u, const vec3 &v) noexcept
{
    return u.x * v.x + u.y * v.y + u.z * v.z;
}

[[nodiscard]] constexpr float cross(const vec2 &u, const vec2 &v) noexcept
{
    return u.x * v.y - u.y * v.x;
}

[[nodiscard]] constexpr vec3 cross(const vec3 &u, const vec3 &v) noexcept
{
    return {
        u.y * v.z - u.z * v.y, u.z * v.x - u.x * v.z, u.x * v.y - u.y * v.x};
}

template <unsigned int N>
[[nodiscard]] inline float norm(const vec<N> &v) noexcept
{
    return std::sqrt(dot(v, v));
}

template <unsigned int N>
[[nodiscard]] inline vec<N> normalize(const vec<N> &v) noexcept
{
    return v * (1.0f / norm(v));
}

#endif
