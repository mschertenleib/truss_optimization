#ifndef VEC2_HPP
#define VEC2_HPP

#include <cmath>

struct vec2
{
    float x;
    float y;

    [[nodiscard]] constexpr bool
    operator==(const vec2 &) const noexcept = default;
    [[nodiscard]] constexpr bool
    operator!=(const vec2 &) const noexcept = default;
};

struct vec3
{
    float x;
    float y;
    float z;

    [[nodiscard]] constexpr bool
    operator==(const vec3 &) const noexcept = default;
    [[nodiscard]] constexpr bool
    operator!=(const vec3 &) const noexcept = default;
};

struct vec4
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

constexpr vec2 &operator+=(vec2 &u, const vec2 &v) noexcept
{
    u = u + v;
    return u;
}

constexpr vec2 &operator-=(vec2 &u, const vec2 &v) noexcept
{
    u = u - v;
    return u;
}

constexpr vec2 &operator*=(vec2 &u, const vec2 &v) noexcept
{
    u = u * v;
    return u;
}

constexpr vec2 &operator/=(vec2 &u, const vec2 &v) noexcept
{
    u = u / v;
    return u;
}

constexpr vec2 &operator+=(vec2 &v, float f) noexcept
{
    v = v + f;
    return v;
}

constexpr vec2 &operator-=(vec2 &v, float f) noexcept
{
    v = v - f;
    return v;
}

constexpr vec2 &operator*=(vec2 &v, float f) noexcept
{
    v = v * f;
    return v;
}

constexpr vec2 &operator/=(vec2 &v, float f) noexcept
{
    v = v / f;
    return v;
}

[[nodiscard]] constexpr float dot(const vec2 &u, const vec2 &v) noexcept
{
    return u.x * v.x + u.y * v.y;
}

[[nodiscard]] constexpr float cross(const vec2 &u, const vec2 &v) noexcept
{
    return u.x * v.y - u.y * v.x;
}

[[nodiscard]] inline float norm(const vec2 &v) noexcept
{
    return std::sqrt(dot(v, v));
}

[[nodiscard]] inline vec2 normalize(const vec2 &v) noexcept
{
    return v * (1.0f / norm(v));
}

#endif // VEC2_HPP
