#ifndef UNIQUE_RESOURCE_HPP
#define UNIQUE_RESOURCE_HPP

#include <concepts>
#include <type_traits>
#include <utility>

// TODO: look at this for reference
// #include <experimental/scope>

template <typename R>
struct default_resource;

template <typename R>
inline constexpr auto default_resource_v = default_resource<R>::value;

template <typename R>
concept resource =
    !std::is_const_v<R> && !std::is_reference_v<R> && std::movable<R> &&
    requires { default_resource_v<R>; } &&
    std::constructible_from<R, decltype(default_resource_v<R>)> &&
    std::equality_comparable_with<R, decltype(default_resource_v<R>)>;

// FIXME: this might be incomplete, especially for reference types
template <typename D, typename R>
concept deleter =
    requires(D &d, R &r) { d(r); } && !std::is_const_v<D> && std::movable<D>;

// TODO: we generally need to probably fix the behaviour when constructing the
// deleter throws

template <resource R, deleter<R> D>
class [[nodiscard]] Unique_resource
{
public:
    constexpr Unique_resource() noexcept(
        std::is_nothrow_constructible_v<R, decltype(default_resource_v<R>)> &&
        std::is_nothrow_default_constructible_v<D>)
        requires std::default_initializable<D>
        : m_resource {default_resource_v<R>}, m_deleter {}
    {
    }

    template <typename RR>
        requires(!std::same_as<std::remove_cvref_t<RR>, Unique_resource>) &&
                    std::constructible_from<R, RR &&> &&
                    std::default_initializable<D>
    explicit constexpr Unique_resource(RR &&resource) noexcept(
        std::is_nothrow_constructible_v<R, RR &&> &&
        std::is_nothrow_default_constructible_v<D>)
        : m_resource {std::forward<RR>(resource)}, m_deleter {}
    {
    }

    template <typename RR, typename DD>
        requires std::constructible_from<R, RR &&> &&
                     std::constructible_from<D, DD &&>
    constexpr Unique_resource(RR &&resource, DD &&deleter) noexcept(
        std::is_nothrow_constructible_v<R, RR &&> &&
        std::is_nothrow_constructible_v<D, DD &&>)
        : m_resource {std::forward<RR>(resource)},
          m_deleter {std::forward<DD>(deleter)}
    {
    }

    // TODO: proper noexcept
    constexpr Unique_resource(Unique_resource &&other) noexcept
        : m_resource {std::exchange(other.m_resource,
                                    R(default_resource_v<R>))},
          m_deleter {std::move(other.m_deleter)}
    {
    }

    // TODO: proper noexcept
    constexpr Unique_resource &operator=(Unique_resource &&other) noexcept
    {
        if (this == &other)
        {
            return *this;
        }

        if (m_resource != default_resource_v<R>)
        {
            m_deleter(m_resource);
        }
        m_resource = std::exchange(other.m_resource, R(default_resource_v<R>));

        if constexpr (std::is_reference_v<D>)
        {
            m_deleter = other.m_deleter;
        }
        else
        {
            m_deleter = std::move(other.m_deleter);
        }

        return *this;
    }

    Unique_resource(const Unique_resource &) = delete;
    Unique_resource &operator=(const Unique_resource &) = delete;

    // TODO: proper noexcept
    constexpr ~Unique_resource() noexcept
    {
        reset();
    }

    // TODO: proper noexcept
    constexpr void reset() noexcept
    {
        if (m_resource != default_resource_v<R>)
        {
            m_deleter(m_resource);
            m_resource = R(default_resource_v<R>);
        }
    }

    // TODO: see cppreference.com for the complete requirements and behaviour
    template <typename RR>
        requires std::constructible_from<R, RR &&>
    constexpr void reset(RR &&resource) noexcept
    {
        if (m_resource != default_resource_v<R>)
        {
            m_deleter(m_resource);
        }
        m_resource = std::forward<RR>(resource);
    }

    // FIXME: return
    constexpr void release() noexcept
    {
        m_resource = R(default_resource_v<R>);
    }

    [[nodiscard]] constexpr const R &get() const noexcept
    {
        return m_resource;
    }

    [[nodiscard]] constexpr const D &get_deleter() const noexcept
    {
        return m_deleter;
    }

    [[nodiscard]] constexpr operator bool() const noexcept
    {
        return m_resource != default_resource_v<R>;
    }

private:
    [[no_unique_address]] R m_resource;
    [[no_unique_address]] D m_deleter;
};

template <class RR, class DD>
Unique_resource(RR &&, DD &&)
    -> Unique_resource<std::remove_cvref_t<RR>, std::remove_cvref_t<DD>>;

#endif