#ifndef UNIQUE_RESOURCE_HPP
#define UNIQUE_RESOURCE_HPP

#include <concepts>
#include <type_traits>
#include <utility>

// TODO: look at this for reference
// #include <experimental/scope>

// TODO: requirements on R and D
template <typename R, typename D>
class [[nodiscard]] Unique_resource
{
public:
    constexpr Unique_resource() noexcept(
        std::is_nothrow_default_constructible_v<R> &&
        std::is_nothrow_default_constructible_v<D>)
        requires(std::is_default_constructible_v<R> &&
                 std::is_default_constructible_v<D>)
        : m_resource {}, m_deleter {}, m_owns_resource {false}
    {
    }

    // TODO: noexcept based on type traits
    // TODO: see cppreference.com for the complete requirements and behaviour
    template <typename RR, typename DD = D>
        requires std::constructible_from<R, RR &&> &&
                     std::constructible_from<D, DD &&> &&
                     (!std::same_as<std::remove_cvref_t<RR>, Unique_resource>)
    explicit constexpr Unique_resource(RR &&resource, DD &&deleter = DD())
        : m_resource {std::forward<RR>(resource)},
          m_deleter {std::forward<DD>(deleter)},
          m_owns_resource {true}
    {
    }

    // TODO: see cppreference.com for the complete requirements
    // TODO: we also need to take into account the possibility that constructing
    // the resource or deleter might throw (see cppreference.com)
    constexpr Unique_resource(Unique_resource &&rhs) noexcept
        : m_resource {std::move(rhs.m_resource)},
          m_deleter {std::move(rhs.m_deleter)},
          m_owns_resource {rhs.m_owns_resource}
    {
        rhs.m_owns_resource = false;
    }

    // TODO: see cppreference.com for the complete requirements and behaviour
    constexpr Unique_resource &operator=(Unique_resource &&rhs) noexcept
    {
        Unique_resource tmp(std::move(rhs));
        swap(tmp);
        return *this;
    }

    Unique_resource(const Unique_resource &) = delete;
    Unique_resource &operator=(const Unique_resource &) = delete;

    // TODO: noexcept?
    constexpr ~Unique_resource() noexcept
    {
        reset();
    }

    [[nodiscard]] constexpr const R &get() const noexcept
    {
        return m_resource;
    }

    [[nodiscard]] constexpr const D &get_deleter() const noexcept
    {
        return m_deleter;
    }

    constexpr void reset() noexcept
    {
        if (m_owns_resource)
        {
            m_deleter(m_resource);
            m_owns_resource = false;
        }
    }

    // TODO: see cppreference.com for the complete requirements and behaviour
    template <typename RR>
        requires std::constructible_from<R, RR &&>
    constexpr void reset(RR &&resource) noexcept
    {
        reset();
        m_resource = std::forward<RR>(resource);
        m_owns_resource = true;
    }

    constexpr void release() noexcept
    {
        m_owns_resource = false;
    }

    // TODO: noexcept based on type traits?
    constexpr void swap(Unique_resource &rhs) noexcept
    {
        using std::swap;
        swap(m_resource, rhs.m_resource);
        swap(m_deleter, rhs.m_deleter);
        swap(m_owns_resource, rhs.m_owns_resource);
    }

    [[nodiscard]] constexpr operator bool() const noexcept
    {
        return m_owns_resource;
    }

private:
    [[no_unique_address]] R m_resource;
    [[no_unique_address]] D m_deleter;
    bool m_owns_resource;
};

template <class RR, class DD>
Unique_resource(RR &&, DD &&)
    -> Unique_resource<std::remove_cvref_t<RR>, std::remove_cvref_t<DD>>;

#endif