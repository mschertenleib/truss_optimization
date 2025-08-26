#ifndef UNIQUE_RESOURCE_HPP
#define UNIQUE_RESOURCE_HPP

#include <concepts>
#include <type_traits>
#include <utility>

// TODO: look at this for reference
// #include <experimental/scope>

template <typename R>
concept resource =
    !std::is_const_v<R> && !std::is_reference_v<R> && std::movable<R>;

// FIXME: this might be incomplete, especially for reference types
template <typename D, typename R>
concept deleter =
    requires(D &d, R &r) { d(r); } && !std::is_const_v<D> && std::movable<D>;

template <resource R, deleter<R> D>
class [[nodiscard]] Unique_resource
{
public:
    constexpr Unique_resource() noexcept(
        std::is_nothrow_default_constructible_v<R> &&
        std::is_nothrow_default_constructible_v<D>)
        requires std::default_initializable<R> && std::default_initializable<D>
        : m_resource {}, m_deleter {}, m_owns_resource {false}
    {
    }

    template <typename RR>
        requires(!std::same_as<std::remove_cvref_t<RR>, Unique_resource>) &&
                    std::constructible_from<R, RR &&> &&
                    std::default_initializable<D>
    explicit constexpr Unique_resource(RR &&resource) noexcept(
        std::is_nothrow_constructible_v<R, RR &&> &&
        std::is_nothrow_default_constructible_v<D>)
        : m_resource {std::forward<RR>(resource)},
          m_deleter {},
          m_owns_resource {true}
    {
    }

    template <typename RR, typename DD>
        requires std::constructible_from<R, RR &&> &&
                     std::constructible_from<D, DD &&>
    constexpr Unique_resource(RR &&resource, DD &&deleter) noexcept(
        std::is_nothrow_constructible_v<R, RR &&> &&
        std::is_nothrow_constructible_v<D, DD &&>)
        : m_resource {std::forward<RR>(resource)},
          m_deleter {std::forward<DD>(deleter)},
          m_owns_resource {true}
    {
    }

    // FIXME: does this work with a reference deleter?
    // TODO: proper noexcept
    constexpr Unique_resource(Unique_resource &&other)
        : m_resource {std::move(other.m_resource)},
          m_deleter {std::move(other.m_deleter)},
          m_owns_resource {std::exchange(other.m_owns_resource, false)}
    {
    }

    // FIXME: check behaviour
    // TODO: proper noexcept
    constexpr Unique_resource &operator=(Unique_resource &&other)
    {
        if (this == &other)
        {
            return *this;
        }

        if (m_owns_resource)
        {
            m_deleter(m_resource);
        }
        m_resource = std::move(other.m_resource);

        if constexpr (std::is_reference_v<D>)
        {
            m_deleter = other.m_deleter;
        }
        else
        {
            m_deleter = std::move(other.m_deleter);
        }

        m_owns_resource = std::exchange(other.m_owns_resource, false);

        return *this;
    }

    Unique_resource(const Unique_resource &) = delete;
    Unique_resource &operator=(const Unique_resource &) = delete;

    constexpr ~Unique_resource() noexcept
    {
        reset();
    }

    constexpr void reset() noexcept
    {
        if (m_owns_resource)
        {
            m_deleter(m_resource);
            m_owns_resource = false;
        }
    }

    // TODO: proper noexcept
    template <typename RR>
        requires std::assignable_from<R &, RR &&>
    constexpr void reset(RR &&resource)
    {
        reset();
        m_resource = std::forward<RR>(resource);
        m_owns_resource = true;
    }

    constexpr void release() noexcept
    {
        m_owns_resource = false;
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