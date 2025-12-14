#ifndef UNIQUE_HANDLE_HPP
#define UNIQUE_HANDLE_HPP

#include <concepts>
#include <type_traits>
#include <utility>

// TODO: allow reference types by wrapping them in a std::reference_wrapper
// TODO: check exception guarantees when intermediate operations throw

template <typename H, typename D, auto null_handle = H {}>
    requires std::equality_comparable_with<H, decltype(null_handle)> &&
             (std::invocable<D &, H &> || std::invocable<D &>) &&
             std::assignable_from<H &, decltype(null_handle)> &&
             (!std::is_reference_v<H>) && (!std::is_reference_v<D>)
class [[nodiscard]] Unique_handle
{
public:
    constexpr Unique_handle() noexcept(
        std::is_nothrow_constructible_v<H, decltype(null_handle)> &&
        std::is_nothrow_default_constructible_v<D>)
        requires std::constructible_from<H, decltype(null_handle)> &&
                     std::default_initializable<D>
        : m_handle {null_handle}, m_deleter {}
    {
    }

    template <typename HH>
        requires(!std::same_as<std::remove_cvref_t<HH>, Unique_handle>) &&
                    std::constructible_from<H, HH &&> &&
                    std::default_initializable<D>
    explicit constexpr Unique_handle(HH &&handle) noexcept(
        std::is_nothrow_constructible_v<H, HH &&> &&
        std::is_nothrow_default_constructible_v<D>)
        : m_handle {std::forward<HH>(handle)}, m_deleter {}
    {
    }

    template <typename HH, typename DD>
        requires std::constructible_from<H, HH &&> &&
                     std::constructible_from<D, DD &&>
    constexpr Unique_handle(HH &&handle, DD &&deleter) noexcept(
        std::is_nothrow_constructible_v<H, HH &&> &&
        std::is_nothrow_constructible_v<D, DD &&>)
        : m_handle {std::forward<HH>(handle)},
          m_deleter {std::forward<DD>(deleter)}
    {
    }

    constexpr Unique_handle(Unique_handle &&other) noexcept(
        std::is_nothrow_move_constructible_v<H> &&
        std::is_nothrow_assignable_v<H &, decltype(null_handle)> &&
        std::is_nothrow_move_constructible_v<D>)
        requires std::move_constructible<H> && std::move_constructible<D>
        : m_handle {std::exchange(other.m_handle, null_handle)},
          m_deleter {std::move(other.m_deleter)}
    {
    }

    constexpr Unique_handle &operator=(Unique_handle &&other) noexcept(
        noexcept(reset()) && noexcept(release()) &&
        std::is_nothrow_move_assignable_v<D>)
        requires std::assignable_from<H &, H> && std::move_constructible<H> &&
                 std::assignable_from<D &, D>
    {
        if (this == &other)
        {
            return *this;
        }

        reset(other.release());

        m_deleter = std::move(other.m_deleter);

        return *this;
    }

    Unique_handle(const Unique_handle &) = delete;
    Unique_handle &operator=(const Unique_handle &) = delete;

    constexpr ~Unique_handle() noexcept(noexcept(reset()))
    {
        reset();
    }

    constexpr void reset() noexcept(
        noexcept(std::declval<const H &>() !=
                 std::declval<decltype(null_handle)>()) &&
        (std::is_invocable_v<D &, H &> ? std::is_nothrow_invocable_v<D &, H &>
                                       : std::is_nothrow_invocable_v<D &>) &&
        std::is_nothrow_assignable_v<H &, decltype(null_handle)>)
    {
        if (m_handle != null_handle)
        {
            if constexpr (std::is_invocable_v<D &, H &>)
            {
                m_deleter(m_handle);
            }
            else
            {
                m_deleter();
            }
            m_handle = null_handle;
        }
    }

    template <typename HH>
        requires std::assignable_from<H &, HH &&>
    constexpr void
    reset(HH &&handle) noexcept(noexcept(reset()) &&
                                std::is_nothrow_assignable_v<H &, HH &&>)
    {
        reset();
        m_handle = std::forward<HH>(handle);
    }

    template <typename HH, typename DD>
        requires std::assignable_from<H &, HH &&> &&
                 std::assignable_from<D &, DD &&>
    constexpr void reset(HH &&handle, DD &&deleter) noexcept(
        noexcept(reset()) && std::is_nothrow_assignable_v<H &, HH &&> &&
        std::is_nothrow_assignable_v<D &, DD &&>)
    {
        reset();
        m_handle = std::forward<HH>(handle);
        m_deleter = std::forward<DD>(deleter);
    }

    constexpr H
    release() noexcept(std::is_nothrow_move_constructible_v<H> &&
                       std::is_nothrow_assignable_v<H &, decltype(null_handle)>)
        requires std::move_constructible<H>
    {
        return std::exchange(m_handle, null_handle);
    }

    [[nodiscard]] constexpr const H &get() const noexcept
    {
        return m_handle;
    }

    [[nodiscard]] constexpr const D &get_deleter() const noexcept
    {
        return m_deleter;
    }

    [[nodiscard]] constexpr bool has_value() const
        noexcept(noexcept(std::declval<const H &>() !=
                          std::declval<decltype(null_handle)>()))
    {
        return m_handle != null_handle;
    }

private:
    [[no_unique_address]] H m_handle;
    [[no_unique_address]] D m_deleter;
};

template <class HH, class DD>
Unique_handle(HH &&, DD &&)
    -> Unique_handle<std::remove_cvref_t<HH>, std::remove_cvref_t<DD>>;

#endif