#ifndef UNIQUE_HANDLE_HPP
#define UNIQUE_HANDLE_HPP

#include <concepts>
#include <type_traits>
#include <utility>

template <typename H, typename D, auto null_handle = H {}>
    requires std::move_constructible<H> && (!std::is_reference_v<H>) &&
             (std::invocable<D &, H &> || std::invocable<D &>) &&
             std::constructible_from<H, decltype(null_handle)> &&
             std::assignable_from<H &, decltype(null_handle)> &&
             requires(const H &h) {
                 { h != null_handle } -> std::convertible_to<bool>;
             }
class [[nodiscard]] Unique_handle
{
public:
    constexpr Unique_handle() noexcept
        requires std::default_initializable<D>
        : m_handle {null_handle}, m_deleter {}
    {
    }

    template <typename HH>
        requires(!std::same_as<std::remove_cvref_t<HH>, Unique_handle>) &&
                    std::constructible_from<H, HH &&> &&
                    std::default_initializable<D>
    explicit constexpr Unique_handle(HH &&handle) noexcept
        : m_handle {std::forward<HH>(handle)}, m_deleter {}
    {
    }

    template <typename HH, typename DD>
        requires std::constructible_from<H, HH &&> &&
                     std::constructible_from<D, DD &&>
    constexpr Unique_handle(HH &&handle, DD &&deleter) noexcept
        : m_handle {std::forward<HH>(handle)},
          m_deleter {std::forward<DD>(deleter)}
    {
    }

    constexpr Unique_handle(Unique_handle &&other) noexcept
        requires std::move_constructible<D>
        : m_handle {std::exchange(other.m_handle, null_handle)},
          m_deleter {std::move(other.m_deleter)}
    {
    }

    constexpr Unique_handle &operator=(Unique_handle &&other) noexcept
        requires std::is_move_assignable_v<D> && (!std::is_reference_v<D>)
    {
        reset(other.release());
        m_deleter = std::move(other.m_deleter);
        return *this;
    }

    Unique_handle(const Unique_handle &) = delete;
    Unique_handle &operator=(const Unique_handle &) = delete;

    constexpr ~Unique_handle() noexcept
    {
        reset();
    }

    constexpr void reset() noexcept
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
    constexpr void reset(HH &&handle) noexcept
    {
        reset();
        m_handle = std::forward<HH>(handle);
    }

    constexpr H release() noexcept
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

    [[nodiscard]] constexpr bool has_value() const noexcept
    {
        return m_handle != null_handle;
    }

private:
    H m_handle;
    [[no_unique_address]] D m_deleter;
};

template <typename HH, typename DD>
Unique_handle(HH &&, DD &&)
    -> Unique_handle<std::remove_cvref_t<HH>, std::remove_cvref_t<DD>>;

#endif