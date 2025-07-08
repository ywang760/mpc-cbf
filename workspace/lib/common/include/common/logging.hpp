#pragma once

// #define SPDLOG_FUNCTION __PRETTY_FUNCTION__ // full signature in %!

#include <spdlog/spdlog.h> // pulls in <fmt/format.h>
#include <sstream>
#include <iomanip>
#include <cstdlib>
#include <string>
#include <algorithm>
#include <cctype>

#include <Eigen/Dense>
#include <ginac/ginac.h>

namespace common
{
    //--------------------------------------------------------------------------
    // Environment variable logging configuration
    //--------------------------------------------------------------------------

    namespace detail
    {
        /**
         * @brief Convert string to lowercase for case-insensitive comparison
         */
        inline std::string toLower(const std::string &str)
        {
            std::string result = str;
            std::transform(result.begin(), result.end(), result.begin(),
                           [](unsigned char c)
                           { return std::tolower(c); });
            return result;
        }

        /**
         * @brief Convert string to spdlog level enum
         */
        inline spdlog::level::level_enum stringToLogLevel(const std::string &level_str)
        {
            std::string lower_str = toLower(level_str);

            if (lower_str == "trace")
                return spdlog::level::trace;
            if (lower_str == "debug")
                return spdlog::level::debug;
            if (lower_str == "info")
                return spdlog::level::info;
            if (lower_str == "warn" || lower_str == "warning")
                return spdlog::level::warn;
            if (lower_str == "error" || lower_str == "err")
                return spdlog::level::err;
            if (lower_str == "critical" || lower_str == "crit")
                return spdlog::level::critical;
            if (lower_str == "off")
                return spdlog::level::off;

            // Default to info if unrecognized
            return spdlog::level::info;
        }
    }

    /**
     * @brief Initialize logging level from environment variables
     *
     * Checks the following environment variables in order:
     * 1. SPDLOG_LEVEL
     * 2. LOG_LEVEL
     * 3. LOGLEVEL
     *
     * Supported values (case-insensitive):
     * - trace, debug, info, warn/warning, error/err, critical/crit, off
     *
     * @param default_level Default log level if no environment variable is set
     */
    inline void initializeLoggingFromEnv(spdlog::level::level_enum default_level = spdlog::level::info)
    {
        spdlog::level::level_enum level = default_level;

        // Check environment variables in order of preference
        const char *env_vars[] = {"SPDLOG_LEVEL", "LOG_LEVEL", "LOGLEVEL"};

        for (const char *env_var : env_vars)
        {
            const char *env_value = std::getenv(env_var);
            if (env_value != nullptr)
            {
                level = detail::stringToLogLevel(env_value);
                break;
            }
        }

        // Always set the level (either from env var or default)
        spdlog::set_level(level);
    }
    /**
     * @brief Initialize logging with a custom pattern and default level
     *
     * Sets the logging pattern and initializes the log level from environment variables.
     * Returns the default logger for further use.
     *
     * @param pattern Custom pattern string (optional)
     * @param default_level Default log level if no environment variable is set
     */
    inline std::shared_ptr<spdlog::logger> initializeLogging(
        const std::string &pattern = "[%H:%M:%S] [%^%l%$] [%s:%# %!] %v",
        spdlog::level::level_enum default_level = spdlog::level::info)
    { 
        spdlog::set_pattern(pattern);
        initializeLoggingFromEnv(default_level);
        auto logger = spdlog::default_logger();

        return logger; // Return the default logger for further use
    }

    //--------------------------------------------------------------------------
    // detail::toString – internal helpers reused by log* and formatters
    //--------------------------------------------------------------------------

    namespace detail
    {
        /*----- Eigen matrix ---------------------------------------------------*/
        template <typename Derived>
        std::string toString(const Eigen::MatrixBase<Derived> &M,
                             int precision = 6)
        {
            std::ostringstream ss;
            ss << "(" << M.rows() << "×" << M.cols() << ")\n";
            ss << std::fixed << std::setprecision(precision);

            for (Eigen::Index i = 0; i < M.rows(); ++i)
            {
                ss << "[ ";
                for (Eigen::Index j = 0; j < M.cols(); ++j)
                {
                    ss << std::setw(6 + precision) << M(i, j);
                    if (j + 1 < M.cols())
                        ss << " , ";
                }
                ss << " ]\n";
            }
            return ss.str();
        }

        /*----- GiNaC matrix ---------------------------------------------------*/
        inline std::string toString(const GiNaC::matrix &M,
                                    int precision = 6)
        {
            std::ostringstream ss;
            ss << "(" << M.rows() << "×" << M.cols() << ")\n";
            ss << std::fixed << std::setprecision(precision);

            for (unsigned i = 0; i < M.rows(); ++i)
            {
                ss << "[ ";
                for (unsigned j = 0; j < M.cols(); ++j)
                {
                    const GiNaC::ex &e = M(i, j);
                    if (e.info(GiNaC::info_flags::numeric))
                    {
                        double val = GiNaC::ex_to<GiNaC::numeric>(e).to_double();
                        ss << std::setw(6 + precision) << val;
                    }
                    else
                    {
                        ss << e; // keep symbolic form
                    }
                    if (j + 1 < M.cols())
                        ss << " , ";
                }
                ss << " ]\n";
            }
            return ss.str();
        }
    } // namespace detail

    //--------------------------------------------------------------------------
    // Public helpers: logMatrix / logVector
    //--------------------------------------------------------------------------

    template <typename Derived>
    inline void logMatrix(std::string_view label,
                          const Eigen::MatrixBase<Derived> &M,
                          spdlog::level::level_enum lvl = spdlog::level::info,
                          int precision = 6)
    {
        auto msg = fmt::format("{} {}", label, detail::toString(M, precision));
        spdlog::log(spdlog::source_loc{__FILE__, __LINE__, SPDLOG_FUNCTION},
                    lvl, "{}", msg);
    }

    inline void logMatrix(std::string_view label,
                          const GiNaC::matrix &M,
                          spdlog::level::level_enum lvl = spdlog::level::info,
                          int precision = 6)
    {
        auto msg = fmt::format("{} {}", label, detail::toString(M, precision));
        spdlog::log(spdlog::source_loc{__FILE__, __LINE__, SPDLOG_FUNCTION},
                    lvl, "{}", msg);
    }

    /*----- Vector helper (Eigen column-vector) --------------------------------*/
    template <typename Derived>
    inline void logVector(std::string_view label,
                          const Eigen::MatrixBase<Derived> &v,
                          spdlog::level::level_enum lvl = spdlog::level::info,
                          int precision = 6)
    {
        static_assert(Derived::ColsAtCompileTime == 1 ||
                          Derived::ColsAtCompileTime == Eigen::Dynamic,
                      "logVector() requires a column-vector");

        auto msg = fmt::format("{} {}", label, detail::toString(v, precision));
        spdlog::log(spdlog::source_loc{__FILE__, __LINE__, SPDLOG_FUNCTION},
                    lvl, "{}", msg);
    }

} // namespace common


//--------------------------------------------------------------------------
// fmt::formatter specialisations – let "{}" print matrices/vectors
//--------------------------------------------------------------------------

namespace fmt
{

    // ------------ Eigen matrices *and* column-vectors ------------------------
    template <typename Derived, typename Char>
    struct formatter<
        Derived, Char,
        std::enable_if_t<std::is_base_of_v<Eigen::EigenBase<Derived>, Derived>, void>>
    {
        int precision_ = 6;

        constexpr auto parse(basic_format_parse_context<Char> &ctx)
        {
            auto it = ctx.begin();
            if (it != ctx.end() && *it == '.') // allow "{:.3}"
            {
                ++it;
                precision_ = 0;
                while (it != ctx.end() && std::isdigit(*it))
                {
                    precision_ = 10 * precision_ + (*it - '0');
                    ++it;
                }
            }
            return it;
        }

        template <typename FormatContext>
        auto format(const Derived &M, FormatContext &ctx) const
        {
            return format_to(ctx.out(), "{}", common::detail::toString(M, precision_));
        }
    };

    // ------------ GiNaC matrices ---------------------------------------------
    template <typename Char>
    struct formatter<GiNaC::matrix, Char>
    {
        int precision_ = 6;

        constexpr auto parse(basic_format_parse_context<Char> &ctx)
        {
            auto it = ctx.begin();
            if (it != ctx.end() && *it == '.')
            {
                ++it;
                precision_ = 0;
                while (it != ctx.end() && std::isdigit(*it))
                {
                    precision_ = 10 * precision_ + (*it - '0');
                    ++it;
                }
            }
            return it;
        }

        template <typename FormatContext>
        auto format(const GiNaC::matrix &M, FormatContext &ctx) const
        {
            return format_to(ctx.out(), "{}", common::detail::toString(M, precision_));
        }
    };

} // namespace fmt
