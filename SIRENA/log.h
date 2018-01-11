
#ifndef LOG_H
#define LOG_H

#undef log_trace
#undef log_debug
#undef log_info
#undef log_warning
#undef log_error
#undef log_fatal

#define LOG

#ifndef LOG

#define log_trace(fmt, ...) ((void)0)
#define log_debug(fmt, ...) ((void)0)
#define log_info(fmt, ...) ((void)0)
#define log_warning(fmt, ...) ((void)0)
#define log_error(fmt, ...) ((void)0)
#define log_fatal(fmt, ...) ((void)0)

#else

#define log_trace log::trace
#define log_debug log::debug
#define log_info log::info
#define log_warning log::warning
#define log_error log::error
#define log_fatal log::fatal

namespace log
{
  enum level
  {
    TRACE = 0,
    DEBUG = 1,
    INFO = 2,
    WARNING = 3,
    ERR = 4,
    FATAL = 5,
  };
  void trace(const char* format, ...);
  void debug(const char* format, ...);
  void info(const char* format, ...);
  void warning(const char* format, ...);
  void error(const char* format, ...);
  void fatal(const char* format, ...);
}

#endif

#endif
