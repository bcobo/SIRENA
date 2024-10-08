/**********************************************************************
*   This software is part of the grant PID2021-122955OB-C41
*   and by 'ERDF A way of making Europe'.
*
***********************************************************************
*                      LOG
*
*  File:       log.h
*  Developers: Beatriz Cobo
* 	           cobo@ifca.unican.es
*              IFCA
*              Maite Ceballos
*              ceballos@ifca.unican.es
*              IFCA
*                                                                     
***********************************************************************/

#ifndef LOG_H
#define LOG_H

#undef log_test
#undef log_trace
#undef log_debug
#undef log_info
#undef log_warning
#undef log_error
#undef log_fatal

#define LOG
#define TEST_LOG

#ifndef TEST_LOG

#define log_test(fmt, ...) ((void)0)

#else

#define log_test tlog::test

namespace tlog 
{
  void test(const char* format, ...);
}

#endif

#ifndef LOG

#define log_trace(fmt, ...) ((void)0)
#define log_debug(fmt, ...) ((void)0)
#define log_info(fmt, ...) ((void)0)
#define log_warning(fmt, ...) ((void)0)
#define log_error(fmt, ...) ((void)0)
#define log_fatal(fmt, ...) ((void)0)

#else

#define log_trace slog::trace
#define log_debug slog::debug
#define log_info slog::info
#define log_warning slog::warning
#define log_error slog::error
#define log_fatal slog::fatal

namespace slog
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
