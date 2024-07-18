#pragma once

#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>
#include <stdexcept>
#include <memory>
#include <assert.h>
#include <iostream>
#include <sstream>
#include <ctime>


struct LogCheckError {
    LogCheckError() : str(nullptr) {}
    explicit LogCheckError(const std::string& str_p) : str(new std::string(str_p)) {}
    ~LogCheckError() {
        if(str != nullptr) delete str;
    }
    //Type conversion
    operator bool() { return str != nullptr; }
    //The error string
    std::string* str;
};

#define DEFINE_CHECK_FUNC(name, op)                      \
template <typename X, typename Y>                               \
inline LogCheckError LogCheck##name(const X& x, const Y& y) {   \
    if (x op y) return LogCheckError();                         \
    std::ostringstream os;                                      \
    os << " (" << x << " vs. " << y << ") ";                    \
    return LogCheckError(os.str());                             \
}                                                               \
inline LogCheckError LogCheck##name(int x, int y) {             \
    return LogCheck##name<int, int>(x, y);                      \
}
DEFINE_CHECK_FUNC(_LT, <)
DEFINE_CHECK_FUNC(_GT, >)
DEFINE_CHECK_FUNC(_LE, <=)
DEFINE_CHECK_FUNC(_GE, >=)
DEFINE_CHECK_FUNC(_EQ, ==)
DEFINE_CHECK_FUNC(_NE, !=)



//Always on checking
#define CHECK(x)                                         \
if(!(x))                                                        \
    LogMessageFatal(__FILE__, __LINE__).stream()    \
    << "Check failed: " #x << " "

#define SURFELWARP_CHECK_BINARY_OP(name, op, x, y)                         \
if(LogCheckError err = LogCheck##name(x, y))   \
    LogMessageFatal(__FILE__, __LINE__).stream()           \
    << "Check failed: " << #x " " #op " " #y << *(err.str)

#define	CHECK_LT(x, y) CHECK_BINARY_OP(_LT, < , x, y)
#define	CHECK_GT(x, y) CHECK_BINARY_OP(_GT, > , x, y)
#define	CHECK_LE(x, y) CHECK_BINARY_OP(_LE, <=, x, y)
#define	CHECK_GE(x, y) CHECK_BINARY_OP(_GE, >=, x, y)
#define	CHECK_EQ(x, y) CHECK_BINARY_OP(_EQ, ==, x, y)
#define	CHECK_NE(x, y) CHECK_BINARY_OP(_NE, !=, x, y)

//The log type for later use
#define LOG_INFO LogMessage(__FILE__, __LINE__)
#define LOG_ERROR LOG_INFO
#define LOG_WARNING LOG_INFO
#define LOG_FATAL LogMessageFatal(__FILE__, __LINE__)
#define LOG_BEFORE_THROW LogMessage().stream()

//For different severity
#define LOG(severity) LOG_##severity.stream()

// The log message
class LogMessage {
public:
    //Constructors
    LogMessage() : log_stream_(std::cout) {}
    LogMessage(const char* file, int line) : log_stream_(std::cout) {
        log_stream_ << file << ":" << line << ": ";
    }
    LogMessage(const LogMessage&) = delete;
    LogMessage& operator=(const LogMessage&) = delete;
    
    //Another line
    ~LogMessage() { log_stream_ << "\n"; }
    
    std::ostream& stream() { return log_stream_; }
protected:
    std::ostream& log_stream_;
};

class LogMessageFatal {
public:
    LogMessageFatal(const char* file, int line) {
        log_stream_ << file << ":" << line << ": ";
    }
    
    //No copy/assign
    LogMessageFatal(const LogMessageFatal&) = delete;
    LogMessageFatal& operator=(LogMessageFatal&) = delete;
    
    //Die the whole system
    ~LogMessageFatal() {
        LOG_BEFORE_THROW << log_stream_.str();
        throw new std::runtime_error(log_stream_.str());
    }
    
    //The output string stream
    std::ostringstream& stream() { return log_stream_; }
protected:
    std::ostringstream log_stream_;
};
	