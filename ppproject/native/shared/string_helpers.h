#pragma once
#include <cstdarg>
#include <cstdlib>
#include <string>

std::string StringFromFormat(const char* fmt, ...);
std::string StringFromFormatV(const char* fmt, va_list ap);
std::string HexDumpString(const void* buf, size_t len);
