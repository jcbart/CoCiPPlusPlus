#ifndef COCIPLOG_H
#define COCIPLOG_H

#include <string_view>

// Writes a message to the log
void CoCiP_LogWrite(const std::string_view msg);

// Logs an error, then exits
void CoCiP_RaiseError(const std::string_view msg, const std::string_view filename, const int line);

#endif