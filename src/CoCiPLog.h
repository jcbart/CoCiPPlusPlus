#ifndef COCIPLOG_H
#define COCIPLOG_H

#include <string>

// Writes a message to the log
void CoCiP_LogWrite(const char* msg);

// Writes a message to the log
void CoCiP_LogWrite(std::string msg);

// Logs an error, then exits
void CoCiP_RaiseError(const char* msg, const char* filename, const int line);

// Logs an error, then exits
void CoCiP_RaiseError(std::string msg, const char* filename, const int line);

#endif