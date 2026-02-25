#include <iostream>
#include <string>
#include <CoCiP++/CoCiPLog.h>

void CoCiP_LogWrite(const char* msg) {
    std::cout << msg << std::endl;
}

void CoCiP_LogWrite(std::string msg) {
    CoCiP_LogWrite(msg.c_str());
}

void CoCiP_RaiseError(const char* msg, const char* filename, const int line) {
    std::string error = std::string("CoCiP error in ") + filename + std::string(", line ")
                        + std::to_string(line) + std::string(": ") + msg;

    std::cerr << error << std::endl;
    exit(EXIT_FAILURE);
}

void CoCiP_RaiseError(std::string msg, const char* filename, const int line) {
    CoCiP_RaiseError(msg.c_str(), filename, line);
}