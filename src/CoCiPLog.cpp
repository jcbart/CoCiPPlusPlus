#include <iostream>
#include <string>
#ifdef ESMF
#include "ESMC.h"
#endif
#include "CoCiPLog.h"

void CoCiP_LogWrite(const char* msg) {
#ifdef ESMF
    int rc = ESMC_LogWrite(msg, ESMC_LOGMSG_INFO);
#else
    std::cout << msg << std::endl;
#endif
}

void CoCiP_LogWrite(std::string msg) {
    CoCiP_LogWrite(msg.c_str());
}

void CoCiP_RaiseError(const char* msg, const char* filename, const int line) {
    std::string error = std::string("Error in ") + filename + std::string(", line ")
                        + std::to_string(line) + std::string(": ") + msg;

#ifdef ESMF
    int rc = ESMC_LogWrite(error.c_str(), ESMC_LOGMSG_ERROR);
#endif
    std::cerr << error << std::endl;
    exit(EXIT_FAILURE);
}

void CoCiP_RaiseError(std::string msg, const char* filename, const int line) {
    CoCiP_RaiseError(msg.c_str(), filename, line);
}