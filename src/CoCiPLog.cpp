#include <iostream>
#include <string_view>
#include <format>
#include <CoCiP++/CoCiPLog.h>

void CoCiP_LogWrite(const std::string_view msg) {
    std::cout << msg << std::endl;
}

void CoCiP_RaiseError(const std::string_view msg, std::string_view filename, const int line) {
    std::string error = std::format("CoCiP error in {}, line {}: {}", filename, line, msg);
    std::cerr << error << std::endl;
    exit(EXIT_FAILURE);
}