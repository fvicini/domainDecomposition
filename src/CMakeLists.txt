list(APPEND domain_decomposition_headers ${CMAKE_CURRENT_SOURCE_DIR}/EllipticProblem.hpp)

list(APPEND domain_decomposition_sources ${CMAKE_CURRENT_SOURCE_DIR}/EllipticProblem.cpp)

list(APPEND domain_decomposition_includes ${CMAKE_CURRENT_SOURCE_DIR})

set(domain_decomposition_sources ${domain_decomposition_sources} PARENT_SCOPE)
set(domain_decomposition_headers ${domain_decomposition_headers} PARENT_SCOPE)
set(domain_decomposition_includes ${domain_decomposition_includes} PARENT_SCOPE)
