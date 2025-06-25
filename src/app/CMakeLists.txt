add_executable(plot-shape-function plot-shape-function.cpp)
target_link_libraries(plot-shape-function common mesh interpolation)

add_executable(plot-rdm-beta plot-rdm-beta.cpp)
target_link_libraries(plot-rdm-beta common mesh interpolation)

add_executable(compute-Vandermonde-cond-number compute-Vandermonde-cond-number.cpp)
target_link_libraries(compute-Vandermonde-cond-number common math mesh interpolation)

add_executable(Runge-phenomenon Runge-phenomenon.cpp)
target_link_libraries(Runge-phenomenon common math mesh interpolation)
