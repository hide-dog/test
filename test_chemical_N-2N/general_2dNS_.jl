# src path
#src_path = "C:\\Users\\hidee\\OneDrive\\デスクトップ\\git\\test_chemical_N-2N\\2d_general\\src\\"
src_path = "C:\\Users\\Koike\\Desktop\\git\\test_chemical_N-2N\\2d_general\\src\\"

# main (変更しないこと)
src_read="read_grid.jl"
include(src_path*src_read)
src_read="read_para.jl"
include(src_path*src_read)
src_read="output.jl"
include(src_path*src_read)

src_read="advection_term.jl"
include(src_path*src_read)
src_read="boundary.jl"
include(src_path*src_read)
src_read="converge.jl"
include(src_path*src_read)
src_read="misc.jl"
include(src_path*src_read)
src_read="rhs.jl"
include(src_path*src_read)
src_read="value_setup.jl"
include(src_path*src_read)
src_read="viscos.jl"
include(src_path*src_read)
src_read="chemical_term.jl"
include(src_path*src_read)
src_read="chemical_constant.jl"

include(src_path*src_read)
src_read="setup.jl"
include(src_path*src_read)
src_read="cal_time_step.jl"
include(src_path*src_read)
src_main="main.jl"
include(src_path*src_main)


