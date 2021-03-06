using Printf
function set_initQbase(xmax, ymax, restart_file, init_rho, init_u, init_v, init_p, init_N2, init_N, dt, 
                        specific_heat_ratio, out_file_front, out_ext, out_dir, restart_num, R, nch, nval)
    Qbase=[]
    cellxmax = xmax - 1
    cellymax = ymax - 1

    init_step = 0
    init_time = 0.0

    restart_check = 0
    try Qbase = setup_restart_value(cellxmax, cellymax, out_dir, restart_file, nch, nval)
        println("Restart "*restart_file)
        restart_check = 2
        init_step = restart_num
        init_time = init_step*dt
    catch 
        restart_check = 1
    end
    
    if restart_check == 1
        Qbase = setup_init_value(cellxmax, cellymax, init_rho, init_u, init_v, init_p, init_N2, init_N, nval)
        println("Start Initial condition")
        restart_num = 0
        Rhat = set_gasconst_hat(Qbase, cellxmax, cellymax, nval, nch, R)
        output_result(0, Qbase, cellxmax, cellymax, specific_heat_ratio, out_file_front, out_ext, out_dir, Rhat, nval)
    end

    return Qbase, cellxmax, cellymax, restart_num, init_step, init_time
end

function setup_init_value(cellxmax, cellymax, init_rho, init_u, init_v, init_p, init_N2, init_N, nval)
    Qbase=zeros(cellxmax, cellymax, nval)
    for i in 1:cellxmax
        for j in 1:cellymax
            Qbase[i,j,1] = init_rho
            Qbase[i,j,2] = init_u
            Qbase[i,j,3] = init_v
            Qbase[i,j,4] = init_p
            Qbase[i,j,5] = init_N2
            Qbase[i,j,6] = init_N
        end
    end
    return Qbase
end

function setup_restart_value(cellxmax, cellymax, out_dir, restart_file, nch, nval)
    Qbase = zeros(cellxmax, cellymax, nval)

    skipnum = 1
    fff = []
    open("result/"*restart_file, "r") do f
        fff = read(f, String)
    end 
    fff = split(fff,"\n", keepempty = false)
    
    for i in 1+skipnum:length(fff)
        fff[i] = replace(fff[i]," \r" => "")
    end
    
    ite = 1
    npre = nval - nch
    for i in 2:cellxmax-1
        for j in 2:cellymax-1
            temp = split(fff[ite+skipnum]," ")
            for l in 1:nval
                Qbase[i,j,l] = parse(Float64,temp[l])
            end
            for l in 1:nch # mass frac
                Qbase[i,j,npre+l] = Qbase[i,j,npre+l]
            end
            ite = ite+1
        end
    end
    return Qbase
end

function check_divrege(Qbase, cellxmax, cellymax, Rhat, fwrite)
    ite = 0
    for i in 1:cellxmax
        for j in 1:cellymax
            T = Qbase[i,j,4]/(Qbase[i,j,1]*Rhat[i,j])
            if T < 0 || isequal(T, NaN) == true
                open( fwrite, "a" ) do f
                    ai = @sprintf("%4.0f", i)
                    aj = @sprintf("%4.0f", j)

                    write(f, "\n")
                    write(f, " thermal diverge ")
                    write(f, "\n")
                    write(f, " i = "*ai)
                    write(f, "\n")
                    write(f, " j = "*aj)
                    write(f, "\n")
                end
                ite = 1
            end            
            if 1+1e-3 < Qbase[i,j,5] || Qbase[i,j,5] < -1e-3
                open( fwrite, "a" ) do f
                    ai = @sprintf("%4.0f", i)
                    aj = @sprintf("%4.0f", j)

                    write(f, "\n")
                    write(f, " chemical diverge ")
                    write(f, "\n")
                    write(f, " i = "*ai)
                    write(f, "\n")
                    write(f, " j = "*aj)
                    write(f, "\n")
                end
                ite = 2
            end
        end
    end

    if ite == 1
        println("\n")
        println(" T < 0 or T = NaN ")
        println(" thermal diverge ")
        println("\n")
        throw(UndefVarError(:x))
    elseif ite == 2
        println("\n")
        println(" N > 1 or N < 0 ")
        println(" chemical diverge ")
        println("\n")
        throw(UndefVarError(:x))
    end
end

function check_small(Qbase, cellxmax, cellymax, nch)
    for i in 1:cellxmax
        for j in 1:cellymax
            for ns in 1:nch
                if Qbase[i,j,4+ns] < 1e-20
                    Qbase[i,j,4+ns] = 0.0
                end
            end
        end
    end
end

function GE(matrix, nval)
    # nval * nval 行列の逆行列を求める。

    for l in 1:nval
        a = 1.0/matrix[l,l]
        matrix[l,l] = 1.0

        for m in 1:nval
            matrix[l,m] = matrix[l,m]*a
        end

        for m in 1:nval
            if m == l
            else
                b = matrix[m,l]
                matrix[m,l] = 0.0
                for n in 1:nval
                    matrix[m,n] = matrix[m,n] - b*matrix[l,n]
                end
            end
        end
    end
    return matrix
end