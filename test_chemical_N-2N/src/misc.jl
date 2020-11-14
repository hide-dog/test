function set_initQbase(xmax, ymax, restart_file, init_rho, init_u, init_v, init_p, init_N2, init_N, 
                        specific_heat_ratio, out_file_front, out_ext, out_dir, restart_num, R, nch, nval)
    Qbase=[]
    cellxmax = xmax - 1
    cellymax = ymax - 1

    restart_check = 0
    try Qbase = setup_restart_value(cellxmax, cellymax, out_dir, restart_file, nch, nval)
        println("Restart "*restart_file)
        restart_check = 2
    catch 
        restart_check = 1
    end

    if restart_check == 1
        Qbase = setup_init_value(cellxmax, cellymax, init_rho, init_u, init_v, init_p, init_N2, init_N, nval)
        println("Start Initial condition")
        restart_num = 0
        Rhat = set_gasconst(Qbase, cellxmax, cellymax, nval, nch, R)
        output_result(0, Qbase, cellxmax, cellymax, specific_heat_ratio, out_file_front, out_ext, out_dir, Rhat, nval)
    end

    return Qbase, cellxmax, cellymax, restart_num
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
            for l in 1:nch # mass frac -> rho
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
            if 1 < Qbase[i,j,5] || Qbase[i,j,5] < 0
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