# -------------------------------------------
# set_initQbase
# 基本量の初期値を設定
# リスタートするか初期値を代入する
# -------------------------------------------

function set_initQbase(xmax, ymax, zmax, restart_file, init_rho, 
                        init_u, init_v, init_w, init_p,
                        specific_heat_ratio, out_file_front,
                        out_ext, out_dir, restartnum, Rd)
    Qbase=[]
    cellxmax = xmax - 1
    cellymax = ymax - 1
    cellzmax = zmax - 1

    restart_check = 0
    try Qbase = setup_restart_value(cellxmax, cellymax, cellzmax, out_dir, restart_file)
        println("Restart "*restart_file)
        restart_check = 2
    catch 
        restart_check = 1
    end

    if restart_check == 1
        Qbase = setup_init_value(cellxmax, cellymax, cellzmax, init_rho, init_u, init_v, init_w, init_p)
        println("Start Initial condition")
        restart_num = 0
        yplus = zeros(cellxmax, cellymax, cellzmax)
        output_result(0, Qbase, cellxmax, cellymax, cellzmax, out_file_front, out_ext, out_dir, Rd, yplus)
    end

    return Qbase, cellxmax, cellymax, cellzmax, restart_num
end

function setup_init_value(cellxmax, cellymax, cellzmax, init_rho, init_u, init_v, init_w, init_p)
    Qbase = zeros(cellxmax, cellymax, cellzmax, 5)
    for i in 1:cellxmax
        for j in 1:cellymax
            for k in 1:cellzmax
                Qbase[i,j,k,1] = init_rho
                Qbase[i,j,k,2] = init_u
                Qbase[i,j,k,3] = init_v
                Qbase[i,j,k,4] = init_w
                Qbase[i,j,k,5] = init_p
            end
        end
    end
    return Qbase
end

function setup_restart_value(cellxmax, cellymax, cellzmax, out_dir, restart_file)
    Qbase = zeros(cellxmax, cellymax, cellzmax, 5)

    skipnum = 1
    fff = []
    open("result/"*restart_file, "r") do f
        fff = read(f,String)
    end 
    fff = split(fff,"\n",keepempty=false)
    
    for i in 1+skipnum:length(fff)
        fff[i] = replace(fff[i]," \r" => "")
    end
    
    ite = 1
    for i in 2:cellxmax-1
        for j in 2:cellymax-1
            for k in 2:cellzmax-1
                temp = split(fff[ite+skipnum]," ")

                for l in 1:5
                    Qbase[i,j,k,l] = parse(Float64,temp[l])
                end
                ite = ite+1
            end
        end
    end
    return Qbase
end

function quicksort(list, first, last, num_cell_data, standard_num)
    # list = zezros(m,n)
    x = list[Int(floor((first+last)/2)),standard_num]
    i = first
    j = last
    
    while true
        while list[i,standard_num] < x
            i = i+1
        end
        while x < list[j,standard_num]
            j = j-1
        end
        
        if (i >= j) 
            break
        else
            for k in 1:num_cell_data
                t = list[i,k]
                list[i,k] = list[j,k]
                list[j,k] = t
            end
            i = i+1
            j = j-1
        end
    end
    
    if (first < i - 1) 
        list = quicksort(list, first, i - 1,num_cell_data,standard_num)
    end
    if (j + 1 < last) 
        list = quicksort(list, j + 1, last,num_cell_data,standard_num)
    end
    return list
end