# ----------------------
# -- read             --
# ----------------------
function read_nodenum(skipnum)
    """ 
    xmax : 仮想セルも含めたnodeのxの数
    ymax : 仮想セルも含めたnodeのyの数
    zmax : 仮想セルも含めたnodeのzの数
    """
    fff=[]
    open("grid/nodesnum", "r") do f
        fff=read(f,String)
    end 
    fff = split(fff,"\n",keepempty=false)
    num_nodes = length(fff)-skipnum

    for i in 1+skipnum:length(fff)
        fff[i] = replace(fff[i]," \r" => "")
    end

    temp = split(fff[2]," ")
    xmax = parse(Int64,temp[1]) 
    ymax = parse(Int64,temp[2])
    zmax = parse(Int64,temp[3])
    
    return xmax, ymax, zmax
end

function read_nodes(skipnum,xmax,ymax,zmax)
    """ 
    nodes[i][j][k][l]
    i : x点の番号
    j : y点の番号
    k : z点の番号
    l=1 : 点のx座標
    l=2 : 点のy座標
    l=3 : 点のz座標
    """

    fff=[]
    open("grid/nodes", "r") do f
        fff = read(f,String)
    end 
    fff = split(fff,"\n",keepempty=false)
    num_nodes = length(fff)-skipnum
    
    for i in 1+skipnum:length(fff)
        fff[i] = replace(fff[i]," \r" => "")
    end

    nodes = zeros(xmax, ymax, zmax, 3)
    for i in 1:num_nodes
        temp = split(fff[i+skipnum]," ")

        xnum = parse(Int64,temp[1])
        ynum = parse(Int64,temp[2])
        znum = parse(Int64,temp[3])

        nodes[xnum,ynum,znum,1] = parse(Float64,temp[4])
        nodes[xnum,ynum,znum,2] = parse(Float64,temp[5]) 
        nodes[xnum,ynum,znum,3] = parse(Float64,temp[6])
    end

    return nodes
end 

function read_vecA(skipnum,xmax,ymax,zmax)
    """ vecAx[i][j][k][l]
        i : x番号
        j : y番号
        k : z番号
        l=1 : x vec
        l=2 : y vec
        l=3 : z vec
    """
    fff=[]
    open("grid/vecAx", "r") do f
        fff = read(f,String)
    end 
    fff = split(fff,"\n",keepempty=false)
    num_cell = length(fff)-skipnum
    
    for i in 1+skipnum:length(fff)
        fff[i] = replace(fff[i]," \r" => "")
    end

    vecAx = zeros(xmax, ymax-1, zmax-1, 3)
    for i in 1:num_cell
        temp = split(fff[i+skipnum]," ")

        x = parse(Int64,temp[1])
        y = parse(Int64,temp[2])
        z = parse(Int64,temp[3])
        vecAx[x,y,z,1] = parse(Float64,temp[4])
        vecAx[x,y,z,2] = parse(Float64,temp[5])
        vecAx[x,y,z,3] = parse(Float64,temp[6])
    end



    fff=[]
    open("grid/vecAy", "r") do f
        fff = read(f,String)
    end 
    fff = split(fff,"\n",keepempty=false)
    num_cell = length(fff)-skipnum
    
    for i in 1+skipnum:length(fff)
        fff[i] = replace(fff[i]," \r" => "")
    end

    vecAy = zeros(xmax-1, ymax, zmax-1, 3)
    for i in 1:num_cell
        temp = split(fff[i+skipnum]," ")

        x = parse(Int64,temp[1])
        y = parse(Int64,temp[2])
        z = parse(Int64,temp[3])
        vecAy[x,y,z,1] = parse(Float64,temp[4])
        vecAy[x,y,z,2] = parse(Float64,temp[5]) 
        vecAy[x,y,z,3] = parse(Float64,temp[6]) 
    end



    fff=[]
    open("grid/vecAz", "r") do f
        fff = read(f,String)
    end 
    fff = split(fff,"\n",keepempty=false)
    num_cell = length(fff)-skipnum
    
    for i in 1+skipnum:length(fff)
        fff[i] = replace(fff[i]," \r" => "")
    end

    vecAz = zeros(xmax-1, ymax-1, zmax, 3)
    for i in 1:num_cell
        temp = split(fff[i+skipnum]," ")

        x = parse(Int64,temp[1])
        y = parse(Int64,temp[2])
        z = parse(Int64,temp[3])
        vecAz[x,y,z,1] = parse(Float64,temp[4])
        vecAz[x,y,z,2] = parse(Float64,temp[5]) 
        vecAz[x,y,z,3] = parse(Float64,temp[6]) 
    end
    return vecAx, vecAy, vecAz
end 

function read_allgrid()
    skip=1
    xmax, ymax, zmax    = read_nodenum(skip)
    nodes               = read_nodes(skip, xmax, ymax, zmax)
    vecAx, vecAy, vecAz = read_vecA(skip, xmax, ymax, zmax)
    println("fin read grid")
    return  xmax, ymax, zmax, nodes, vecAx, vecAy, vecAz
end
