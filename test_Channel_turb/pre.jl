using Printf

function main()
    xy_or_r = 1          # 1:x,y 2:r,t
    xnum = 20           #       2:rnum
    ynum = 20            #       2:tnum
    znum = 20            # z方向に伸ばすだけ（準2次元）    
    lenx = 0.1          #       2:inr
    leny = 0.1          #       2:outr
    lenz = 0.1          #       2:outr
    st_ang  = -1/2*pi       # 2:angle
    en_ang  = -3/2*pi          # 2:angle

    outdir = "grid"
    make_dir(outdir)
    result = "result"
    make_dir(result)
    result = "post_result"
    make_dir(result)
    nodes, xnum_max, ynum_max, znum_max= mk_gird(xnum,ynum,znum,lenx,leny,lenz,xy_or_r,st_ang,en_ang,outdir)
    vecA(nodes,xnum_max,ynum_max,znum_max,outdir)
end

function mk_gird(xnum,ynum,znum,lenx,leny,lenz,xy_or_r,st_ang,en_ang,outdir)
    """
    nodes[i,j,k]
    i   : x,r方向の番号
    j   : y,theta方向の番号
    k=1 : 点のx座標
    k=2 : 点のy座標

    nodes[1,:]      : x方向境界
    nodes[:,1]      : y方向境界
    nodes[xnum+2,:] : x方向境界
    nodes[:,ynum+2] : y方向境界
    """
    xnum_max = xnum+1+2
    ynum_max = ynum+1+2
    znum_max = znum+1+2
    nodes = zeros(xnum_max, ynum_max, znum_max,3)
    if xy_or_r == 1
        for i in 2:xnum_max-1
            for j in 2:ynum_max-1
                for k in 2:znum_max-1
                    x = lenx/(xnum)*(i-2)
                    y = leny/(ynum)*(j-2)
                    z = lenz/(znum)*(k-2)
                    nodes[i,j,k,1] = x
                    nodes[i,j,k,2] = y
                    nodes[i,j,k,3] = z
                end
            end 
        end
    elseif xy_or_r == 2
        for i in 2:xnum_max-1
            for j in 2:ynum_max-1
                for k in 2:znum_max-1
                    r = lenx + (leny-lenx)/ynum*(j-2)
                    theta = st_ang - (st_ang-en_ang)/xnum*(i-2)
                    x = r * cos(theta)
                    y = r * sin(theta)
                    z = lenz/(znum)*(k-2)
                    nodes[i,j,k,1]=x
                    nodes[i,j,k,2]=y
                    nodes[i,j,k,3] = z
                end
            end 
        end
    end

    #=
    仮想セルの作成：現在あるセル境界線を延長して作成
    ベクトルで考えれば下記のような計算になる．たぶん
    また、[1,1]等の角にはダミー値が入っている.
    =#
    for j in 1:ynum_max
        for k in 1:znum_max
            for l in 1:3
                nodes[1,j,k,l]        =  2.0*nodes[2,j,k,l] - nodes[3,j,k,l]
                nodes[xnum_max,j,k,l] =  2.0*nodes[xnum_max-1,j,k,l] - nodes[xnum_max-2,j,k,l]
            end
        end
    end
    for k in 1:znum_max
        for i in 1:xnum_max
            for l in 1:3
                nodes[i,1,k,l]        =  2.0*nodes[i,2,k,l] - nodes[i,3,k,l]
                nodes[i,ynum_max,k,l] =  2.0*nodes[i,ynum_max-1,k,l] - nodes[i,ynum_max-2,k,l]
            end
        end
    end
    for i in 1:xnum_max
        for j in 1:ynum_max
            for l in 1:3
                nodes[i,j,1,l]        =  2.0*nodes[i,j,2,l] - nodes[i,j,3,l]
                nodes[i,j,znum_max,l] =  2.0*nodes[i,j,znum_max-1,l] - nodes[i,j,znum_max-2,l]
            end
        end
    end

    fff=outdir*"/nodes"
    open(fff,"w") do f
        write(f,"nodes: xnum, ynum , x, y\n")
        for i in 1:xnum_max
            for j in 1:ynum_max
                for k in 1:znum_max
                    x = @sprintf("%8.8e", nodes[i,j,k,1])
                    y = @sprintf("%8.8e", nodes[i,j,k,2])
                    z = @sprintf("%8.8e", nodes[i,j,k,3])
                    write(f,string(i)*" "*string(j)*" "*string(k)*" "*x*" "*y*" "*z*"\n")
                end
            end
        end
    end
    println("write "*fff)

    # nodes_forvtk
    nodes_num = zeros(Int,xnum_max, ynum_max, znum_max)
    fff=outdir*"/nodes_forvtk"
    open(fff,"w") do f
        write(f,"nodes: xnum, ynum , x, y\n")
        l=1
        for i in 2:xnum_max-1
            for j in 2:ynum_max-1
                for k in 2:znum_max-1
                    x = @sprintf("%8.8e", nodes[i,j,k,1])
                    y = @sprintf("%8.8e", nodes[i,j,k,2])
                    z = @sprintf("%8.8e", nodes[i,j,k,3])
                    write(f,string(l)*" "*x*" "*y*" "*z*"\n")
                    nodes_num[i,j,k] = l
                    l = l+1
                end
            end
        end
    end
    println("write "*fff)

    fff=outdir*"/nodesnum"
    open(fff,"w") do f
        write(f,"nodesnum: xnum_max, ynum_max\n")
        write(f,string(xnum_max)*" "*string(ynum_max)*" "*string(znum_max)*"\n")
    end

    ### element ###
    fff=outdir*"/element_forvtk"
    open(fff,"w") do f
        write(f,"elements:cell_xnum, lup,rup,ldown,rdown \n")
        l=1
        for i in 2:xnum_max-2
            for j in 2:ynum_max-2
                for k in 2:znum_max-2
                    d1 = @sprintf("%1.0f", nodes_num[i,j,k])
                    d2 = @sprintf("%1.0f", nodes_num[i,j+1,k])
                    d3 = @sprintf("%1.0f", nodes_num[i+1,j+1,k])
                    d4 = @sprintf("%1.0f", nodes_num[i+1,j,k])
                    d5 = @sprintf("%1.0f", nodes_num[i,j,k+1])
                    d6 = @sprintf("%1.0f", nodes_num[i,j+1,k+1])
                    d7 = @sprintf("%1.0f", nodes_num[i+1,j+1,k+1])
                    d8 = @sprintf("%1.0f", nodes_num[i+1,j,k+1])
                    write(f,string(l)*" "*d1*" "*d2*" "*d3*" "*d4*" "*d5*" "*d6*" "*d7*" "*d8*"\n")
                    l = l+1
                end
            end
        end
    end
    println("write "*fff)

    return  nodes,xnum_max,ynum_max,znum_max
end

function vecA(nodes,xnum_max,ynum_max,znum_max,outdir)
    a = zeros(3)
    b = zeros(3)

    vecAx=zeros(xnum_max, ynum_max-1, znum_max-1, 3)
    for i in 1:xnum_max
        for j in 1:ynum_max-1
            for k in 1:znum_max-1
                # 3dim Ax=(a2b3-a3b2,a3b1-a1b3,a1b2-a2b1)
                # vec_a = (a1,a2,a3) = nodes[i,j+1,k+1] -nodes[i,j,k] 
                # vec_b = (b1,b2,b3) = nodes[i,j,k+1] -nodes[i,j+1,k] 
                for l in 1:3
                    a[l] = nodes[i,j+1,k+1,l]-nodes[i,j,k,l]
                    b[l] = nodes[i,j,k+1,l]-nodes[i,j+1,k,l]
                end
                vecAx[i,j,k,1] = 0.5*(a[2]*b[3]-a[3]*b[2])
                vecAx[i,j,k,2] = 0.5*(a[3]*b[1]-a[1]*b[3])
                vecAx[i,j,k,3] = 0.5*(a[1]*b[2]-a[2]*b[1])
            end
        end
    end

    vecAy=zeros(xnum_max-1, ynum_max, znum_max-1, 3)
    for i in 1:xnum_max-1
        for j in 1:ynum_max
            for k in 1:znum_max-1
                # 3dim Ay=(a2b3-a3b2,a3b1-a1b3,a1b2-a2b1)
                # vec_a = (a1,a2,a3) = nodes[i+1,j,k+1] -nodes[i,j,k] 
                # vec_b = (b1,b2,b3) = nodes[i+1,j,k] -nodes[i,j,k+1] 
                for l in 1:3
                    a[l] = nodes[i+1,j,k+1,l]-nodes[i,j,k,l]
                    b[l] = nodes[i+1,j,k,l]-nodes[i,j,k+1,l]
                end
                vecAy[i,j,k,1] = 0.5*(a[2]*b[3]-a[3]*b[2])
                vecAy[i,j,k,2] = 0.5*(a[3]*b[1]-a[1]*b[3])
                vecAy[i,j,k,3] = 0.5*(a[1]*b[2]-a[2]*b[1])
            end
        end
    end

    vecAz=zeros(xnum_max-1, ynum_max-1, znum_max, 3)
    for i in 1:xnum_max-1
        for j in 1:ynum_max-1
            for k in 1:znum_max
                # 3dim Az=(a2b3-a3b2,a3b1-a1b3,a1b2-a2b1)
                # vec_a = (a1,a2,a3) = nodes[i+1,j+1,k] -nodes[i,j,k] 
                # vec_b = (b1,b2,b3) = nodes[i,j+1,k] -nodes[i+1,j,k] 
                for l in 1:3
                    a[l] = nodes[i+1,j+1,k,l]-nodes[i,j,k,l]
                    b[l] = nodes[i,j+1,k,l]-nodes[i+1,j,k,l]
                end
                vecAz[i,j,k,1] = 0.5*(a[2]*b[3]-a[3]*b[2])
                vecAz[i,j,k,2] = 0.5*(a[3]*b[1]-a[1]*b[3])
                vecAz[i,j,k,3] = 0.5*(a[1]*b[2]-a[2]*b[1])
            end
        end
    end

    fff=outdir*"/vecAx"
    open(fff,"w") do f
        write(f,"vecAx: xnum, ynum, znum, x vec, y vec, z vec\n")
        for i in 1:xnum_max
            for j in 1:ynum_max-1
                for k in 1:znum_max-1
                    x = @sprintf("%8.8e", vecAx[i,j,k,1])
                    y = @sprintf("%8.8e", vecAx[i,j,k,2])
                    z = @sprintf("%8.8e", vecAx[i,j,k,3])
                    write(f,string(i)*" "*string(j)*" "*string(k)*" "*x*" "*y*" "*z*"\n")
                end
            end
        end
    end
    println("write "*fff)

    fff=outdir*"/vecAy"
    open(fff,"w") do f
        write(f,"vecAy: xnum, ynum, znum, x vec, y vec, z vec\n")
        for i in 1:xnum_max-1
            for j in 1:ynum_max
                for k in 1:znum_max-1
                    x = @sprintf("%8.8e", vecAy[i,j,k,1])
                    y = @sprintf("%8.8e", vecAy[i,j,k,2])
                    z = @sprintf("%8.8e", vecAy[i,j,k,3])
                    write(f,string(i)*" "*string(j)*" "*string(k)*" "*x*" "*y*" "*z*"\n")
                end
            end
        end
    end
    println("write "*fff)

    fff=outdir*"/vecAz"
    open(fff,"w") do f
        write(f,"vecAz: xnum, ynum, znum, x vec, y vec, z vec\n")
        for i in 1:xnum_max-1
            for j in 1:ynum_max-1
                for k in 1:znum_max
                    x = @sprintf("%8.8e", vecAz[i,j,k,1])
                    y = @sprintf("%8.8e", vecAz[i,j,k,2])
                    z = @sprintf("%8.8e", vecAz[i,j,k,3])
                    write(f,string(i)*" "*string(j)*" "*string(k)*" "*x*" "*y*" "*z*"\n")
                end
            end
        end
    end
    println("write "*fff)
end

function make_dir(outdir)
    k = 0
    try rm(outdir,recursive=true)
    catch
        mkdir(outdir)
        k = 1
    end

    if k == 0
        mkdir(outdir)
    end
end

# ---------------------------------
main() 