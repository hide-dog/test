function set_bdx(Qbase, rho, u, v, p, ch, nch, icell, j)
    Qbase[icell, j, 1] = rho
    Qbase[icell, j, 2] = u
    Qbase[icell, j, 3] = v
    Qbase[icell, j, 4] = p
    for ns in 1:nch
        Qbase[icell, j, 4+ns] = ch[ns]
    end
    return Qbase
end

function set_bdy(Qbase, rho, u, v, p, ch, nch, jcell, i)
    Qbase[i, jcell, 1] = rho
    Qbase[i, jcell, 2] = u
    Qbase[i, jcell, 3] = v
    Qbase[i, jcell, 4] = p
    for ns in 1:nch
        Qbase[i, jcell, 4+ns] = ch[ns]
    end
    return Qbase
end

function set_boundary(Qbase, cellxmax, cellymax, vecAx, vecAy, bdcon, Rhat, nval, nch)
    """bdcon[i][j]
        i : 境界番号(x-,x+ ,y-,y+)
        j=1-6 : "bd1_con":"2",
                "bd1_rho":"1.0",
                "bd1_u"  :"300.0",
                "bd1_v"  :"0.0",
                "bd1_p"  :"1.0",
                "bd1_N2" :"1.0",
                "bd1_N"  :"0.0",
                "bd1_T"  :"300.0",
    """
    rho = 0.0
    u   = 0.0
    v   = 0.0
    p   = 0.0
    ch  = zeros(nch)
    
    # bd1 = x-
    # bd2 = x+
    bdx         = [1, cellxmax]
    bdx_nei     = [2, cellxmax-1]
    bdx_nei_rev = [cellxmax-1, 2]
    for dim in 1:2
        if Int(bdcon[dim][1]) == 00
            for j in 1:cellymax
                rho = bdcon[dim][2]
                u   = bdcon[dim][3]
                v   = bdcon[dim][4]
                p   = bdcon[dim][5]
                for ns in 1:nch
                    ch[ns] =  bdcon[dim][5+ns]
                end
                set_bdx(Qbase, rho, u, v, p, ch, nch, bdx[dim], j)
            end
        elseif Int(bdcon[dim][1]) == 01
            for j in 1:cellymax
                rho = bdcon[dim][2]
                u   = bdcon[dim][3]
                v   = bdcon[dim][4]
                T   = bdcon[dim][6+nch]
                p   = rho*Rhat[i,j]*T
                for ns in 1:nch
                    ch[ns] =  bdcon[dim][5+ns]
                end
                set_bdx(Qbase, rho, u, v, p, ch, nch, bdx[dim], j)
            end
        elseif Int(bdcon[dim][1]) == 20
            for j in 1:cellymax
                rho = Qbase[bdx_nei[dim],j,1]
                u   = Qbase[bdx_nei[dim],j,2]
                v   = Qbase[bdx_nei[dim],j,3]
                p   = Qbase[bdx_nei[dim],j,4]
                for ns in 1:nch
                    ch[ns] =  Qbase[bdx_nei[dim],j,4+ns]
                end
                set_bdx(Qbase, rho, u, v, p, ch, nch, bdx[dim], j)
            end
        elseif Int(bdcon[dim][1]) == 30
            for j in 1:cellymax
                rho = Qbase[bdx_nei[dim],j,1]
                u   = -Qbase[bdx_nei[dim],j,2]
                v   = -Qbase[bdx_nei[dim],j,3]
                p   = Qbase[bdx_nei[dim],j,4]
                for ns in 1:nch
                    ch[ns] =  Qbase[bdx_nei[dim],j,4+ns]
                end
                set_bdx(Qbase, rho, u, v, p, ch, nch, bdx[dim], j)
            end
        elseif Int(bdcon[dim][1]) == 31
            for j in 1:cellymax
                u   = -Qbase[bdx_nei[dim],j,2]
                v   = -Qbase[bdx_nei[dim],j,3]
                p   = Qbase[bdx_nei[dim],j,4]
                T   = bdcon[dim][6+nch]
                rho = p/(Rhat[i,j]*T)
                for ns in 1:nch
                    ch[ns] =  Qbase[bdx_nei[dim],j,4+ns]
                end
                set_bdx(Qbase, rho, u, v, p, ch, nch, bdx[dim], j)
            end
        elseif Int(bdcon[dim][1]) == 40
            for j in 1:cellymax
                rho = Qbase[bdx_nei_rev[dim],j,1]
                u   = Qbase[bdx_nei_rev[dim],j,2]
                v   = Qbase[bdx_nei_rev[dim],j,3]
                p   = Qbase[bdx_nei_rev[dim],j,4]
                for ns in 1:nch
                    ch[ns] =  Qbase[bdx_nei_rev[dim],j,4+ns]
                end
                set_bdx(Qbase, rho, u, v, p, ch, nch, bdx[dim], j)
            end
        end
    end

    

    # bd3 = y-
    # bd4 = y+
    bdy     = [1, cellymax]
    bdy_nei     = [2, cellymax-1]
    bdy_nei_rev = [cellymax-1, 2]
    for dim in 3:4
        if Int(bdcon[dim][1]) == 00
            for i in 1:cellxmax
                rho = bdcon[dim][2]
                u   = bdcon[dim][3]
                v   = bdcon[dim][4]
                p   = bdcon[dim][5]
                for ns in 1:nch
                    ch[ns] =  bdcon[dim][5+ns]
                end
                set_bdy(Qbase, rho, u, v, p, ch, nch, bdy[dim-2],i)
            end
        elseif Int(bdcon[dim][1]) == 01
            for i in 1:cellxmax
                rho = bdcon[dim][2]
                u   = bdcon[dim][3]
                v   = bdcon[dim][4]
                T   = bdcon[dim][6+nch]
                p   = rho*Rhat[i,j]*T
                for ns in 1:nch
                    ch[ns] =  bdcon[dim][5+ns]
                end
                set_bdy(Qbase, rho, u, v, p, ch, nch, bdy[dim-2],i)
            end
        elseif Int(bdcon[dim][1]) == 02
            temp1 = Int64(round(cellxmax/4)+1) # 1/4流出
            temp2 = Int64((temp1-1)*3+1)       # 1/4~3/4流入
            for i in 1:temp1
                rho = Qbase[i,bdy_nei[dim-2],1]
                u   = Qbase[i,bdy_nei[dim-2],2]
                v   = Qbase[i,bdy_nei[dim-2],3]
                p   = Qbase[i,bdy_nei[dim-2],4]
                for ns in 1:nch
                    ch[ns] =  Qbase[i,bdy_nei[dim-2],4+ns]
                end
                set_bdy(Qbase, rho, u, v, p, ch, nch, bdy[dim-2],i)
            end
            for i in temp1+1:temp2
                rho = bdcon[dim][2]
                u   = bdcon[dim][3]
                v   = bdcon[dim][4]
                p   = bdcon[dim][5]
                for ns in 1:nch
                    ch[ns] =  bdcon[dim][5+ns]
                end
                set_bdy(Qbase, rho, u, v, p, ch, nch, bdy[dim-2],i)
            end
            for i in temp2+1:cellxmax
                rho = Qbase[i,bdy_nei[dim-2],1]
                u   = Qbase[i,bdy_nei[dim-2],2]
                v   = Qbase[i,bdy_nei[dim-2],3]
                p   = Qbase[i,bdy_nei[dim-2],4]
                for ns in 1:nch
                    ch[ns] =  Qbase[i,bdy_nei[dim-2],4+ns]
                end
                set_bdy(Qbase, rho, u, v, p, ch, nch, bdy[dim-2],i)
            end
        elseif Int(bdcon[dim][1]) == 03
            temp1 = Int64(round(cellxmax/4)+1) # 1/4流出
            temp2 = Int64((temp1-1)*3+1)       # 1/4~3/4流入
            for i in 1:temp1
                rho = Qbase[i,bdy_nei[dim-2],1]
                u   = Qbase[i,bdy_nei[dim-2],2]
                v   = Qbase[i,bdy_nei[dim-2],3]
                p   = Qbase[i,bdy_nei[dim-2],4]
                for ns in 1:nch
                    ch[ns] =  Qbase[i,bdy_nei[dim-2],4+ns]
                end
                set_bdy(Qbase, rho, u, v, p, ch, nch, bdy[dim-2],i)
            end
            for i in temp1+1:temp2
                rho = bdcon[dim][2]
                u   = bdcon[dim][3]
                v   = bdcon[dim][4]
                T   = bdcon[dim][6+nch]
                p   = rho*Rhat[i,j]*T
                for ns in 1:nch
                    ch[ns] =  bdcon[dim][5+ns]
                end
                set_bdy(Qbase, rho, u, v, p, ch, nch, bdy[dim-2],i)
            end
            for i in temp2+1:cellxmax
                rho = Qbase[i,bdy_nei[dim-2],1]
                u   = Qbase[i,bdy_nei[dim-2],2]
                v   = Qbase[i,bdy_nei[dim-2],3]
                p   = Qbase[i,bdy_nei[dim-2],4]
                for ns in 1:nch
                    ch[ns] =  Qbase[i,bdy_nei[dim-2],4+ns]
                end
                set_bdy(Qbase, rho, u, v, p, ch, nch, bdy[dim-2],i)
            end
        elseif Int(bdcon[dim][1]) == 20
            for i in 1:cellxmax
                rho = Qbase[i,bdy_nei[dim-2],1]
                u   = Qbase[i,bdy_nei[dim-2],2]
                v   = Qbase[i,bdy_nei[dim-2],3]
                p   = Qbase[i,bdy_nei[dim-2],4]
                for ns in 1:nch
                    ch[ns] =  Qbase[i,bdy_nei[dim-2],4+ns]
                end
                set_bdy(Qbase, rho, u, v, p, ch, nch, bdy[dim-2],i)
            end
        elseif Int(bdcon[dim][1]) == 30
            for i in 1:cellxmax
                rho = Qbase[i,bdy_nei[dim-2],1]
                u   = -Qbase[i,bdy_nei[dim-2],2]
                v   = -Qbase[i,bdy_nei[dim-2],3]
                p   = Qbase[i,bdy_nei[dim-2],4]
                for ns in 1:nch
                    ch[ns] =  Qbase[i,bdy_nei[dim-2],4+ns]
                end
                set_bdy(Qbase, rho, u, v, p, ch, nch, bdy[dim-2],i)
            end
        elseif Int(bdcon[dim][1]) == 31
            for i in 1:cellxmax
                u   = -Qbase[i,bdy_nei[dim-2],2]
                v   = -Qbase[i,bdy_nei[dim-2],3]
                p   = Qbase[i,bdy_nei[dim-2],4]
                T   = bdcon[dim][6+nch]
                rho = p/(Rhat[i,j]*T)
                for ns in 1:nch
                    ch[ns] =  Qbase[i,bdy_nei[dim-2],4+ns]
                end
                set_bdy(Qbase, rho, u, v, p, ch, nch, bdy[dim-2],i)
            end
        elseif Int(bdcon[dim][1]) == 40
            for i in 1:cellxmax
                rho = Qbase[i,bdy_nei_rev[dim-2],1]
                u   = Qbase[i,bdy_nei_rev[dim-2],2]
                v   = Qbase[i,bdy_nei_rev[dim-2],3]
                p   = Qbase[i,bdy_nei_rev[dim-2],4]
                for ns in 1:nch
                    ch[ns] = Qbase[i,bdy_nei_rev[dim-2],4+ns]
                end
                set_bdy(Qbase, rho, u, v, p, ch, nch, bdy[dim-2],i)
            end
        end
    end

    return Qbase
end

function check_bd(bdcon)
    bd = -1
    for l in 1:4
        if     Int(bdcon[l][1]) == 00
            bd = 00
        elseif Int(bdcon[l][1]) == 01
            bd = 01
        elseif Int(bdcon[l][1]) == 02
            bd = 02
        elseif Int(bdcon[l][1]) == 03
            bd = 03

        elseif Int(bdcon[l][1]) == 20
            bd = 20

        elseif Int(bdcon[l][1]) == 30
            bd = 30
        elseif Int(bdcon[l][1]) == 31
            bd = 31
            
        elseif Int(bdcon[l][1]) == 40
            bd = 40
        else
            println("\n check boundary condition ! \n")
            throw(UndefVarError(:x))
        end
        println("bd condition")
        println(bd)
    end
end