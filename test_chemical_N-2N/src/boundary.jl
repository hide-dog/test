function set_boundary(Qbase, cellxmax, cellymax, vecAx, vecAy, bdcon, Rhat, nval)
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

    # bd1 = x-
    if Int(bdcon[1][1]) == 0
        for j in 1:cellymax
            for l in 1:nval
                Qbase[1,j,l] = bdcon[1][l+1]
            end
        end
    elseif Int(bdcon[1][1]) == 1
        for j in 1:cellymax
            for l in 1:nval
                Qbase[1,j,l] = Qbase[2,j,l]
            end
        end
    elseif Int(bdcon[1][1]) == 2
        for j in 1:cellymax
            for l in 1:nval
                Qbase[1,j,l] = Qbase[2,j,l]
            end

            xvec = vecAx[2,j,1]
            yvec = vecAx[2,j,2]
            u = Qbase[2,j,2]
            v = Qbase[2,j,3]

            Qbase[1,j,2] = ((-xvec^2+yvec^2)*u-2*xvec*yvec*v)/(xvec^2+yvec^2)
            Qbase[1,j,3] = (-2*xvec*yvec*u+(xvec^2-yvec^2)*v)/(xvec^2+yvec^2)            
        end
    elseif Int(bdcon[1][1]) == 3
        for j in 1:cellymax
            for l in 1:nval
                Qbase[1,j,l] = Qbase[2,j,l]
            end

            u = Qbase[2,j,2]
            v = Qbase[2,j,3]

            Qbase[1,j,2] = -u
            Qbase[1,j,3] = -v
        end
    elseif Int(bdcon[1][1]) == 4
        for j in 1:cellymax
            for l in 1:nval
                Qbase[1,j,l] = Qbase[cellxmax-1,j,l]
            end
        end
    end

    # bd2 = x+
    if Int(bdcon[2][1]) == 0
        for j in 1:cellymax
            for l in 1:nval
                Qbase[cellxmax,j,l] = bdcon[2][l+1]
            end
        end
    elseif Int(bdcon[2][1]) == 1
        for j in 1:cellymax
            for l in 1:nval
                Qbase[cellxmax,j,l] = Qbase[cellxmax-1,j,l]
            end
        end
    elseif Int(bdcon[2][1]) == 2
        for j in 1:cellymax
            for l in 1:nval
                Qbase[cellxmax,j,l] = Qbase[cellxmax-1,j,l]
            end

            xvec = vecAx[cellymax,j,1]
            yvec = vecAx[cellymax,j,2]
            u = Qbase[cellxmax-1,j,2]
            v = Qbase[cellxmax-1,j,3]

            Qbase[cellxmax,j,2] = ((-xvec^2+yvec^2)*u-2*xvec*yvec*v)/(xvec^2+yvec^2)
            Qbase[cellxmax,j,3] = (-2*xvec*yvec*u+(xvec^2-yvec^2)*v)/(xvec^2+yvec^2)
        end
    elseif Int(bdcon[2][1]) == 3
        for j in 1:cellymax
            for l in 1:nval
                Qbase[cellxmax,j,l] = Qbase[cellxmax-1,j,l]
            end
            u = Qbase[cellxmax-1,j,2]
            v = Qbase[cellxmax-1,j,3]

            Qbase[cellxmax,j,2] = -u
            Qbase[cellxmax,j,3] = -v
        end
    elseif Int(bdcon[2][1]) == 4
        for j in 1:cellymax
            for l in 1:nval
                Qbase[cellxmax,j,l] = Qbase[2,j,l]
            end
        end
    end

    # bd3 = y-
    if Int(bdcon[3][1]) == 0
        for i in 1:cellxmax
            for l in 1:nval
                Qbase[i,1,l] = bdcon[3][l+1]
            end
        end
    elseif Int(bdcon[3][1]) == 1
        for i in 1:cellxmax
            for l in 1:nval
                Qbase[i,1,l] = Qbase[i,2,l]
            end
        end
    elseif Int(bdcon[3][1]) == 2
        for i in 1:cellxmax
            for l in 1:nval
                Qbase[i,1,l] = Qbase[i,2,l]
            end

            xvec = vecAy[i,2,1]
            yvec = vecAy[i,2,2]
            u = Qbase[i,2,2]
            v = Qbase[i,2,3]
            
            Qbase[i,1,2] = ((-xvec^2+yvec^2)*u-2*xvec*yvec*v)/(xvec^2+yvec^2)
            Qbase[i,1,3] = (-2*xvec*yvec*u+(xvec^2-yvec^2)*v)/(xvec^2+yvec^2)
        end
    elseif Int(bdcon[3][1]) == 3
        for i in 1:cellxmax
            for l in 1:nval
                Qbase[i,1,l] = Qbase[i,2,l]
            end

            u = Qbase[i,2,2]
            v = Qbase[i,2,3]

            Qbase[i,1,2] = -u
            Qbase[i,1,3] = -v
        end
    elseif Int(bdcon[3][1]) == 4
        for i in 1:cellxmax
            for l in 1:nval
                Qbase[i,1,l] = Qbase[i,cellymax-1,l]
            end 
        end
    end

    # bd4 = y+
    if Int(bdcon[4][1]) == 0
        for i in 1:cellxmax
            for l in 1:nval
                Qbase[i,cellymax,l] = bdcon[4][l+1]
            end
        end
    elseif Int(bdcon[4][1]) == 1
        for i in 1:cellxmax
            for l in 1:nval
                Qbase[i,cellymax,l] = Qbase[i,cellymax-1,l]
            end
        end
    elseif Int(bdcon[4][1]) == 2
        for i in 1:cellxmax
            for l in 1:nval
                Qbase[i,cellymax,l] = Qbase[i,cellymax-1,l]
            end

            xvec = vecAy[i,cellymax,1]
            yvec = vecAy[i,cellymax,2]
            u = Qbase[i,cellymax-1,2]
            v = Qbase[i,cellymax-1,3]

            Qbase[i,cellymax,2] = ((-xvec^2+yvec^2)*u-2*xvec*yvec*v)/(xvec^2+yvec^2)
            Qbase[i,cellymax,3] = (-2*xvec*yvec*u+(xvec^2-yvec^2)*v)/(xvec^2+yvec^2)
        end
    elseif Int(bdcon[4][1]) == 3
        for i in 1:cellxmax
            for l in 1:nval
                Qbase[i,cellymax,l] = Qbase[i,cellymax-1,l]
            end

            u = Qbase[i,cellymax-1,2]
            v = Qbase[i,cellymax-1,3]

            Qbase[i,cellymax,2] = -u
            Qbase[i,cellymax,3] = -v
        end
    elseif Int(bdcon[4][1]) == 4
        for i in 1:cellxmax
            for l in 1:nval
                Qbase[i,cellymax,l] = Qbase[i,2,l]
            end 
        end
    elseif Int(bdcon[4][1]) == 5 
        temp1 = Int64(round(cellxmax/4)+1) # 1/4流出
        temp2 = Int64((temp1-1)*3+1)       # 1/4~3/4流入
        for i in 1:temp1
            for l in 1:nval
                Qbase[i,cellymax,l] = Qbase[i,cellymax-1,l]
            end
        end
        for i in temp1+1:temp2
            for l in 1:nval
                Qbase[i,cellymax,l] = bdcon[4][l+1]
            end
        end
        for i in temp2+1:cellxmax
            for l in 1:nval
                Qbase[i,cellymax,l] = Qbase[i,cellymax-1,l]
            end
        end
    elseif Int(bdcon[4][1]) == 6
        temp1 = Int64(round(cellxmax/4)+1) # 1/4流出
        temp2 = Int64((temp1-1)*3+1)       # 1/4~3/4流入
        for i in 1:temp1
            for l in 1:nval
                Qbase[i,cellymax,l] = Qbase[i,cellymax-1,l]
            end
        end
        for i in temp1+1:temp2
            for l in 1:nval
                Qbase[i,cellymax,l] = bdcon[4][l+1]
            end
            T = bdcon[4][nval+2]
            p = (Qbase[i,cellymax,1]*Rhat[i,cellymax]) * T
            Qbase[i,cellymax,4] = p
        end
        for i in temp2+1:cellxmax
            for l in 1:nval
                Qbase[i,cellymax,l] = Qbase[i,cellymax-1,l]
            end
        end
    end

    return Qbase
end

function check_bd(bdcon)
    for l in 1:4
        if Int(bdcon[l][1]) == 0
        elseif Int(bdcon[l][1]) == 1
        elseif Int(bdcon[l][1]) == 2
        elseif Int(bdcon[l][1]) == 3
        elseif Int(bdcon[l][1]) == 4
        elseif Int(bdcon[l][1]) == 5
        elseif Int(bdcon[l][1]) == 6
        else
            println("\n check boundary condition ! \n")
            throw(UndefVarError(:x))
        end
    end
end