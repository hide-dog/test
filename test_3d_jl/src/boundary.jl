function set_boundary(Qbase, cellxmax, cellymax, cellzmax, vecAx, vecAy, vecAz, bdcon)
    """bdcon[i][j]
        i : 境界番号(x-,x+ ,y-,y+)
        j=1-6 : "bd1_con":"2",
                "bd1_rho":"1.0",
                "bd1_u":"300.0",
                "bd1_v":"0.0",
                "bd1_p":"1.0",
                "bd1_T":"300.0",
    """

    #println(bdcon)
    
    # ---------------------------------
    # bd1 = x-
    # ---------------------------------
    if Int(bdcon[1][1]) == 0
        for j in 1:cellymax
            for k in 1:cellzmax
                for l in 1:5
                    Qbase[1,j,k,l] = bdcon[1][l+1]
                end
            end
        end
    elseif Int(bdcon[1][1]) == 1
        for j in 1:cellymax
            for k in 1:cellzmax
                for l in 1:5
                    Qbase[1,j,k,l] = Qbase[2,j,k,l]
                end
            end
        end
    elseif Int(bdcon[1][1]) == 2
        println("\n We are preparing now! \n")
        throw(UndefVarError(:x))
        for j in 1:cellymax
            for k in 1:cellzmax
                xvec = vecAx[2,j,k,1]
                yvec = vecAx[2,j,k,2]
                zvec = vecAx[2,j,k,3]
                u = Qbase[2,j,k,2]
                v = Qbase[2,j,k,3]
                w = Qbase[2,j,k,4]

                Qbase[1,j,k,1] = Qbase[2,j,k,1]
                Qbase[1,j,k,2] = ((-xvec^2+yvec^2)*u-2*xvec*yvec*v)/(xvec^2+yvec^2)
                Qbase[1,j,k,3] = (-2*xvec*yvec*u+(xvec^2-yvec^2)*v)/(xvec^2+yvec^2)
                Qbase[1,j,k,4] = (-2*xvec*yvec*u+(xvec^2-yvec^2)*v)/(xvec^2+yvec^2)
                Qbase[1,j,k,5] = Qbase[2,j,k,5]
            end
        end
    elseif Int(bdcon[1][1]) == 3
        for j in 1:cellymax
            for k in 1:cellzmax
                u = Qbase[2,j,k,2]
                v = Qbase[2,j,k,3]
                w = Qbase[2,j,k,4]

                Qbase[1,j,k,1] = Qbase[2,j,k,1]
                Qbase[1,j,k,2] = -u
                Qbase[1,j,k,3] = -v
                Qbase[1,j,k,4] = -w
                Qbase[1,j,k,5] = Qbase[2,j,k,5]
            end
        end
    elseif Int(bdcon[1][1]) == 4
        for j in 1:cellymax
            for k in 1:cellzmax
                Qbase[1,j,k,1] = Qbase[cellxmax-1,j,k,1]
                Qbase[1,j,k,2] = Qbase[cellxmax-1,j,k,2]
                Qbase[1,j,k,3] = Qbase[cellxmax-1,j,k,3]
                Qbase[1,j,k,4] = Qbase[cellxmax-1,j,k,4]
                Qbase[1,j,k,5] = Qbase[cellxmax-1,j,k,5]
            end
        end
    elseif Int(bdcon[1][1]) == 4 && Int(bdcon[2][1]) != 4
        println("\n check boundary condition error 4! \n")
        throw(UndefVarError(:x))
    else
        println("\n check boundary condition ! \n")
        throw(UndefVarError(:x))
    end
    
    # ---------------------------------
    # bd2 = x+
    # ---------------------------------
    if Int(bdcon[2][1]) == 0
        for j in 1:cellymax
            for k in 1:cellzmax
                for l in 1:5
                    Qbase[cellxmax,j,k,l] = bdcon[2][l+1]
                end
            end
        end
    elseif Int(bdcon[2][1]) == 1
        for j in 1:cellymax
            for k in 1:cellzmax
                for l in 1:5
                    Qbase[cellxmax,j,k,l] = Qbase[cellxmax-1,j,k,l]
                end
            end
        end
    elseif Int(bdcon[2][1]) == 2
        println("\n We are preparing now! \n")
        throw(UndefVarError(:x))
        for j in 1:cellymax
            xvec = vecAx[cellymax,j,k,1]
            yvec = vecAx[cellymax,j,k,2]
            u = Qbase[cellxmax-1,j,k,2]
            v = Qbase[cellxmax-1,j,k,3]

            Qbase[cellxmax,j,k,1] = Qbase[cellxmax-1,j,k,1]
            Qbase[cellxmax,j,k,2] = ((-xvec^2+yvec^2)*u-2*xvec*yvec*v)/(xvec^2+yvec^2)
            Qbase[cellxmax,j,k,3] = (-2*xvec*yvec*u+(xvec^2-yvec^2)*v)/(xvec^2+yvec^2)
            Qbase[cellxmax,j,k,4] = Qbase[cellxmax-1,j,k,4]  
        end
    elseif Int(bdcon[2][1]) == 3
        for j in 1:cellymax
            for k in 1:cellzmax
                u = Qbase[cellxmax-1,j,k,2]
                v = Qbase[cellxmax-1,j,k,3]
                w = Qbase[cellxmax-1,j,k,4]

                Qbase[cellxmax,j,k,1] = Qbase[cellxmax-1,j,k,1]
                Qbase[cellxmax,j,k,2] = -u
                Qbase[cellxmax,j,k,3] = -v
                Qbase[cellxmax,j,k,4] = -w
                Qbase[cellxmax,j,k,5] = Qbase[cellxmax-1,j,k,5]
            end
        end
    elseif Int(bdcon[2][1]) == 4
        for j in 1:cellymax
            for k in 1:cellzmax
                Qbase[cellxmax,j,k,1] = Qbase[2,j,k,1]
                Qbase[cellxmax,j,k,2] = Qbase[2,j,k,2]
                Qbase[cellxmax,j,k,3] = Qbase[2,j,k,3]
                Qbase[cellxmax,j,k,4] = Qbase[2,j,k,4]
                Qbase[cellxmax,j,k,5] = Qbase[2,j,k,5]
            end
        end
    elseif Int(bdcon[1][1]) != 4 && Int(bdcon[2][1]) == 4
        println("\n check boundary condition error 4! \n")
        throw(UndefVarError(:x))
    else
        println("\n check boundary condition ! \n")
        throw(UndefVarError(:x))
    end

    # ---------------------------------
    # bd3 = y-
    # ---------------------------------
    if Int(bdcon[3][1]) == 0
        for i in 1:cellxmax
            for k in 1:cellzmax
                for l in 1:5
                    Qbase[i,1,k,l] = bdcon[3][l+1]
                end
            end
        end
    elseif Int(bdcon[3][1]) == 1
        for i in 1:cellxmax
            for k in 1:cellzmax
                for l in 1:5
                    Qbase[i,1,k,l] = Qbase[i,2,k,l]
                end
            end
        end
    elseif Int(bdcon[3][1]) == 2
        println("\n We are preparing now! \n")
        throw(UndefVarError(:x))
        for i in 1:cellxmax
            for k in 1:cellzmax
                xvec = vecAy[i,2,1]
                yvec = vecAy[i,2,2]
                u = Qbase[i,2,2]
                v = Qbase[i,2,3]

                Qbase[i,1,k,1] = Qbase[i,2,1]
                Qbase[i,1,k,2] = ((-xvec^2+yvec^2)*u-2*xvec*yvec*v)/(xvec^2+yvec^2)
                Qbase[i,1,k,3] = (-2*xvec*yvec*u+(xvec^2-yvec^2)*v)/(xvec^2+yvec^2)
                Qbase[i,1,k,4] = Qbase[i,2,4]  
            end
        end
    elseif Int(bdcon[3][1]) == 3
        for i in 1:cellxmax
            for k in 1:cellzmax
                u = Qbase[i,2,k,2]
                v = Qbase[i,2,k,3]
                w = Qbase[i,2,k,4]

                Qbase[i,1,k,1] = Qbase[i,2,k,1]
                Qbase[i,1,k,2] = -u
                Qbase[i,1,k,3] = -v
                Qbase[i,1,k,4] = -w
                Qbase[i,1,k,5] = Qbase[i,2,k,5]  
            end
        end
    elseif Int(bdcon[3][1]) == 4
        for i in 1:cellxmax
            for k in 1:cellzmax
                Qbase[i,1,k,1] = Qbase[i,cellymax-1,k,1]  
                Qbase[i,1,k,2] = Qbase[i,cellymax-1,k,2]  
                Qbase[i,1,k,3] = Qbase[i,cellymax-1,k,3]  
                Qbase[i,1,k,4] = Qbase[i,cellymax-1,k,4]
            end
        end
    elseif Int(bdcon[3][1]) == 4 && Int(bdcon[4][1]) != 4
        println("\n check boundary condition error 4! \n")
        throw(UndefVarError(:x))
    else
        println("\n check boundary condition ! \n")
        throw(UndefVarError(:x))
    end
    
    # ---------------------------------
    # bd4 = y+
    # ---------------------------------
    if Int(bdcon[4][1]) == 0
        for i in 1:cellxmax
            for k in 1:cellzmax
                for l in 1:5
                    Qbase[i,cellymax,k,l] = bdcon[4][l+1]
                end
            end
        end
    elseif Int(bdcon[4][1]) == 1
        for i in 1:cellxmax
            for k in 1:cellzmax
                for l in 1:5
                    Qbase[i,cellymax,k,l] = Qbase[i,cellxmax-1,k,l]
                end
            end
        end
    elseif Int(bdcon[4][1]) == 2
        println("\n We are preparing now! \n")
        throw(UndefVarError(:x))
        for i in 1:cellxmax
            xvec = vecAy[i,cellymax,1]
            yvec = vecAy[i,cellymax,2]
            u = Qbase[i,cellymax-1,2]
            v = Qbase[i,cellymax-1,3]

            Qbase[i,cellymax,1] = Qbase[i,cellymax-1,1]
            Qbase[i,cellymax,2] = ((-xvec^2+yvec^2)*u-2*xvec*yvec*v)/(xvec^2+yvec^2)
            Qbase[i,cellymax,3] = (-2*xvec*yvec*u+(xvec^2-yvec^2)*v)/(xvec^2+yvec^2)
            Qbase[i,cellymax,4] = Qbase[i,cellymax-1,4]  
        end
    elseif Int(bdcon[4][1]) == 3
        for i in 1:cellxmax
            for k in 1:cellzmax
                u = Qbase[i,cellymax-1,k,2]
                v = Qbase[i,cellymax-1,k,3]
                w = Qbase[i,cellymax-1,k,4]

                Qbase[i,cellymax,k,1] = Qbase[i,cellymax-1,k,1]
                Qbase[i,cellymax,k,2] = -u
                Qbase[i,cellymax,k,3] = -v
                Qbase[i,cellymax,k,4] = -w
                Qbase[i,cellymax,k,5] = Qbase[i,cellymax-1,k,5]  
            end
        end
    elseif Int(bdcon[4][1]) == 4
        for i in 1:cellxmax
            for k in 1:cellzmax
                Qbase[i,cellymax,k,1] = Qbase[i,2,k,1]
                Qbase[i,cellymax,k,2] = Qbase[i,2,k,2]
                Qbase[i,cellymax,k,3] = Qbase[i,2,k,3]
                Qbase[i,cellymax,k,4] = Qbase[i,2,k,4]
                Qbase[i,cellymax,k,5] = Qbase[i,2,k,5]
            end
        end
    elseif Int(bdcon[3][1]) != 4 && Int(bdcon[4][1]) == 4
        println("\n check boundary condition error 4! \n")
        throw(UndefVarError(:x))
    else
        println("\n check boundary condition ! \n")
        throw(UndefVarError(:x))
    end
    
    # ---------------------------------
    # bd5 = z-
    # ---------------------------------
    if Int(bdcon[5][1]) == 0
        for i in 1:cellxmax
            for j in 1:cellymax
                for l in 1:5
                    Qbase[i,j,1,l] = bdcon[3][l+1]
                end
            end
        end
    elseif Int(bdcon[5][1]) == 1
        for i in 1:cellxmax
            for j in 1:cellymax
                for l in 1:5
                    Qbase[i,j,1,l] = Qbase[i,j,2,l]
                end
            end
        end
    elseif Int(bdcon[5][1]) == 2
        println("\n We are preparing now! \n")
        throw(UndefVarError(:x))
        for i in 1:cellxmax
            for k in 1:cellzmax

            end
        end
    elseif Int(bdcon[5][1]) == 3
        for i in 1:cellxmax
            for j in 1:cellymax
                u = Qbase[i,j,2,2]
                v = Qbase[i,j,2,3]
                w = Qbase[i,j,2,4]

                Qbase[i,j,1,1] = Qbase[i,j,2,1]
                Qbase[i,j,1,2] = -u
                Qbase[i,j,1,3] = -v
                Qbase[i,j,1,4] = -w
                Qbase[i,j,1,5] = Qbase[i,j,2,5]  
            end
        end
    elseif Int(bdcon[5][1]) == 4
        for i in 1:cellxmax
            for j in 1:cellymax
                Qbase[i,j,1,1] = Qbase[i,j,cellzmax-1,1]  
                Qbase[i,j,1,2] = Qbase[i,j,cellzmax-1,2]  
                Qbase[i,j,1,3] = Qbase[i,j,cellzmax-1,3]  
                Qbase[i,j,1,4] = Qbase[i,j,cellzmax-1,4]
                Qbase[i,j,1,5] = Qbase[i,j,cellzmax-1,5]
            end
        end
    elseif Int(bdcon[5][1]) == 4 && Int(bdcon[6][1]) != 4
        println("\n check boundary condition error 4! \n")
        throw(UndefVarError(:x))
    else
        println("\n check boundary condition ! \n")
        throw(UndefVarError(:x))
    end
    
    # ---------------------------------
    # bd6 = z+
    # ---------------------------------
    if Int(bdcon[6][1]) == 0
        for i in 1:cellxmax
            for j in 1:cellymax
                for l in 1:5
                    Qbase[i,j,cellzmax,l] = bdcon[4][l+1]
                end
            end
        end
    elseif Int(bdcon[6][1]) == 1
        for i in 1:cellxmax
            for j in 1:cellymax
                for l in 1:5
                    Qbase[i,j,cellzmax,l] = Qbase[i,cellxmax-1,k,l]
                end
            end
        end
    elseif Int(bdcon[6][1]) == 2
        println("\n We are preparing now! \n")
        throw(UndefVarError(:x))
        for i in 1:cellxmax
            xvec = vecAy[i,cellymax,1]
            yvec = vecAy[i,cellymax,2]
            u = Qbase[i,cellymax-1,2]
            v = Qbase[i,cellymax-1,3]

            Qbase[i,cellymax,1] = Qbase[i,cellymax-1,1]
            Qbase[i,cellymax,2] = ((-xvec^2+yvec^2)*u-2*xvec*yvec*v)/(xvec^2+yvec^2)
            Qbase[i,cellymax,3] = (-2*xvec*yvec*u+(xvec^2-yvec^2)*v)/(xvec^2+yvec^2)
            Qbase[i,cellymax,4] = Qbase[i,cellymax-1,4]  
        end
    elseif Int(bdcon[6][1]) == 3
        for i in 1:cellxmax
            for j in 1:cellymax
                u = Qbase[i,j,cellzmax-1,2]
                v = Qbase[i,j,cellzmax-1,3]
                w = Qbase[i,j,cellzmax-1,4]

                Qbase[i,j,cellzmax,1] = Qbase[i,j,cellzmax-1,1]
                Qbase[i,j,cellzmax,2] = -u
                Qbase[i,j,cellzmax,3] = -v
                Qbase[i,j,cellzmax,4] = -w
                Qbase[i,j,cellzmax,5] = Qbase[i,j,cellzmax-1,5]  
            end
        end
    elseif Int(bdcon[6][1]) == 4
        for i in 1:cellxmax
            for j in 1:cellymax
                Qbase[i,j,cellzmax,1] = Qbase[i,j,2,1]
                Qbase[i,j,cellzmax,2] = Qbase[i,j,2,2]
                Qbase[i,j,cellzmax,3] = Qbase[i,j,2,3]
                Qbase[i,j,cellzmax,4] = Qbase[i,j,2,4]
                Qbase[i,j,cellzmax,5] = Qbase[i,j,2,5]
            end
        end
    elseif Int(bdcon[5][1]) != 4 && Int(bdcon[6][1]) == 4
        println("\n check boundary condition error 4! \n")
        throw(UndefVarError(:x))
    else
        println("\n check boundary condition ! \n")
        throw(UndefVarError(:x))
    end

    return Qbase
end