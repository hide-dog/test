function set_volume(nodes, vecAx, vecAy, vecAz, cellxmax, cellymax, cellzmax)
    volume = zeros(cellxmax, cellymax, cellzmax)
    for i in 1:cellxmax
        for j in 1:cellymax
            for k in 1:cellzmax
                vec_rx = nodes[i+1,j+1,k+1,1] - nodes[i,j,k,1]
                vec_ry = nodes[i+1,j+1,k+1,2] - nodes[i,j,k,2]
                vec_rz = nodes[i+1,j+1,k+1,3] - nodes[i,j,k,3]

                Ax = vecAx[i,j,k,1] + vecAy[i,j,k,1] + vecAz[i,j,k,1]
                Ay = vecAx[i,j,k,2] + vecAy[i,j,k,2] + vecAz[i,j,k,2]
                Az = vecAx[i,j,k,3] + vecAy[i,j,k,3] + vecAz[i,j,k,3]

                volume[i,j,k] = abs(Ax*vec_rx + Ay*vec_ry + Az*vec_rz) /3
            end
        end
    end
    return volume
end 

function set_cellcenter(nodes, cellxmax, cellymax, cellzmax)
    cellcenter = zeros(cellxmax, cellymax, cellzmax, 3)
    for i in 1:cellxmax
        for j in 1:cellymax
            for k in 1:cellzmax
                for l in 1:3
                    xl = 0.125*(nodes[i+1,j+1,k+1,l] + nodes[i+1,j,k+1,l] + nodes[i,j+1,k+1,l] + nodes[i,j,k+1,l] +
                                nodes[i+1,j+1,k,l] + nodes[i+1,j,k,l] + nodes[i,j+1,k,l] + nodes[i,j,k,l])
                    
                    cellcenter[i,j,k,l] = xl
                end
            end
        end
    end
    return cellcenter
end

function set_mu(Qbase, cellxmax, cellymax, cellzmax, specific_heat_ratio, Rd)
    mu = zeros(cellxmax, cellymax, cellzmax)

    # サザーランドの式
    # https://www.jstage.jst.go.jp/article/jsam1937/37/4/37_4_694/_pdf/-char/ja
    mu0 = 1.82e-5 # 基準粘度[Pa s]
    T0 = 293.15   # 基準温度[K]
    C = 117       # サザーランド定数[K]

    for i in 1:cellxmax
        for j in 1:cellymax
            for k in 1:cellzmax
                # T = p/(rho Rd )
                T = Qbase[i,j,k,5]/(Qbase[i,j,k,1]*Rd)

                if T <0
                    println("T is munus ( 'A' ) !!")
                    println("i : "*string(i)*"  j : "*string(j)*"  k : "*string(k))
                    println(Qbase[i,j,k,:])
                end
                
                mu[i,j,k] = mu0 * (T/T0)^1.5 * (T0+C)/(T+C)
            end
        end
    end
    
    return mu
end


function set_lambda(Qbase, cellxmax, cellymax, cellzmax, mu, specific_heat_ratio, Rd)
    
    lambda = zeros(cellxmax, cellymax, cellzmax)
    
    # サザーランドの式
    lam0 = 22.3*10^(-3)  # 基準熱伝導率　[W/mK]
    T0 = 273.15       # 基準温度[K]
    C = 125       # サザーランド定数[K]

    for i in 1:cellxmax
        for j in 1:cellymax
            for k in 1:cellzmax
                T = Qbase[i,j,k,5]/(Qbase[i,j,k,1]*Rd)
                lambda[i,j,k] = lam0*((T0+C)/(T+C))*(T/T0)^1.5
            end
        end
    end

    return lambda
end