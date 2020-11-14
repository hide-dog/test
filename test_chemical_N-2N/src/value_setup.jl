function set_volume(nodes, cellxmax, cellymax)
    volume = zeros(cellxmax, cellymax)
    for i in 1:cellxmax
        for j in 1:cellymax
            vec_r1x = nodes[i+1,j+1,1] - nodes[i,j,1]
            vec_r1y = nodes[i+1,j+1,2] - nodes[i,j,2]
            vec_r2x = nodes[i,j+1,1] - nodes[i+1,j,1]
            vec_r2y = nodes[i,j+1,2] - nodes[i+1,j,2]

            volume[i,j] = abs(vec_r1x*vec_r2y - vec_r1y*vec_r2x) /2
        end
    end
    return volume
end 

function set_dx_lts(nodes, cellxmax, cellymax)
    # 境界で定義
    dx = zeros(cellxmax+1, cellymax)
    dy = zeros(cellxmax, cellymax+1)
    
    for i in 2:cellxmax+1 -1
        for j in 2:cellymax -1
            x1 = nodes[i,j,1]
            y1 = nodes[i,j,2]
            x2 = nodes[i,j+1,1]
            y2 = nodes[i,j+1,2]

            a  = (y2-y1) / (x2-x1)
            b  = y1 - a*x1

            ccx = 0.125 * (nodes[i,j,1] + nodes[i+1,j,1] + nodes[i,j+1,1] + nodes[i+1,j+1,1])
            ccy = 0.125 * (nodes[i,j,2] + nodes[i+1,j,2] + nodes[i,j+1,2] + nodes[i+1,j+1,2])
            
            dx1 = abs(-a*ccx + ccy -b)/abs(a)

            ccx = 0.125 * (nodes[i-1,j,1] + nodes[i,j,1] + nodes[i-1,j+1,1] + nodes[i,j+1,1])
            ccy = 0.125 * (nodes[i-1,j,2] + nodes[i,j,2] + nodes[i-1,j+1,2] + nodes[i,j+1,2])
            
            dx2 = abs(-a*ccx + ccy -b)/abs(a)

            dx[i,j] = 0.5 * (dx1+dx2)
        end
    end

    for i in 2:cellxmax -1
        for j in 2:cellymax+1 -1
            x1 = nodes[i,j,1]
            y1 = nodes[i,j,2]
            x2 = nodes[i+1,j,1]
            y2 = nodes[i+1,j,2]

            a  = (y2-y1) / (x2-x1)
            b  = y1 - a*x1

            ccx = 0.125 * (nodes[i,j,1] + nodes[i+1,j,1] + nodes[i,j+1,1] + nodes[i+1,j+1,1])
            ccy = 0.125 * (nodes[i,j,2] + nodes[i+1,j,2] + nodes[i,j+1,2] + nodes[i+1,j+1,2])
            
            dy1 = abs(-a*ccx + ccy -b)/abs(a)

            ccx = 0.125 * (nodes[i,j,1] + nodes[i+1,j,1] + nodes[i,j-1,1] + nodes[i+1,j-1,1])
            ccy = 0.125 * (nodes[i,j,2] + nodes[i+1,j,2] + nodes[i,j-1,2] + nodes[i+1,j-1,2])
            
            dy2 = abs(-a*ccx + ccy -b)/abs(a)

            dy[i,j] = 0.5 * (dy1+dy2)
        end
    end

    return dx, dy
end

function set_lts(Qbase, cellxmax, cellymax, mu, dx, dy, vecAx, vecAy, volume, specific_heat_ratio, cfl)
    dtau = zeros(cellxmax, cellymax)
    lambda_facex = zeros(cellxmax+1, cellymax)
    lambda_facey = zeros(cellxmax, cellymax+1)
    g = specific_heat_ratio

    for i in 2:cellxmax+1 -1
        for j in 2:cellymax -1

            rho_av = 0.5 * (Qbase[i,j,1] + Qbase[i-1,j,1])
            u_av   = 0.5 * (Qbase[i,j,2] + Qbase[i-1,j,2])
            v_av   = 0.5 * (Qbase[i,j,3] + Qbase[i-1,j,3])
            mu_av  = 0.5 * (   mu[i,j]   +    mu[i-1,j]  )
            
            ap = (g * Qbase[  i,j,4] / Qbase[  i,j,1])^0.5
            am = (g * Qbase[i-1,j,4] / Qbase[i-1,j,1])^0.5
            a_av = 0.5 * (ap + am)

            U   = u_av*vecAx[i,j,1] + v_av*vecAx[i,j,2]
            lambda_facex[i,j] = abs(U) + a_av + 2*mu_av/(rho_av*dx[i,j])
        end
    end
    for i in 2:cellxmax -1
        for j in 2:cellymax+1 -1

            rho_av = 0.5 * (Qbase[i,j,1] + Qbase[i,j-1,1])
            u_av   = 0.5 * (Qbase[i,j,2] + Qbase[i,j-1,2])
            v_av   = 0.5 * (Qbase[i,j,3] + Qbase[i,j-1,3])
            mu_av  = 0.5 * (   mu[i,j]   +    mu[i,j-1]  )
            
            ap = (g * Qbase[  i,j,4] / Qbase[  i,j,1])^0.5
            am = (g * Qbase[i,j-1,4] / Qbase[i,j-1,1])^0.5
            a_av = 0.5 * (ap + am)

            V   = u_av*vecAy[i,j,1] + v_av*vecAy[i,j,2]
            lambda_facey[i,j] = abs(V) + a_av + 2*mu_av/(rho_av*dy[i,j])
        end
    end

    for i in 2:cellxmax-1
        for j in 2:cellymax-1
            a1 = lambda_facex[  i,  j]
            a2 = lambda_facex[i+1,  j]
            a3 = lambda_facey[  i,  j]
            a4 = lambda_facey[  i,j+1]
            lmax = maximum([a1,a2,a3,a4])

            dtau[i,j] = cfl * volume[i,j] / lmax
        end
    end

    return dtau
end


function set_mu(Qbase,cellxmax,cellymax,specific_heat_ratio,Rd)
    mu = zeros(cellxmax,cellymax)

    # サザーランドの式
    # https://www.jstage.jst.go.jp/article/jsam1937/37/4/37_4_694/_pdf/-char/ja
    mu0 = 1.82e-5 # 基準粘度[Pa s]
    T0 = 293.15       # 基準温度[K]
    C = 117       # サザーランド定数[K]

    for i in 1:cellxmax
        for j in 1:cellymax
            # T = p/(rho Rd )
            T = Qbase[i,j,4]/(Qbase[i,j,1]*Rd)

            if T <0
                println("T is munus ( 'A' ) !!")
                println("i : "*string(i)*"  j : "*string(j))
                println(Qbase[i,j,:])
            end
            
            mu[i,j] = mu0 * (T/T0)^1.5 * (T0+C)/(T+C)
        end
    end

    
    return mu
end


function set_lambda(Qbase,cellxmax,cellymax,mu,specific_heat_ratio,Rd)
    
    lambda = zeros(cellxmax,cellymax)
    
    # サザーランドの式
    lam0 = 22.3*10^(-3)  # 基準熱伝導率　[W/mK]
    T0 = 273.15       # 基準温度[K]
    C = 125       # サザーランド定数[K]

    for i in 1:cellxmax
        for j in 1:cellymax
            T = Qbase[i,j,4]/(Qbase[i,j,1]*Rd)
            lambda[i,j] = lam0*((T0+C)/(T+C))*(T/T0)^1.5
        end
    end

    return lambda
end

function set_Minf(bdcon, specific_heat_ratio)
    rho = 0 
    u   = 0 
    v   = 0 
    p   = 0 
    Rd  = 287.0

    for i in 1:4
        if Int(bdcon[i][1]) == 0 || Int(bdcon[i][1]) == 5
            rho = bdcon[i][2]
            u   = bdcon[i][3]
            v   = bdcon[i][4]
            p   = bdcon[i][5]
        elseif Int(bdcon[i][1]) == 6
            rho = bdcon[i][2]
            u   = bdcon[i][3]
            v   = bdcon[i][4]
            T   = bdcon[i][8]
            p   = rho*Rd*T
        end
    end

    a = (specific_heat_ratio * p / rho)^0.5
    u = (u^2 + v^2)^0.5
    M = u/a

    return M
end

function set_gasconst_hat(Qbase,cellxmax,cellymax,nval,nch,R)
    Rhat = zeros(cellxmax, cellymax)
    mw   = set_mw()
    npre = nval - nch

    for i in 1:cellxmax
        for j in 1:cellymax
            for s in 1:nch
                Rs = R/mw[s]
                Rhat[i,j] += Qbase[i,j,npre+s] * Rs
            end
        end
    end
    return Rhat
end