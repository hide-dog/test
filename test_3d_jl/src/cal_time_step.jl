function time_integration_explicit(dt, Qcon_hat, RHS, cellxmax, cellymax, cellzmax)
    Threads.@threads for i in 2:cellxmax-1
        for j in 2:cellymax-1
            for k in 2:cellzmax
                for l in 1:5
                    Qcon_hat[i,j,k,l] = Qcon_hat[i,j,k,l] + dt*RHS[i,j,k,l]
                end
            end
        end
    end
    return Qcon_hat
end 

function one_wave(Qbase, Qcon, cellxmax, cellymax, cellzmax, vecAx, vecAy, vecAz, specific_heat_ratio, volume)
    # セル中心のonewave
    A_adv_hat_p = zeros(cellxmax,cellymax,cellzmax,5,5)
    A_adv_hat_m = zeros(cellxmax,cellymax,cellzmax,5,5)
    B_adv_hat_p = zeros(cellxmax,cellymax,cellzmax,5,5)
    B_adv_hat_m = zeros(cellxmax,cellymax,cellzmax,5,5)
    C_adv_hat_p = zeros(cellxmax,cellymax,cellzmax,5,5)
    C_adv_hat_m = zeros(cellxmax,cellymax,cellzmax,5,5)
    A_beta_shig = zeros(cellxmax,cellymax,cellzmax)
    B_beta_shig = zeros(cellxmax,cellymax,cellzmax)
    C_beta_shig = zeros(cellxmax,cellymax,cellzmax)
    beta = 1.1

    for i in 2:cellxmax
        for j in 2:cellymax
            for k in 2:cellzmax
                for n in 1:3 #A,B,C
                    jacob_temp = zeros(5,5)
                    if n == 1
                        kx_av = 0.5*(vecAx[i,j,k,1]+vecAx[i+1,j,k,1])
                        ky_av = 0.5*(vecAx[i,j,k,2]+vecAx[i+1,j,k,2])
                        kz_av = 0.5*(vecAx[i,j,k,3]+vecAx[i+1,j,k,3])
                    elseif n == 2
                        kx_av = 0.5*(vecAy[i,j,k,1]+vecAy[i,j+1,k,1])
                        ky_av = 0.5*(vecAy[i,j,k,2]+vecAy[i,j+1,k,2])
                        kz_av = 0.5*(vecAy[i,j,k,3]+vecAy[i,j+1,k,3])
                    elseif n == 3
                        kx_av = 0.5*(vecAz[i,j,k,1]+vecAz[i,j,k+1,1])
                        ky_av = 0.5*(vecAz[i,j,k,2]+vecAz[i,j,k+1,2])
                        kz_av = 0.5*(vecAz[i,j,k,3]+vecAz[i,j,k+1,3])
                    end

                    rho = Qbase[i,j,k,1]
                    u = Qbase[i,j,k,2]
                    v = Qbase[i,j,k,3]
                    w = Qbase[i,j,k,4]
                    e = Qcon[i,j,k,5]
                    p = Qbase[i,j,k,5]

                    g = specific_heat_ratio
                    
                    Z = kx_av*u + ky_av*v + kz_av*w
                    q2   = 0.5*(u^2+v^2)
                    b1c2 = 0.5*q2*(g-1)
                    gebyrho = g*e/rho
                    theta = Z

                    jacob_temp[1,1] = 0.0
                    jacob_temp[1,2] = kx_av
                    jacob_temp[1,3] = ky_av
                    jacob_temp[1,4] = kz_av
                    jacob_temp[1,5] = 0.0
                
                    jacob_temp[2,1] = -u*theta + kx_av*b1c2
                    jacob_temp[2,2] = theta    - kx_av*(g-2)*u
                    jacob_temp[2,3] = ky_av*u  - kx_av*(g-1)*v
                    jacob_temp[2,4] = kz_av*u  - kx_av*(g-1)*w
                    jacob_temp[2,5] = kx_av*(g-1)
                
                    jacob_temp[3,1] = -v*theta + ky_av*b1c2
                    jacob_temp[3,2] = kx_av*v  - ky_av*(g-1)*u
                    jacob_temp[3,3] = theta    - ky_av*(g-2)*v
                    jacob_temp[3,4] = kz_av*v  - ky_av*(g-1)*w
                    jacob_temp[3,5] = ky_av*(g-1)
                
                    jacob_temp[4,1] = -w*theta + kx_av*b1c2
                    jacob_temp[4,2] = kx_av*w  - kz_av*(g-1)*u
                    jacob_temp[4,3] = ky_av*w  - kz_av*(g-1)*v
                    jacob_temp[4,4] = theta    - kz_av*(g-2)*w
                    jacob_temp[4,5] = kz_av*(g-1)

                    jacob_temp[5,1] = Z*(-gebyrho + 2*b1c2)
                    jacob_temp[5,2] = kx_av*(gebyrho-b1c2) - (g-1)*u*theta
                    jacob_temp[5,3] = ky_av*(gebyrho-b1c2) - (g-1)*v*theta
                    jacob_temp[5,4] = kz_av*(gebyrho-b1c2) - (g-1)*w*theta
                    jacob_temp[5,5] = g*theta

                    c = (g*rho/p)^0.5
                    shigma = abs(Z) + c*(kx_av^2+ky_av^2)^0.5

                    I_temp = zeros(5,5)
                    for l in 1:5
                        I_temp[l,l] = beta * shigma
                    end

                    if k == 1
                        for l in 1:5
                            for m in 1:5
                                A_adv_hat_p[i,j,k,l,m] = 0.5*(jacob_temp[l,m] + I_temp[l,m])
                                A_adv_hat_m[i,j,k,l,m] = 0.5*(jacob_temp[l,m] - I_temp[l,m])
                            end
                        end
                        A_beta_shig[i,j,k] = beta * shigma

                    elseif k == 2
                        for l in 1:5
                            for m in 1:5
                                B_adv_hat_p[i,j,k,l,m] = 0.5*(jacob_temp[l,m] + I_temp[l,m])
                                B_adv_hat_m[i,j,k,l,m] = 0.5*(jacob_temp[l,m] - I_temp[l,m])
                            end
                        end
                        B_beta_shig[i,j,k] = beta * shigma

                    elseif k == 3
                        for l in 1:5
                            for m in 1:5
                                C_adv_hat_p[i,j,k,l,m] = 0.5*(jacob_temp[l,m] + I_temp[l,m])
                                C_adv_hat_m[i,j,k,l,m] = 0.5*(jacob_temp[l,m] - I_temp[l,m])
                            end
                        end
                        C_beta_shig[i,j,k] = beta * shigma
                    end
                end
            end
        end
    end

    return A_adv_hat_p, A_adv_hat_m, B_adv_hat_p, B_adv_hat_m, C_adv_hat_p, C_adv_hat_m, A_beta_shig, B_beta_shig, C_beta_shig
end

function central_diff_jacobian(Qbase, Qcon, cellxmax, cellymax, cellzmax, mu, lambda, vecAx, vecAy, vecAz, specific_heat_ratio, volume)
    jalphaP = zeros(cellxmax, cellymax, cellzmax)
    jbetaP  = zeros(cellxmax, cellymax, cellzmax)
    jgammaP = zeros(cellxmax, cellymax, cellzmax)

    # 角セルの処理
    tempx = [1 cellxmax]
    tempy = [1 cellymax]
    tempz = [1 cellzmax]
    for i in 1:2
        for j in 1:2 
            for k in 1:2
                Qbase[tempx[i],tempy[j],tempz[k],1] = 1.0
                volume[tempx[i],tempy[j],tempz[k]]  = 1.0
            end
        end
    end
    
    for i in 1:cellxmax
        for j in 1:cellymax
            for k in 1:cellzmax
                rho = Qbase[i,j,k,1]

                xi_x_av = 0.5*(vecAx[i,j,k,1] + vecAx[i+1,j,k,1]) 
                xi_y_av = 0.5*(vecAx[i,j,k,2] + vecAx[i+1,j,k,2]) 
                xi_z_av = 0.5*(vecAx[i,j,k,3] + vecAx[i+1,j,k,3]) 
                #alpha = (xi_x_av^2+xi_y_av^2)^0.5
                alpha = xi_x_av^2 + xi_y_av^2 + xi_z_av^2
                jalphaP[i,j,k] = alpha / volume[i,j,k] * (2*mu[i,j,k]/rho)

                eta_x_av = 0.5*(vecAy[i,j,k,1] + vecAy[i,j+1,k,1])
                eta_y_av = 0.5*(vecAy[i,j,k,2] + vecAy[i,j+1,k,2])
                eta_z_av = 0.5*(vecAy[i,j,k,3] + vecAy[i,j+1,k,3])
                #beta = (xi_x_av^2+xi_y_av^2)^0.5
                beta = eta_x_av^2 + eta_y_av^2 + eta_z_av^2
                jbetaP[i,j,k] = beta / volume[i,j,k] * (2*mu[i,j,k]/rho)

                zeta_x_av = 0.5*(vecAz[i,j,k,1] + vecAz[i,j,k+1,1])
                zeta_y_av = 0.5*(vecAz[i,j,k,2] + vecAz[i,j,k+1,2])
                zeta_z_av = 0.5*(vecAz[i,j,k,3] + vecAz[i,j,k+1,3])
                gamma = zeta_x_av^2 + zeta_y_av^2 + zeta_z_av^2
                jgammaP[i,j,k] = gamma / volume[i,j,k] * (2*mu[i,j,k]/rho)
            end
        end
    end
    return jalphaP, jbetaP, jgammaP
end

function lusgs(dt, Qcon_hat, Qconn_hat, delta_Q, A_adv_hat_p, A_adv_hat_m, B_adv_hat_p, B_adv_hat_m, C_adv_hat_p, C_adv_hat_m,
    A_beta_shig, B_beta_shig, C_beta_shig, jalphaP, jbetaP, jgammaP, RHS, cellxmax, cellymax,cellzmax, volume)
    
    D  = zeros(cellxmax, cellymax, cellzmax)
    Lx = zeros(cellxmax, cellymax, cellzmax, 5, 5)
    Ly = zeros(cellxmax, cellymax, cellzmax, 5, 5)
    Lz = zeros(cellxmax, cellymax, cellzmax, 5, 5)
    Ux = zeros(cellxmax, cellymax, cellzmax, 5, 5)
    Uy = zeros(cellxmax, cellymax, cellzmax, 5, 5)
    Uz = zeros(cellxmax, cellymax, cellzmax, 5, 5)
    
    I = zeros(5, 5)
    for i in 1:5
        I[i,i] = 1.0
    end
    
    # L,U
    for i in 1:cellxmax
        for j in 1:cellymax
            for k in 1:cellzmax
                for l in 1:5
                    for m in 1:5
                        Lx[i,j,k,l,m] = dt*(A_adv_hat_p[i,j,k,l,m] + jalphaP[i,j,k]*I[l,m])
                        Ly[i,j,k,l,m] = dt*(B_adv_hat_p[i,j,k,l,m] + jbetaP[i,j,k]*I[l,m])
                        Lz[i,j,k,l,m] = dt*(C_adv_hat_p[i,j,k,l,m] + jgammaP[i,j,k]*I[l,m])
                        Ux[i,j,k,l,m] = dt*(A_adv_hat_m[i,j,k,l,m] - jalphaP[i,j,k]*I[l,m])
                        Uy[i,j,k,l,m] = dt*(B_adv_hat_m[i,j,k,l,m] - jbetaP[i,j,k]*I[l,m])
                        Uz[i,j,k,l,m] = dt*(C_adv_hat_m[i,j,k,l,m] - jgammaP[i,j,k]*I[l,m])
                    end
                end
            end
        end
    end

    LdQ = zeros(cellxmax, cellymax, cellzmax, 5)
    UdQ = zeros(cellxmax, cellymax, cellzmax, 5)
    for i in 2:cellxmax-1
        for j in 2:cellymax-1
            for k in 2:cellzmax-1
                for l in 1:5
                    for m in 1:5
                        LdQ[i,j,k,l] += Lx[i-1,j,k,l,m]*delta_Q[i-1,j,k,m] +
                                        Ly[i,j-1,k,l,m]*delta_Q[i,j-1,k,m] +
                                        Lz[i,j,k-1,l,m]*delta_Q[i,j,k-1,m]
                        UdQ[i,j,k,l] += Ux[i+1,j,k,l,m]*delta_Q[i+1,j,k,m] +
                                        Uy[i,j+1,k,l,m]*delta_Q[i,j+1,k,m] +
                                        Uy[i,j,k+1,l,m]*delta_Q[i,j,k+1,m]
                    end
                end
            end
        end
    end

    # diagonal
    for i in 2:cellxmax-1
        for j in 2:cellymax-1
            for k in 2:cellzmax-1
                D[i,j,k] = 2.0 + dt*(A_beta_shig[i,j,k]+2*jalphaP[i,j,k] + B_beta_shig[i,j,k]+2*jbetaP[i,j,k])
            end
            #D[i,j,k] = volume[i,j,k]*2 + dt*(A_beta_shig[i,j,k]+2*jalphaP[i,j,k] + B_beta_shig[i,j,k]+2*jbetaP[i,j,k])
        end
    end
    
    
    #println(volume[50,50])
    #println(A_beta_shig[50,50])
    #println(B_beta_shig[50,50])
    #println(jalphaP[50,50])
    #println(jbetaP[50,50])
    
    
    RHS_temp = zeros(cellxmax, cellymax, cellzmax, 5)
    # RHS
    for i in 2:cellxmax-1
        for j in 2:cellymax-1
            for k in 2:cellzmax-1
                for l in 1:5
                    RHS_temp[i,j,k,l] = - (Qcon_hat[i,j,k,l]-Qconn_hat[i,j,k,l])*1.0 + dt*RHS[i,j,k,l]
                end
            end
        end
    end

    #println(RHS_temp[:,50,:])
    #throw(UndefVarError(:x))

    # lower sweep
    for i in 2:cellxmax-1
        for j in 2:cellymax-1
            for k in 2:cellzmax-1
                for l in 1:5
                    delta_Q[i,j,k,l] = D[i,j,k]^(-1) * (RHS_temp[i,j,k,l]+LdQ[i,j,k,l])
                end
            end
        end
    end             

    # upepr sweep
    for i in 2:cellxmax-1
        for j in 2:cellymax-1
            for k in 2:cellzmax-1
                for l in 1:5
                    delta_Q[i,j,k,l] = delta_Q[i,j,k,l] - D[i,j,k]^(-1) * UdQ[i,j,k,l]
                end
            end
        end
    end
    
    #println(Qcon_hat[:,50,:])
    #println(delta_Q[:,50,:])
    #throw(UndefVarError(:x))

    return delta_Q
end

