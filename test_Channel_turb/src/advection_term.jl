function AUSM(Qbase,Qcon,cellxmax,cellymax,vecAx,vecAy,specific_heat_ratio)
    E_adv_hat = zeros(cellxmax+1,cellymax)
    F_adv_hat = zeros(cellxmax,cellymax+1)
    
    for i in 2:cellxmax+1 -1
        for j in 2:cellymax -1
            # i-1セル
            cell_vecAx = zeros(2)
            cell_vecAx[1] = 0.5*(vecAx[i,j,1]+vecAx[i+1,j,1])
            cell_vecAx[2] = 0.5*(vecAx[i,j,2]+vecAx[i+1,j,2])

            rho = Qbase[i-1,j,1]
            U = Qbase[i-1,j,2]*cell_vecAx[1] + Qbase[i-1,j,3]*cell_vecAx[2]
            p = Qbase[i-1,j,4]
            a = (specific_heat_ratio * rho / p)^0.5
            M = u/a

            if M <= 1
                u_p = 0.25 * M*a * (M+1)^2
                p_p = 0.25 *p * (1+M)^2 *(2-M)
            elseif M > 1
                u_p = 0.5 * a * (M+abs(M))^2
                p_p = 0.5 *p * (M+abs(M))/M
            end

            # iセル
            cell_vecAx = zeros(2)
            cell_vecAx[1] = 0.5*(vecAx[i+1,j,1]+vecAx[i+2,j,1])
            cell_vecAx[2] = 0.5*(vecAx[i+1,j,2]+vecAx[i+2,j,2])

            rho = Qbase[i,j,1]
            U = Qbase[i,j,2]*cell_vecAx[1] + Qbase[i,j,3]*cell_vecAx[2]
            p = Qbase[i,j,4]
            a = (specific_heat_ratio * rho / p)^0.5
            M = u/a

            if M <= 1
                u_m = -0.25 * M*a * (M-1)^2
                p_m = 0.25 *p * (1-M)^2 *(2+M)
            elseif M > 1
                u_m = 0.5 * a * (M-abs(M))^2
                p_m = 0.5 *p * (M-abs(M))/M
            end

            # bound
            u_half = u_p + u_m
            p_half = p_p + p_m
            temp_vecX = zeros(4)
            temp_vecX[2] = vecAx[i,j,1]
            temp_vecX[3] = vecAx[i,j,2]

            for k in 1:4
                #----------------------------------------
                E_adv_hat[i,j,k] = 0.5 * ((u_half+abs(u_half))*Qcon_hat[i-1,j,k]+(u_half-abs(u_half))*Qcon_hat[i,j,k]) + p_half * temp_vecX[k]
            end
        end
    end

    for i in 2:cellxmax -1
        for j in 2:cellymax+1 -1
            # j-1セル
            cell_vecAy = zeros(2)
            cell_vecAy[1] = 0.5*(vecAy[i,j,1]+vecAy[i,j+1,1])
            cell_vecAy[2] = 0.5*(vecAy[i,j,2]+vecAy[i,j+1,2])

            rho = Qbase[i,j-1,1]
            #-------------
            V = Qbase[i,j-1,2]*cell_vecAy[1] + Qbase[i,j-1,3]*cell_vecAy[2]
            p = Qbase[i,j-1,4]
            a = (specific_heat_ratio * rho / p)^0.5
            M = u/a

            if M <= 1
                v_p = 0.25 * M*a * (M+1)^2
                p_p = 0.25 *p * (1+M)^2 *(2-M)
            elseif M > 1
                v_p = 0.5 * a * (M+abs(M))^2
                p_p = 0.5 *p * (M+abs(M))/M
            end

            # iセル
            cell_vecAy = zeros(2)
            cell_vecAy[1] = 0.5*(vecAy[i,j+1,1]+vecAy[i,j+2,1])
            cell_vecAy[2] = 0.5*(vecAy[i,j+1,2]+vecAy[i,j+2,2])

            rho = Qbase[i,j,1]
            V = Qbase[i,j,2]*cell_vecAy[1] + Qbase[i,j,3]*cell_vecAy[2]
            p = Qbase[i,j,4]
            a = (specific_heat_ratio * rho / p)^0.5
            M = u/a

            if M <= 1
                v_m = -0.25 * M*a * (M-1)^2
                p_m = 0.25 *p * (1-M)^2 *(2+M)
            elseif M > 1
                v_m = 0.5 * a * (M-abs(M))^2
                p_m = 0.5 *p * (M-abs(M))/M
            end

            # bound
            v_half = v_p + v_m
            p_half = p_p + p_m
            temp_vecX = zeros(4)
            temp_vecX[2] = vecAy[i,j,1]
            temp_vecX[3] = vecAy[i,j,2]

            for k in 1:4
                #----------------------------------------
                F_adv_hat[i,j,k] = 0.5 * ((v_half+abs(v_half))*Qcon_hat[i,j-1,k]+(v_half-abs(v_half))*Qcon_hat[i,j,k]) + p_half * temp_vecX[k]
            end
        end
    end

    return E_adv_hat, F_adv_hat
end

function setup_cell_flux_hat(Qcon, cellxmax, cellymax, cellzmax, cell_Ahat_plas, cell_Ahat_minus, cell_Bhat_plas, cell_Bhat_minus, cell_Chat_plas, cell_Chat_minus)
    cell_E_hat_plas  = zeros(cellxmax, cellymax, cellzmax, 5)
    cell_E_hat_minus = zeros(cellxmax, cellymax, cellzmax, 5)
    cell_F_hat_plas  = zeros(cellxmax, cellymax, cellzmax, 5)
    cell_F_hat_minus = zeros(cellxmax, cellymax, cellzmax, 5)
    cell_G_hat_plas  = zeros(cellxmax, cellymax, cellzmax, 5)
    cell_G_hat_minus = zeros(cellxmax, cellymax, cellzmax, 5)

    for i in 1:cellxmax
        for j in 1:cellymax
            for k in 1:cellzmax
                for l in 1:5
                    for m in 1:5
                        # Ehat = Ahat * Qhat = Ahat * Q/J = 1/J*Ahat * Q
                        # 1/J*Ahatの計算をvecAを使ったため，ここではQconを掛けてる
                        cell_E_hat_plas[i,j,k,l]  += cell_Ahat_plas[i,j,k,l,m]*Qcon[i,j,k,m]
                        cell_E_hat_minus[i,j,k,l] += cell_Ahat_minus[i,j,k,l,m]*Qcon[i,j,k,m]
                        cell_F_hat_plas[i,j,k,l]  += cell_Bhat_plas[i,j,k,l,m]*Qcon[i,j,k,m]
                        cell_F_hat_minus[i,j,k,l] += cell_Bhat_minus[i,j,k,l,m]*Qcon[i,j,k,m]
                        cell_G_hat_plas[i,j,k,l]  += cell_Chat_plas[i,j,k,l,m]*Qcon[i,j,k,m]
                        cell_G_hat_minus[i,j,k,l] += cell_Chat_minus[i,j,k,l,m]*Qcon[i,j,k,m]
                    end
                end
            end
        end
    end

    return cell_E_hat_plas, cell_E_hat_minus, cell_F_hat_plas, cell_F_hat_minus, cell_G_hat_plas, cell_G_hat_minus
end

function cal_jacobi(Qbase, Qcon, cellxmax, cellymax, cellzmax, specific_heat_ratio, vecAx, vecAy, vecAz, volume)
    """
    Qbase=[rho,u,v,p]
    Qcon=[rho,rhou,rhov,e]
    
    q=u^2+v^2
    c=(rp/rho)^0.5
    theta=Z=U= A^xi_x * u   +   A^xi_y * v
    H=(e+p)/rho
    b1=q^2/2 * (r-1)/(c^2)
    b2=(r-1)/(c^2)

    xi_x_bar= xi_x / (xi_x^2+xi_y^2)^0.5
    xi_y_bar= xi_y / (xi_x^2+xi_y^2)^0.5
    U_bar=Z_bar=theta_bar=


    A = Rleft * Lambda * Rright
    """

    cell_Ahat_plas  = zeros(cellxmax, cellymax, cellzmax, 5, 5)
    cell_Ahat_minus = zeros(cellxmax, cellymax, cellzmax, 5, 5)
    cell_Bhat_plas  = zeros(cellxmax, cellymax, cellzmax, 5, 5)
    cell_Bhat_minus = zeros(cellxmax, cellymax, cellzmax, 5, 5)
    cell_Chat_plas  = zeros(cellxmax, cellymax, cellzmax, 5, 5)
    cell_Chat_minus = zeros(cellxmax, cellymax, cellzmax, 5, 5)
    r = specific_heat_ratio

    for i in 1:cellxmax
        for j in 1:cellymax
            for k in 1:cellzmax
                q = (Qbase[i,j,k,2]^2+Qbase[i,j,k,3]^2+Qbase[i,j,k,4]^2)^0.5
                c = 0
                try c = (r*Qbase[i,j,k,5]/Qbase[i,j,k,1])^0.5
                catch
                    println("\n"*string(i)*","*string(j)*","*string(k)*" Qbase error")
                    println(Qbase[i,j,k,:])
                    println("\n")
                    throw(UndefVarError(:x))
                end

                H  = (Qcon[i,j,k,5]+Qbase[i,j,k,5])/Qbase[i,j,k,1]
                b1 = q^2/2 * (r-1)/(c^2)
                b2 = (r-1)/(c^2)
                u  = Qbase[i,j,k,2]
                v  = Qbase[i,j,k,3]
                w  = Qbase[i,j,k,4]

                # cell_Ahat
                k_x = 0.5*(vecAx[i,j,k,1]+vecAx[i+1,j,k,1])
                k_y = 0.5*(vecAx[i,j,k,2]+vecAx[i+1,j,k,2])
                k_z = 0.5*(vecAx[i,j,k,3]+vecAx[i+1,j,k,3])
                
                Rleft_E, Lambda_E, Rright_E = Eigenvalue_vector(q,c,H,b1,b2,u,v,w,k_x,k_y,k_z)
                A_plas, A_minus = RLmbdaR(Rleft_E,Lambda_E,Rright_E)

                for l in 1:5
                    for m in 1:5
                        cell_Ahat_plas[i,j,k,l,m]  = A_plas[l,m]
                        cell_Ahat_minus[i,j,k,l,m] = A_minus[l,m]
                    end
                end
                
                # cell_Bhat
                k_x = 0.5*(vecAy[i,j,k,1]+vecAy[i,j+1,k,1])
                k_y = 0.5*(vecAy[i,j,k,2]+vecAy[i,j+1,k,2])
                k_z = 0.5*(vecAy[i,j,k,3]+vecAy[i,j+1,k,3])

                Rleft_E, Lambda_E, Rright_E = Eigenvalue_vector(q,c,H,b1,b2,u,v,w,k_x,k_y,k_z)
                A_plas, A_minus = RLmbdaR(Rleft_E,Lambda_E,Rright_E)

                for l in 1:5
                    for m in 1:5
                        cell_Bhat_plas[i,j,k,l,m]  = A_plas[l,m]
                        cell_Bhat_minus[i,j,k,l,m] = A_minus[l,m]
                    end
                end

                # cell_Chat
                k_x = 0.5*(vecAz[i,j,k,1]+vecAz[i,j,k+1,1])
                k_y = 0.5*(vecAz[i,j,k,2]+vecAz[i,j,k+1,2])
                k_z = 0.5*(vecAz[i,j,k,3]+vecAz[i,j,k+1,3])

                Rleft_E, Lambda_E, Rright_E = Eigenvalue_vector(q,c,H,b1,b2,u,v,w,k_x,k_y,k_z)
                A_plas, A_minus = RLmbdaR(Rleft_E,Lambda_E,Rright_E)

                for l in 1:5
                    for m in 1:5
                        cell_Chat_plas[i,j,k,l,m]  = A_plas[l,m]
                        cell_Chat_minus[i,j,k,l,m] = A_minus[l,m]
                    end
                end
            end
        end
    end

    return cell_Ahat_plas, cell_Ahat_minus, cell_Bhat_plas, cell_Bhat_minus, cell_Chat_plas, cell_Chat_minus
end


function Eigenvalue_vector(q,c,H,b1,b2,u,v,w,k_x,k_y,k_z)
    Z       = k_x*u + k_y*v + k_z*w
    k_x_bar = k_x / (k_x^2 + k_y^2 + k_z^2)^0.5
    k_y_bar = k_y / (k_x^2 + k_y^2 + k_z^2)^0.5
    k_z_bar = k_z / (k_x^2 + k_y^2 + k_z^2)^0.5
    Z_bar   =   Z / (k_x^2 + k_y^2 + k_z^2)^0.5
    
    Lambda = [Z - c*(k_x^2+k_y^2+k_z^2)^0.5 
              Z
              Z
              Z
              Z + c*(k_x^2+k_y^2+k_z^2)^0.5]
    

    temp1 = k_y_bar*q^2/2+k_z_bar*u*c-k_x_bar*w*c
    temp2 = k_z_bar*q^2/2-k_y_bar*u*c+k_x_bar*v*c
    temp3 = k_x_bar*q^2/2-k_z_bar*v*c+k_y_bar*w*c

    Rleft = [[1           k_y_bar             k_z_bar             k_x_bar             1          ]
             [u-k_x_bar*c k_y_bar*u+k_z_bar*c k_z_bar*u-k_y_bar*c k_x_bar*u           u+k_x_bar*c]
             [v-k_y_bar*c k_y_bar*v           k_z_bar*v+k_x_bar*c k_x_bar*v-k_z_bar*c v+k_y_bar*c]
             [w-k_z_bar*c k_y_bar*w-k_x_bar*c k_z_bar*w           k_x_bar*w+k_y_bar*c w+k_z_bar*c]
             [H-c*Z_bar   temp1               temp2               temp3               H+c*Z_bar  ]]

    temp1 = k_y_bar*(1-b1)+(k_x_bar*w-k_z_bar*u)/c
    temp2 = k_z_bar*(1-b1)+(k_y_bar*u-k_x_bar*v)/c
    temp3 = k_x_bar*(1-b1)+(k_z_bar*v-k_y_bar*w)/c

    Rright = [[(b1+Z_bar/c)/2 -(k_x_bar/c+b2*u)/2    -(k_y_bar/c+b2*v)/2    -(k_z_bar/c+b2*w)/2    b2/2       ]
              [temp1          k_y_bar*b2*u+k_z_bar/c k_y_bar*b2*v           k_y_bar*b2*w-k_x_bar/c -k_y_bar*b2]
              [temp2          k_z_bar*b2*u-k_y_bar/c k_z_bar*b2*v+k_x_bar/c k_z_bar*b2*w           -k_z_bar*b2]
              [temp3          k_x_bar*b2*u           k_x_bar*b2*v-k_z_bar/c k_x_bar*b2*w+k_z_bar/c -k_x_bar*b2]
              [(b1-Z_bar/c)/2 (k_x_bar/c-b2*u)/2     (k_y_bar/c-b2*v)/2     (k_z_bar/c-b2*w)/2     b2/2       ]]
    
    return Rleft, Lambda, Rright
end

function RLmbdaR(R, Lam, Rm)
    Lamp = zeros(5,5)
    Lamm = zeros(5,5)
    for i in 1:length(Lam)
        Lamp[i,i] = 0.5*(Lam[i]+abs(Lam[i]))
        Lamm[i,i] = 0.5*(Lam[i]-abs(Lam[i]))
    end
    Ap = nn_inner_product(nn_inner_product(R,Lamp),Rm)
    Am = nn_inner_product(nn_inner_product(R,Lamm),Rm)
    return Ap, Am
end

function nn_inner_product(a, b)
    temp = zeros(size(a)[1],size(a)[1])
    for i in 1:size(a)[1]
        for j in 1:size(a)[1]
            for k in 1:size(a)[1]
                temp[i,j] += a[i,k]*b[k,j]
            end
        end
    end
    return temp    
end


function FVS(cell_E_hat_plas, cell_E_hat_minus, cell_F_hat_plas, cell_F_hat_minus, cell_G_hat_plas, cell_G_hat_minus, cellxmax, cellymax, cellzmax)
    E_hat = zeros(cellxmax+1,   cellymax,   cellzmax, 5)
    F_hat = zeros(  cellxmax, cellymax+1,   cellzmax, 5)
    G_hat = zeros(  cellxmax,   cellymax, cellzmax+1, 5)

    Threads.@threads for i in 2:cellxmax+1-1
        for j in 2:cellymax-1
            for k in 2:cellzmax-1
                for l in 1:5
                    E_hat[i,j,k,l] = cell_E_hat_plas[i-1,j,k,l] + cell_E_hat_minus[i,j,k,l]
                end
            end
        end
    end
    Threads.@threads for i in 2:cellxmax-1
        for j in 2:cellymax+1-1
            for k in 2:cellzmax-1
                for l in 1:5
                    F_hat[i,j,k,l] = cell_F_hat_plas[i,j-1,k,l] + cell_F_hat_minus[i,j,k,l]
                end
            end
        end
    end
    Threads.@threads for i in 2:cellxmax-1
        for j in 2:cellymax-1
            for k in 2:cellzmax+1-1
                for l in 1:5
                    G_hat[i,j,k,l] = cell_G_hat_plas[i,j,k-1,l] + cell_G_hat_minus[i,j,k,l]
                end
            end
        end
    end

    return E_hat, F_hat, G_hat
end