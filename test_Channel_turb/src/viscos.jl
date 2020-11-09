function central_diff(Qbase, Qcon, cellxmax, cellymax, cellzmax, mu, lambda, vecAx, vecAy, vecAz, specific_heat_ratio, volume, Rd)
    E_vis_hat = zeros(cellxmax+1,   cellymax,   cellzmax, 5)
    F_vis_hat = zeros(  cellxmax, cellymax+1,   cellzmax, 5)
    G_vis_hat = zeros(  cellxmax,   cellymax, cellzmax+1, 5)
    
    for i in 2:cellxmax+1 -1
        for j in 2:cellymax -1
            for k in 2:cellzmax-1

                u_av = 0.5*(Qbase[i,j,k,2]+Qbase[i-1,j,k,2])
                v_av = 0.5*(Qbase[i,j,k,3]+Qbase[i-1,j,k,3])
                w_av = 0.5*(Qbase[i,j,k,4]+Qbase[i-1,j,k,4])
                
                mu_av     = 0.5*(mu[i-1,j,k]+mu[i,j,k]) 
                lambda_av = (0.5*(1/lambda[i-1,j,k]+1/lambda[i,j,k]))^(-1)

                volume_av = 0.5*(volume[i,j,k]+volume[i-1,j,k])
        
                dudxi = Qbase[i,j,k,2] - Qbase[i-1,j,k,2]
                dvdxi = Qbase[i,j,k,3] - Qbase[i-1,j,k,3]
                dwdxi = Qbase[i,j,k,4] - Qbase[i-1,j,k,4]
                dTdxi = Qbase[i,j,k,5]/(Qbase[i,j,k,1]*Rd) - Qbase[i-1,j,k,5]/(Qbase[i-1,j,k,1]*Rd)
        
                dudeta = 0.25*(Qbase[i,j+1,k,2] - Qbase[i,j-1,k,2] + Qbase[i-1,j+1,k,2] - Qbase[i-1,j-1,k,2])
                dvdeta = 0.25*(Qbase[i,j+1,k,3] - Qbase[i,j-1,k,3] + Qbase[i-1,j+1,k,3] - Qbase[i-1,j-1,k,3])
                dwdeta = 0.25*(Qbase[i,j+1,k,4] - Qbase[i,j-1,k,4] + Qbase[i-1,j+1,k,4] - Qbase[i-1,j-1,k,4])
                dTdeta = 0.25*( Qbase[i,j+1,k,5] /   (Qbase[i,j+1,k,1]*Rd) - 
                                Qbase[i,j-1,k,5] /   (Qbase[i,j-1,k,1]*Rd) +
                                Qbase[i-1,j+1,k,5] / (Qbase[i-1,j+1,k,1]*Rd) - 
                                Qbase[i-1,j-1,k,5] / (Qbase[i-1,j-1,k,1]*Rd) )
                
                dudzeta = 0.25*(Qbase[i,j,k+1,2] - Qbase[i,j,k-1,2] + Qbase[i-1,j,k+1,2] - Qbase[i-1,j,k-1,2])
                dvdzeta = 0.25*(Qbase[i,j,k+1,3] - Qbase[i,j,k-1,3] + Qbase[i-1,j,k+1,3] - Qbase[i-1,j,k-1,3])
                dwdzeta = 0.25*(Qbase[i,j,k+1,4] - Qbase[i,j,k-1,4] + Qbase[i-1,j,k+1,4] - Qbase[i-1,j,k-1,4])
                dTdzeta = 0.25*(Qbase[i,j,k+1,5] /   (Qbase[i,j,k+1,1]*Rd) - 
                                Qbase[i,j,k-1,5] /   (Qbase[i,j,k-1,1]*Rd) +
                                Qbase[i-1,j,k+1,5] / (Qbase[i-1,j,k+1,1]*Rd) - 
                                Qbase[i-1,j,k-1,5] / (Qbase[i-1,j,k-1,1]*Rd) )
                
                vecAx_xav = vecAx[i,j,k,1]
                vecAx_yav = vecAx[i,j,k,2]
                vecAx_zav = vecAx[i,j,k,3]
                
                vecAy_xav = 0.25*( vecAy[i,j,k,1] + vecAy[i,j+1,k,1] + vecAy[i-1,j,k,1] + vecAy[i-1,j+1,k,1] )
                vecAy_yav = 0.25*( vecAy[i,j,k,2] + vecAy[i,j+1,k,2] + vecAy[i-1,j,k,2] + vecAy[i-1,j+1,k,2] )
                vecAy_zav = 0.25*( vecAy[i,j,k,3] + vecAy[i,j+1,k,3] + vecAy[i-1,j,k,3] + vecAy[i-1,j+1,k,3] )

                vecAz_xav = 0.25*( vecAz[i,j,k,1] + vecAz[i,j,k+1,1] + vecAz[i-1,j,k,1] + vecAz[i-1,j,k+1,1] )
                vecAz_yav = 0.25*( vecAz[i,j,k,2] + vecAz[i,j,k+1,2] + vecAz[i-1,j,k,2] + vecAz[i-1,j,k+1,2] )
                vecAz_zav = 0.25*( vecAz[i,j,k,3] + vecAz[i,j,k+1,3] + vecAz[i-1,j,k,3] + vecAz[i-1,j,k+1,3] )

                dudx = vecAx_xav*dudxi + vecAy_xav*dudeta +vecAz_xav*dudzeta
                dudy = vecAx_yav*dudxi + vecAy_yav*dudeta +vecAz_yav*dudzeta
                dudz = vecAx_zav*dudxi + vecAy_zav*dudeta +vecAz_zav*dudzeta
                
                dvdx = vecAx_xav*dvdxi + vecAy_xav*dvdeta +vecAz_xav*dvdzeta
                dvdy = vecAx_yav*dvdxi + vecAy_yav*dvdeta +vecAz_yav*dvdzeta
                dvdz = vecAx_zav*dvdxi + vecAy_zav*dvdeta +vecAz_zav*dvdzeta
                
                dwdx = vecAx_xav*dwdxi + vecAy_xav*dwdeta +vecAz_xav*dwdzeta
                dwdy = vecAx_yav*dwdxi + vecAy_yav*dwdeta +vecAz_yav*dwdzeta
                dwdz = vecAx_zav*dwdxi + vecAy_zav*dwdeta +vecAz_zav*dwdzeta
        

                tau_xx = 2/3*mu_av*(2*dudx - dvdy - dwdz)
                tau_yy = 2/3*mu_av*(2*dvdy - dwdz - dudx)
                tau_zz = 2/3*mu_av*(2*dwdz - dudx - dvdy)
                        
                tau_xy = mu_av*(dudy + dvdx)
                tau_yz = mu_av*(dvdz + dwdy)
                tau_zx = mu_av*(dwdx + dudz)
                
                dTdx = vecAx_xav*dTdxi + vecAy_xav*dTdeta + vecAz_xav*dTdzeta
                dTdy = vecAx_yav*dTdxi + vecAy_yav*dTdeta + vecAz_yav*dTdzeta
                dTdz = vecAx_zav*dTdxi + vecAy_zav*dTdeta + vecAz_zav*dTdzeta
                
                betax = u_av*tau_xx + v_av*tau_xy + w_av*tau_zx + lambda_av * dTdx
                betay = u_av*tau_xy + v_av*tau_yy + w_av*tau_yz + lambda_av * dTdy
                betaz = u_av*tau_zx + v_av*tau_yz + w_av*tau_zz + lambda_av * dTdz

                E_vis_hat[i,j,k,1] = 0.0
                E_vis_hat[i,j,k,2] = (vecAx_xav*tau_xx + vecAx_yav*tau_xy + vecAx_zav*tau_zx) / volume_av
                E_vis_hat[i,j,k,3] = (vecAx_xav*tau_xy + vecAx_yav*tau_yy + vecAx_zav*tau_yz) / volume_av
                E_vis_hat[i,j,k,4] = (vecAx_xav*tau_zx + vecAx_yav*tau_yz + vecAx_zav*tau_zz) / volume_av
                E_vis_hat[i,j,k,5] = (vecAx_xav*betax  + vecAx_yav*betay  + vecAx_zav*betaz )  / volume_av
            end
        end
    end
    
    for i in 2:cellxmax -1
        for j in 2:cellymax+1 -1
            for k in 2:cellzmax -1
                u_av = 0.5*(Qbase[i,j,k,2]+Qbase[i,j-1,k,2])
                v_av = 0.5*(Qbase[i,j,k,3]+Qbase[i,j-1,k,3])
                w_av = 0.5*(Qbase[i,j,k,4]+Qbase[i,j-1,k,4])

                mu_av     = 0.5*(mu[i,j-1,k]+mu[i,j,k])
                lambda_av = (0.5*(1/lambda[i,j-1,k]+1/lambda[i,j,k]))^(-1)

                volume_av = 0.5*(volume[i,j,k]+volume[i,j-1,k])
        
                dudxi = 0.25*(Qbase[i+1,j,k,2] - Qbase[i-1,j,k,2] + Qbase[i+1,j-1,k,2] - Qbase[i-1,j-1,k,2])
                dvdxi = 0.25*(Qbase[i+1,j,k,3] - Qbase[i-1,j,k,3] + Qbase[i+1,j-1,k,3] - Qbase[i-1,j-1,k,3])
                dwdxi = 0.25*(Qbase[i+1,j,k,4] - Qbase[i-1,j,k,4] + Qbase[i+1,j-1,k,4] - Qbase[i-1,j-1,k,4])
                dTdxi = 0.25*(  Qbase[i+1,j,k,4]/  (Qbase[i+1,j,k,1]*Rd)   - 
                                Qbase[i-1,j,k,4]/  (Qbase[i-1,j,k,1]*Rd)   +
                                Qbase[i+1,j-1,k,4]/(Qbase[i+1,j-1,k,1]*Rd) - 
                                Qbase[i-1,j-1,k,4]/(Qbase[i-1,j-1,k,1]*Rd) )

                dudeta = Qbase[i,j,k,2]-Qbase[i,j-1,k,2]
                dvdeta = Qbase[i,j,k,3]-Qbase[i,j-1,k,3]
                dwdeta = Qbase[i,j,k,4]-Qbase[i,j-1,k,4]
                dTdeta = (Qbase[i,j,k,5]/(Qbase[i,j,k,1]*Rd) - Qbase[i,j-1,k,5]/(Qbase[i,j-1,k,1]*Rd))

                dudzeta = 0.25*(Qbase[i,j,k+1,2] - Qbase[i,j,k-1,2] + Qbase[i,j-1,k+1,2] - Qbase[i,j-1,k-1,2])
                dvdzeta = 0.25*(Qbase[i,j,k+1,3] - Qbase[i,j,k-1,3] + Qbase[i,j-1,k+1,3] - Qbase[i,j-1,k-1,3])
                dwdzeta = 0.25*(Qbase[i,j,k+1,4] - Qbase[i,j,k-1,4] + Qbase[i,j-1,k+1,4] - Qbase[i,j-1,k-1,4])
                dTdzeta = 0.25*(Qbase[i,j,k+1,4]/  (Qbase[i,j,k+1,1]*Rd)   - 
                                Qbase[i,j,k-1,4]/  (Qbase[i,j,k-1,1]*Rd)   +
                                Qbase[i,j-1,k+1,4]/(Qbase[i,j-1,k+1,1]*Rd) - 
                                Qbase[i,j-1,k-1,4]/(Qbase[i,j-1,k-1,1]*Rd) )
                        
                vecAx_xav = 0.25*( vecAx[i,j,k,1] + vecAx[i+1,j,k,1] + vecAx[i,j-1,k,1] + vecAx[i+1,j-1,k,1] )
                vecAx_yav = 0.25*( vecAx[i,j,k,2] + vecAx[i+1,j,k,2] + vecAx[i,j-1,k,2] + vecAx[i+1,j-1,k,2] )
                vecAx_zav = 0.25*( vecAx[i,j,k,3] + vecAx[i+1,j,k,3] + vecAx[i,j-1,k,3] + vecAx[i+1,j-1,k,3] )

                vecAy_xav = vecAy[i,j,k,1]
                vecAy_yav = vecAy[i,j,k,2]
                vecAy_zav = vecAy[i,j,k,3]

                vecAz_xav = 0.25*( vecAz[i,j,k,1] + vecAz[i,j,k+1,1] + vecAz[i,j-1,k,1] + vecAz[i,j-1,k+1,1] )
                vecAz_yav = 0.25*( vecAz[i,j,k,2] + vecAz[i,j,k+1,2] + vecAz[i,j-1,k,2] + vecAz[i,j-1,k+1,2] )
                vecAz_zav = 0.25*( vecAz[i,j,k,3] + vecAz[i,j,k+1,3] + vecAz[i,j-1,k,3] + vecAz[i,j-1,k+1,3] )

                dudx = vecAx_xav*dudxi + vecAy_xav*dudeta +vecAz_xav*dudzeta
                dudy = vecAx_yav*dudxi + vecAy_yav*dudeta +vecAz_yav*dudzeta
                dudz = vecAx_zav*dudxi + vecAy_zav*dudeta +vecAz_zav*dudzeta
                
                dvdx = vecAx_xav*dvdxi + vecAy_xav*dvdeta +vecAz_xav*dvdzeta
                dvdy = vecAx_yav*dvdxi + vecAy_yav*dvdeta +vecAz_yav*dvdzeta
                dvdz = vecAx_zav*dvdxi + vecAy_zav*dvdeta +vecAz_zav*dvdzeta
                
                dwdx = vecAx_xav*dwdxi + vecAy_xav*dwdeta +vecAz_xav*dwdzeta
                dwdy = vecAx_yav*dwdxi + vecAy_yav*dwdeta +vecAz_yav*dwdzeta
                dwdz = vecAx_zav*dwdxi + vecAy_zav*dwdeta +vecAz_zav*dwdzeta
        

                tau_xx = 2/3*mu_av*(2*dudx - dvdy - dwdz)
                tau_yy = 2/3*mu_av*(2*dvdy - dwdz - dudx)
                tau_zz = 2/3*mu_av*(2*dwdz - dudx - dvdy)
                        
                tau_xy = mu_av*(dudy + dvdx)
                tau_yz = mu_av*(dvdz + dwdy)
                tau_zx = mu_av*(dwdx + dudz)
                
                dTdx = vecAx_xav*dTdxi + vecAy_xav*dTdeta + vecAz_xav*dTdzeta
                dTdy = vecAx_yav*dTdxi + vecAy_yav*dTdeta + vecAz_yav*dTdzeta
                dTdz = vecAx_zav*dTdxi + vecAy_zav*dTdeta + vecAz_zav*dTdzeta
                
                betax = u_av*tau_xx + v_av*tau_xy + w_av*tau_zx + lambda_av * dTdx
                betay = u_av*tau_xy + v_av*tau_yy + w_av*tau_yz + lambda_av * dTdy
                betaz = u_av*tau_zx + v_av*tau_yz + w_av*tau_zz + lambda_av * dTdz

                F_vis_hat[i,j,k,1] = 0.0
                F_vis_hat[i,j,k,2] = (vecAy_xav*tau_xx + vecAy_yav*tau_xy + vecAy_zav*tau_zx) / volume_av
                F_vis_hat[i,j,k,3] = (vecAy_xav*tau_xy + vecAy_yav*tau_yy + vecAy_zav*tau_yz) / volume_av
                F_vis_hat[i,j,k,4] = (vecAy_xav*tau_zx + vecAy_yav*tau_yz + vecAy_zav*tau_zz) / volume_av
                F_vis_hat[i,j,k,5] = (vecAy_xav*betax  + vecAy_yav*betay  + vecAy_zav*betaz ) / volume_av
            end
        end
    end


    for i in 2:cellxmax -1
        for j in 2:cellymax -1
            for k in 2:cellzmax+1 -1
                u_av = 0.5*(Qbase[i,j,k,2]+Qbase[i,j,k-1,2])
                v_av = 0.5*(Qbase[i,j,k,3]+Qbase[i,j,k-1,3])
                w_av = 0.5*(Qbase[i,j,k,4]+Qbase[i,j,k-1,4])

                mu_av     = 0.5*(mu[i,j,k-1]+mu[i,j,k])
                lambda_av = (0.5*(1/lambda[i,j,k-1]+1/lambda[i,j,k]))^(-1)

                volume_av = 0.5*(volume[i,j,k]+volume[i,j,k-1])
        
                dudxi = 0.25*(Qbase[i+1,j,k,2] - Qbase[i-1,j,k,2] + Qbase[i+1,j,k-1,2] - Qbase[i-1,j,k-1,2])
                dvdxi = 0.25*(Qbase[i+1,j,k,3] - Qbase[i-1,j,k,3] + Qbase[i+1,j,k-1,3] - Qbase[i-1,j,k-1,3])
                dwdxi = 0.25*(Qbase[i+1,j,k,4] - Qbase[i-1,j,k,4] + Qbase[i+1,j,k-1,4] - Qbase[i-1,j,k-1,4])
                dTdxi = 0.25*(  Qbase[i+1,j,k,4]/  (Qbase[i+1,j,k,1]*Rd)   - 
                                Qbase[i-1,j,k,4]/  (Qbase[i-1,j,k,1]*Rd)   +
                                Qbase[i+1,j,k-1,4]/(Qbase[i+1,j,k-1,1]*Rd) - 
                                Qbase[i-1,j,k-1,4]/(Qbase[i-1,j,k-1,1]*Rd) )
                
                dudeta = 0.25*(Qbase[i,j+1,k,2] - Qbase[i,j-1,k,2] + Qbase[i,j+1,k-1,2] - Qbase[i,j-1,k-1,2])
                dvdeta = 0.25*(Qbase[i,j+1,k,3] - Qbase[i,j-1,k,3] + Qbase[i,j+1,k-1,3] - Qbase[i,j-1,k-1,3])
                dwdeta = 0.25*(Qbase[i,j+1,k,4] - Qbase[i,j-1,k,4] + Qbase[i,j+1,k-1,4] - Qbase[i,j-1,k-1,4])
                dTdeta = 0.25*( Qbase[i,j+1,k,4]/  (Qbase[i,j+1,k,1]*Rd)   - 
                                Qbase[i,j-1,k,4]/  (Qbase[i,j-1,k,1]*Rd)   +
                                Qbase[i,j+1,k-1,4]/(Qbase[i,j+1,k-1,1]*Rd) - 
                                Qbase[i,j-1,k-1,4]/(Qbase[i,j-1,k-1,1]*Rd) )

                dudzeta = Qbase[i,j,k,2]-Qbase[i,j,k-1,2]
                dvdzeta = Qbase[i,j,k,3]-Qbase[i,j,k-1,3]
                dwdzeta = Qbase[i,j,k,4]-Qbase[i,j,k-1,4]
                dTdzeta = (Qbase[i,j,k,5]/(Qbase[i,j,k,1]*Rd) - Qbase[i,j,k-1,5]/(Qbase[i,j,k-1,1]*Rd))
                        
                vecAx_xav = 0.25*( vecAx[i,j,k,1] + vecAx[i+1,j,k,1] + vecAx[i,j,k-1,1] + vecAx[i+1,j,k-1,1] )
                vecAx_yav = 0.25*( vecAx[i,j,k,2] + vecAx[i+1,j,k,2] + vecAx[i,j,k-1,2] + vecAx[i+1,j,k-1,2] )
                vecAx_zav = 0.25*( vecAx[i,j,k,3] + vecAx[i+1,j,k,3] + vecAx[i,j,k-1,3] + vecAx[i+1,j,k-1,3] )

                vecAy_xav = 0.25*( vecAy[i,j,k,1] + vecAy[i,j+1,k,1] + vecAy[i,j,k-1,1] + vecAy[i,j+1,k-1,1] )
                vecAy_yav = 0.25*( vecAy[i,j,k,2] + vecAy[i,j+1,k,2] + vecAy[i,j,k-1,2] + vecAy[i,j+1,k-1,2] )
                vecAy_zav = 0.25*( vecAy[i,j,k,3] + vecAy[i,j+1,k,3] + vecAy[i,j,k-1,3] + vecAy[i,j+1,k-1,3] )

                vecAz_xav = vecAz[i,j,k,1]
                vecAz_yav = vecAz[i,j,k,2]
                vecAz_zav = vecAz[i,j,k,3]


                dudx = vecAx_xav*dudxi + vecAy_xav*dudeta +vecAz_xav*dudzeta
                dudy = vecAx_yav*dudxi + vecAy_yav*dudeta +vecAz_yav*dudzeta
                dudz = vecAx_zav*dudxi + vecAy_zav*dudeta +vecAz_zav*dudzeta
                
                dvdx = vecAx_xav*dvdxi + vecAy_xav*dvdeta +vecAz_xav*dvdzeta
                dvdy = vecAx_yav*dvdxi + vecAy_yav*dvdeta +vecAz_yav*dvdzeta
                dvdz = vecAx_zav*dvdxi + vecAy_zav*dvdeta +vecAz_zav*dvdzeta
                
                dwdx = vecAx_xav*dwdxi + vecAy_xav*dwdeta +vecAz_xav*dwdzeta
                dwdy = vecAx_yav*dwdxi + vecAy_yav*dwdeta +vecAz_yav*dwdzeta
                dwdz = vecAx_zav*dwdxi + vecAy_zav*dwdeta +vecAz_zav*dwdzeta
        

                tau_xx = 2/3*mu_av*(2*dudx - dvdy - dwdz)
                tau_yy = 2/3*mu_av*(2*dvdy - dwdz - dudx)
                tau_zz = 2/3*mu_av*(2*dwdz - dudx - dvdy)
                        
                tau_xy = mu_av*(dudy + dvdx)
                tau_yz = mu_av*(dvdz + dwdy)
                tau_zx = mu_av*(dwdx + dudz)
                
                dTdx = vecAx_xav*dTdxi + vecAy_xav*dTdeta + vecAz_xav*dTdzeta
                dTdy = vecAx_yav*dTdxi + vecAy_yav*dTdeta + vecAz_yav*dTdzeta
                dTdz = vecAx_zav*dTdxi + vecAy_zav*dTdeta + vecAz_zav*dTdzeta
                
                betax = u_av*tau_xx + v_av*tau_xy + w_av*tau_zx + lambda_av * dTdx
                betay = u_av*tau_xy + v_av*tau_yy + w_av*tau_yz + lambda_av * dTdy
                betaz = u_av*tau_zx + v_av*tau_yz + w_av*tau_zz + lambda_av * dTdz

                G_vis_hat[i,j,k,1] = 0.0
                G_vis_hat[i,j,k,2] = (vecAz_xav*tau_xx + vecAz_yav*tau_xy + vecAz_zav*tau_zx) / volume_av
                G_vis_hat[i,j,k,3] = (vecAz_xav*tau_xy + vecAz_yav*tau_yy + vecAz_zav*tau_yz) / volume_av
                G_vis_hat[i,j,k,4] = (vecAz_xav*tau_zx + vecAz_yav*tau_yz + vecAz_zav*tau_zz) / volume_av
                G_vis_hat[i,j,k,5] = (vecAz_xav*betax  + vecAz_yav*betay  + vecAz_zav*betaz ) / volume_av
            end
        end
    end


    return E_vis_hat, F_vis_hat, G_vis_hat
end
