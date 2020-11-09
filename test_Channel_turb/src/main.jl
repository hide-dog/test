using ProgressMeter

function main()
    out_dir="result"
    PARAMDAT="PARAMDAT.json"

    xmax, ymax, zmax, nodes, vecAx, vecAy, vecAz = read_allgrid()

    out_file_front, out_ext, restartnum, restart_file, init_small, norm_ok,
    time_integ, nt, dt, every_outnum, in_nt, dtau,
    init_rho, init_u, init_v, init_w, init_p, specific_heat_ratio, Rd, bdcon = input_para(PARAMDAT)

    Qbase, cellxmax, cellymax, cellzmax, restartnum = set_initQbase(xmax, ymax, zmax, restart_file, init_rho, 
                                                                    init_u, init_v, init_w, init_p,
                                                                    specific_heat_ratio, out_file_front,
                                                                    out_ext, out_dir, restartnum, Rd)
    
    # init Delta_Qcon_hat
    volume         = set_volume(nodes, vecAx, vecAy, vecAz, cellxmax, cellymax, cellzmax)
    cellcenter     = set_cellcenter(nodes, cellxmax, cellymax, cellzmax)

    wallpoint, nop = set_wallpoint(nodes, bdcon, cellxmax, cellymax, cellzmax)
    wally          = set_wally(cellcenter, wallpoint, nop, cellxmax, cellymax, cellzmax)

    # main loop
    print("threads num : ")
    println(Threads.nthreads())
    
    #throw(UndefVarError(:x))
    # main loop
    prog = Progress(nt,1)
    for t in 1:nt
        next!(prog)
        
        evalnum = t+restartnum
        if time_integ == "1"
            # exlicit scheme
            Qbase    = set_boundary(Qbase, cellxmax, cellymax, cellzmax, vecAx, vecAy, vecAz, bdcon)
            Qcon     = base_to_conservative(Qbase, cellxmax, cellymax, cellzmax, specific_heat_ratio)
            Qcon_hat = setup_Qcon_hat(Qcon, cellxmax, cellymax, cellzmax, volume)            

            # initial_setup
            mu     = set_mu(Qbase, cellxmax, cellymax, cellzmax, specific_heat_ratio, Rd)
            lambda = set_lambda(Qbase, cellxmax, cellymax, cellzmax, mu, specific_heat_ratio, Rd)
            
            # RHS
            # advection_term
            # E_adv_hat, F_adv_hat = AUSM(Qbase,Qcon,cellxmax, cellymax, cellzmax,vecAx,vecAy,specific_heat_ratio)
            cell_Ahat_plas, cell_Ahat_minus, cell_Bhat_plas, cell_Bhat_minus, cell_Chat_plas, cell_Chat_minus = 
                cal_jacobi(Qbase, Qcon, cellxmax, cellymax, cellzmax, specific_heat_ratio, vecAx, vecAy, vecAz, volume)
            cell_E_hat_plas, cell_E_hat_minus, cell_F_hat_plas, cell_F_hat_minus, cell_G_hat_plas, cell_G_hat_minus =
                setup_cell_flux_hat(Qcon, cellxmax, cellymax, cellzmax, cell_Ahat_plas, cell_Ahat_minus, cell_Bhat_plas, cell_Bhat_minus, cell_Chat_plas, cell_Chat_minus)
            E_adv_hat, F_adv_hat, G_adv_hat =
                FVS(cell_E_hat_plas, cell_E_hat_minus, cell_F_hat_plas, cell_F_hat_minus, cell_G_hat_plas, cell_G_hat_minus, cellxmax, cellymax, cellzmax)

            # viscos_term
            E_vis_hat, F_vis_hat, G_vis_hat = central_diff(Qbase, Qcon, cellxmax, cellymax, cellzmax, mu, lambda, vecAx, vecAy, vecAz, specific_heat_ratio, volume, Rd)
            #=
            E_vis_hat = zeros(cellxmax+1, cellymax+1, cellzmax+1, 5)
            F_vis_hat = zeros(cellxmax+1, cellymax+1, cellzmax+1, 5)
            G_vis_hat = zeros(cellxmax+1, cellymax+1, cellzmax+1, 5)
            =#

            # turblence_term
            yplus = cal_yplus(Qbase, wally, mu, cellxmax, cellymax, cellzmax)
            
            # RHS
            RHS = setup_RHS(cellxmax, cellymax, cellzmax, E_adv_hat, F_adv_hat, G_adv_hat, E_vis_hat, F_vis_hat, G_vis_hat)
            
            # time integral
            Qcon_hat = time_integration_explicit(dt, Qcon_hat, RHS, cellxmax, cellymax, cellzmax)

            # reset
            Qcon  = Qhat_to_Q(Qcon_hat, cellxmax, cellymax, cellzmax, volume)
            Qbase = conservative_to_base(Qcon, cellxmax, cellymax, cellzmax, specific_heat_ratio)

        elseif time_integ == "2"
            # Qcn
            Qbasen    = copy(Qbase)
            Qconn     = base_to_conservative(Qbasen, cellxmax, cellymax, cellzmax, specific_heat_ratio)
            Qconn_hat = setup_Qcon_hat(Qconn, cellxmax, cellymax, cellzmax, volume)

            Qbasem = copy(Qbase)

            for tau in 1:in_nt
                # LHS (A_adv_hat=jacobian)
                Qbasem   = set_boundary(Qbasem, cellxmax, cellymax, cellzmax, vecAx, vecAy, vecAz, bdcon)  
                Qcon     = base_to_conservative(Qbasem, cellxmax, cellymax, cellzmax, specific_heat_ratio)
                Qcon_hat = setup_Qcon_hat(Qcon, cellxmax, cellymax, cellzmax, volume)

                # initial_setup
                mu     = set_mu(Qbasem, cellxmax, cellymax, cellzmax, specific_heat_ratio, Rd)
                lambda = set_lambda(Qbasem, cellxmax, cellymax, cellzmax, mu, specific_heat_ratio, Rd)

                # turb
                # tauwall = set_tauwall(Qbasem,wally,mu,cellxmax, cellymax, cellzmax)
                                
                # RHS
                #advection_term
                #E_adv_hat, F_adv_hat = AUSM(Qbasem,Qcon,cellxmax, cellymax, cellzmax,vecAx,vecAy,specific_heat_ratio)
                cell_Ahat_plas, cell_Ahat_minus, cell_Bhat_plas, cell_Bhat_minus, cell_Chat_plas, cell_Chat_minus = 
                    cal_jacobi(Qbasem, Qcon, cellxmax, cellymax, cellzmax, specific_heat_ratio, vecAx, vecAy, vecAz, volume)
                cell_E_hat_plas, cell_E_hat_minus, cell_F_hat_plas, cell_F_hat_minus, cell_G_hat_plas, cell_G_hat_minus =
                    setup_cell_flux_hat(Qcon, cellxmax, cellymax, cellzmax, cell_Ahat_plas, cell_Ahat_minus, cell_Bhat_plas, cell_Bhat_minus, cell_Chat_plas, cell_Chat_minus)
                E_adv_hat, F_adv_hat, G_adv_hat =
                    FVS(cell_E_hat_plas, cell_E_hat_minus, cell_F_hat_plas, cell_F_hat_minus, cell_G_hat_plas, cell_G_hat_minus, cellxmax, cellymax, cellzmax)


                # viscos_term
                E_vis_hat, F_vis_hat, G_vis_hat = central_diff(Qbasem, Qcon, cellxmax, cellymax, cellzmax, mu, lambda, vecAx, vecAy, vecAz, specific_heat_ratio, volume, Rd)
                
                RHS = setup_RHS(cellxmax, cellymax, cellzmax, E_adv_hat, F_adv_hat, G_adv_hat, E_vis_hat, F_vis_hat, G_vis_hat)

                #println(E_adv_hat[:,50,:])
                #println(E_adv_hat[:,50,:])

                #println(RHS[:,50,:])

                #throw(UndefVarError(:x))


                # lusgs_advection_term
                A_adv_hat_p, A_adv_hat_m, B_adv_hat_p, B_adv_hat_m, C_adv_hat_p, C_adv_hat_m, A_beta_shig, B_beta_shig, C_beta_shig = 
                    one_wave(Qbase, Qcon, cellxmax, cellymax, cellzmax, vecAx, vecAy, vecAz, specific_heat_ratio, volume)
                # lusgs_viscos_term
                jalphaP, jbetaP, jgammaP = 
                    central_diff_jacobian(Qbasem, Qcon, cellxmax, cellymax, cellzmax, mu, lambda, vecAx, vecAy, vecAz, specific_heat_ratio, volume)
                #jalphaP = zeros(cellxmax, cellymax)
                #jbetaP = zeros(cellxmax, cellymax)

                # LUSGS
                delta_Q = zeros(cellxmax, cellymax, cellzmax,5)

                ite=0
                while true
                    delta_Q_temp = copy(delta_Q)
                    delta_Q = lusgs(dt, Qcon_hat, Qconn_hat, delta_Q, A_adv_hat_p, A_adv_hat_m, B_adv_hat_p, B_adv_hat_m, C_adv_hat_p, C_adv_hat_m,
                                    A_beta_shig, B_beta_shig, C_beta_shig, jalphaP, jbetaP, jgammaP, RHS, cellxmax, cellymax,cellzmax, volume)
                    
                    #println(RHS[2,:,:])
                    #println(delta_Q[2,:,:])

                    res   = set_res(delta_Q, delta_Q_temp, cellxmax, cellymax, cellzmax)
                    norm2 = check_converge(res, RHS, cellxmax, cellymax,cellzmax, init_small)

                    ite += 1
                    if norm2[1] < norm_ok && norm2[2] < norm_ok && norm2[3] < norm_ok && norm2[4] < norm_ok && norm2[5] < norm_ok
                        break
                    end
                    if ite % 100 ==0
                        println(norm2)
                    end
                end
                #println(ite)
                #println(delta_Q[:,50,:])
                
                for i in 2:cellxmax-1
                    for j in 2:cellymax-1
                        for k in 2:cellzmax-1
                            for l in 1:5
                                Qcon_hat[i,j,k,l] = Qcon_hat[i,j,k,l] + delta_Q[i,j,k,l]
                            end
                        end
                    end
                end

                Qcon   = Qhat_to_Q(Qcon_hat, cellxmax, cellymax, cellzmax, volume)
                Qbasem = conservative_to_base(Qcon, cellxmax, cellymax, cellzmax, specific_heat_ratio)

                #println("xxx")
            end            
            Qbase = copy(Qbasem)
        end
        
        if round(evalnum) % every_outnum == 0
            println("nt_______________________________"*string(round(evalnum)))
            output_result(evalnum, Qbase, cellxmax, cellymax, cellzmax, out_file_front, out_ext, out_dir, Rd, yplus)
        end

        if isequal(Qbase[2,2,2,1], NaN) == true
            println("\n NaN !! \n")
            println(" check !!")

            throw(UndefVarError(:x))
            break
        end
        if isequal(yplus[2,2,2], NaN) == true
            println("\n NaN yplus!! \n")
            println(" check yplus!!")

            throw(UndefVarError(:x))
            break
        end

        for i in 2:cellxmax-1
            for j in 2:cellymax-1
                for k in 2:cellzmax-1
                    if Qbase[i,j,k,2] < 1.0e-12
                        Qbase[i,j,k,2] = 0.0
                    end
                end
            end
        end
    end
end

# -- main --
main()