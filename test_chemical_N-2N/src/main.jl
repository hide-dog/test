using ProgressMeter
using Dates

function main()
    start_t = now()

    out_dir  = "result"
    PARAMDAT = "PARAMDAT.json"
    fwrite   = "write"
    
    nval = 6    # rho,u,v,p, N2, N
    nch  = 2    # N2, N
    nre  = 2    # N2+N2 -> 2N+N2
                # N2+N  -> 2N+N
    R    = set_gas_const()

    xmax, ymax, nodes, vecAx, vecAy = read_allgrid()
    out_file_front, out_ext, restartnum, restart_file, init_small, norm_ok,
    time_integ, nt, dt, every_outnum, in_nt, dtau, cfl,
    init_rho, init_u, init_v, init_p, init_N2, init_N, specific_heat_ratio, Rd, bdcon = input_para(PARAMDAT)

    Qbase, cellxmax, cellymax, restartnum, init_step, init_time = 
        set_initQbase(xmax, ymax, restart_file, init_rho, init_u, init_v, init_p, init_N2, init_N, dt,
        specific_heat_ratio, out_file_front, out_ext, out_dir, restartnum, R, nch, nval)

    # init Delta_Qcon_hat
    volume = set_volume(nodes, cellxmax, cellymax)
    Minf   = set_Minf(bdcon, specific_heat_ratio)
    dx, dy = set_dx_lts(nodes, cellxmax, cellymax)
    reset_write(fwrite)
    
    # main loop
    print("threads num : ")
    println(Threads.nthreads())

    # check bd
    check_bd(bdcon)
    
    #throw(UndefVarError(:x))

    # main loop
    prog = Progress(nt,1)
    for t in 1:nt
        next!(prog)
                
        evalnum = t+restartnum
        if time_integ == "1"
            Rhat = set_gasconst_hat(Qbase,cellxmax,cellymax,nval,nch,R)

            # exlicit scheme
            Qbase    = set_boundary(Qbase, cellxmax, cellymax, vecAx, vecAy, bdcon, Rhat, nval, nch)
            Qcon     = base_to_conservative(Qbase, cellxmax, cellymax, specific_heat_ratio, nval, nch)
            Qcon_hat = setup_Qcon_hat(Qcon, cellxmax, cellymax, volume, nval, nch)
            
            # initial_setup
            # 境界条件分のRhatを定義
            Rhat = set_gasconst_hat(Qbase,cellxmax,cellymax,nval,nch,R)
            chmu, chlambda_tr, chD, chmolef = chemical_value(Qbase, cellxmax, cellymax, Rhat, nval, nch)
            
            # RHS
            # source term
            wdot = set_wdot(Qbase, cellxmax, cellymax, Rhat, nch, nre, nval)
            
            # advection_term           
            #E_adv_hat, F_adv_hat = AUSM_plusup(Qbase, Qcon, cellxmax, cellymax, vecAx, vecAy, specific_heat_ratio, Minf, volume, nval)
            E_adv_hat, F_adv_hat = AUSM_plus(Qbase, Qcon, cellxmax, cellymax, vecAx, vecAy, specific_heat_ratio, volume, nval)
            
            # viscos_term
            E_vis_hat, F_vis_hat = central_diff(Qbase, Qcon, cellxmax, cellymax, chmu, chlambda_tr, chD, chmolef,
                                                vecAx, vecAy, specific_heat_ratio, volume, Rhat, nval, nch)                
            
            
            
            RHS = setup_RHS(cellxmax, cellymax, E_adv_hat, F_adv_hat, E_vis_hat, F_vis_hat, wdot, nval, volume)
            
            println(wdot[50,10,:])

            # time integral
            Qcon_hat = time_integration_explicit(dt, Qcon_hat, RHS, cellxmax, cellymax, nval)

            Qcon  = Qhat_to_Q(Qcon_hat, cellxmax, cellymax, volume, nval, nch)

            println(Qbase[50,10,:])
            Qbase = conservative_to_base(Qcon, cellxmax, cellymax, specific_heat_ratio, nval, nch)
            println(Qbase[50,10,:])

        elseif time_integ == "2"

            Rhat  = set_gasconst_hat(Qbase,cellxmax,cellymax,nval,nch,R)
            Qbase = set_boundary(Qbase, cellxmax, cellymax, vecAx, vecAy, bdcon, Rhat, nval, nch)

            # Qcn
            Qbasen    = copy(Qbase)
            Qconn     = base_to_conservative(Qbasen, cellxmax, cellymax, specific_heat_ratio, nval, nch)
            Qconn_hat = setup_Qcon_hat(Qconn, cellxmax, cellymax, volume, nval, nch)
            Qbasem    = copy(Qbase)

            for tau in 1:in_nt
                # LHS (A_adv_hat=jacobian)
                Rhat = set_gasconst_hat(Qbasem, cellxmax, cellymax, nval, nch, R)

                Qbasem   = set_boundary(Qbasem, cellxmax, cellymax, vecAx, vecAy, bdcon, Rhat, nval, nch)
                Qcon     = base_to_conservative(Qbasem, cellxmax, cellymax, specific_heat_ratio, nval, nch)
                Qcon_hat = setup_Qcon_hat(Qcon, cellxmax, cellymax, volume, nval, nch)

                # initial_setup
                # 境界条件分のRhatを定義
                Rhat = set_gasconst_hat(Qbasem, cellxmax, cellymax, nval, nch, R)
                chmu, chlambda_tr, chD, chmolef = chemical_value(Qbasem, cellxmax, cellymax, Rhat, nval, nch)
                dtau   = set_lts(Qbasem, cellxmax, cellymax, chmu, dx, dy, vecAx, vecAy, volume, specific_heat_ratio, cfl)

                # RHS
                # sorce term
                wdot, H_hat = set_wdot_for_implicit(Qbasem, cellxmax, cellymax, Rhat, nch, nre, nval, volume)
                println(wdot[50,10,:])

                # advection_term
                #E_adv_hat, F_adv_hat = AUSM_plusup(Qbasem, Qcon, cellxmax, cellymax, vecAx, vecAy, specific_heat_ratio, Minf, volume, nval, nch)
                E_adv_hat, F_adv_hat = AUSM_plus(Qbasem, Qcon, cellxmax, cellymax, vecAx, vecAy, specific_heat_ratio, volume, nval)

                # viscos_term
                E_vis_hat, F_vis_hat = central_diff(Qbasem, Qcon, cellxmax, cellymax, chmu, chlambda_tr, chD, chmolef,
                                                vecAx, vecAy, specific_heat_ratio, volume, Rhat, nval, nch)

                RHS = setup_RHS(cellxmax, cellymax, E_adv_hat, F_adv_hat, E_vis_hat, F_vis_hat, wdot, nval, volume)

                # lusgs_advection_term
                A_adv_hat_p, A_adv_hat_m, B_adv_hat_p, B_adv_hat_m, A_beta_shig, B_beta_shig = 
                    one_wave(Qbasem, Qcon, cellxmax, cellymax, vecAx, vecAy, specific_heat_ratio, volume, nval, nch)
                # lusgs_viscos_term
                jalphaP, jbetaP = 
                    central_diff_jacobian(Qbasem, Qcon, cellxmax, cellymax, chmu, chlambda_tr, vecAx, vecAy, specific_heat_ratio, volume)
                
                # LUSGS
                delta_Q = zeros(cellxmax, cellymax, nval)

                ite=0
                norm2 = zeros(nval)
                while true
                    delta_Q_temp = copy(delta_Q)
                    delta_Q      = lusgs_point(dt, dtau, Qcon_hat, Qconn_hat, delta_Q, 
                                                A_adv_hat_p, A_adv_hat_m, B_adv_hat_p, B_adv_hat_m,
                                                A_beta_shig, B_beta_shig, jalphaP, jbetaP, 
                                                H_hat, RHS, cellxmax, cellymax, volume, nval)
                    
                    res   = set_res(delta_Q, delta_Q_temp, cellxmax, cellymax, nval)
                    norm2 = check_converge(res, RHS, cellxmax, cellymax, init_small, nval)

                    ite += 1
                    if norm2[1] < norm_ok && norm2[2] < norm_ok && norm2[3] < norm_ok && norm2[4] < norm_ok
                        break
                    end
                    if ite % 100 ==0
                        println(" now cal norm2 ")
                        println(norm2)
                    end
                end
                output_innertime(fwrite, tau, norm2, nval)

                #println(delta_Q[50,10,:])

                for i in 2:cellxmax-1
                    for j in 2:cellymax-1
                        for l in 1:nval
                            Qcon_hat[i,j,l] = Qcon_hat[i,j,l] + delta_Q[i,j,l]
                        end
                    end
                end
                
                Qcon   = Qhat_to_Q(Qcon_hat, cellxmax, cellymax, volume, nval, nch)
                Qbasem = conservative_to_base(Qcon, cellxmax, cellymax, specific_heat_ratio, nval, nch)
            end
            println(Qbase[50,10,:])
            Qbase = copy(Qbasem)
            println(Qbase[50,10,:])
        end
        #throw(UndefVarError(:x))
        
        output_physicaltime(fwrite, t, dt, init_step, init_time)

        check_divrege(Qbase, cellxmax, cellymax, Rhat, fwrite)
        check_small(Qbase, cellxmax, cellymax, nch)
        
        if round(evalnum) % every_outnum == 0
            println("\n")
            println("nt_______________________________"*string(round(evalnum)))
            Rhat = set_gasconst_hat(Qbase,cellxmax,cellymax,nval,nch,R)
            output_result(evalnum, Qbase, cellxmax, cellymax, specific_heat_ratio, out_file_front, out_ext, out_dir, Rhat, nval)
        end
    end
    end_t = now()
    output_fin(fwrite, start_t, end_t, nt, dt, in_nt, cellxmax, cellymax)
end

# -- main --
main()