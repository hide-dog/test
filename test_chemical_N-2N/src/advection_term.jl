function AUSM_plusup(Qbase, Qcon, cellxmax, cellymax, vecAx, vecAy, specific_heat_ratio, Minf, volume, nval)
    E_adv_hat = zeros(cellxmax+1,   cellymax, nval)
    F_adv_hat = zeros(  cellxmax, cellymax+1, nval)
    g         = specific_heat_ratio

    for i in 2:cellxmax+1 -1
        for j in 2:cellymax -1
            # i-1セル
            cell_vecAx = zeros(2)
            cell_vecAx[1] = 0.5*(vecAx[i-1,j,1]+vecAx[i,j,1])
            cell_vecAx[2] = 0.5*(vecAx[i-1,j,2]+vecAx[i,j,2])
            
            rhoL = Qbase[i-1,j,1]
            UL = Qcon[i-1,j,2] / Qcon[i-1,j,1]*cell_vecAx[1] + Qcon[i-1,j,3] / Qcon[i-1,j,1]*cell_vecAx[2]
            pL = Qbase[i-1,j,4]
            
            # iセル
            cell_vecAx = zeros(2)
            cell_vecAx[1] = 0.5*(vecAx[i,j,1]+vecAx[i+1,j,1])
            cell_vecAx[2] = 0.5*(vecAx[i,j,2]+vecAx[i+1,j,2])
            
            rhoR = Qbase[i,j,1]
            #UR = Qbase[i,j,2]*cell_vecAx[1] + Qbase[i,j,3]*cell_vecAx[2]
            UR = Qcon[i,j,2] / Qcon[i,j,1]*cell_vecAx[1] + Qcon[i,j,3] / Qcon[i,j,1]*cell_vecAx[2]
            pR = Qbase[i,j,4]
            
            # AUSM+up
            mdot, ph = AUSM_plusup_half(rhoL, rhoR, UL, UR, pL, pR, Minf, g)
#=
            if rhoL == 0
                println(" L  ")
                println(i)
                println(j)
            end
            if rhoR == 0
                println("   ")
                println(i)
                println(j)
            end
=#
            # flux half
            temp_vecX = zeros(nval)
            temp_vecX[2] = vecAx[i,j,1]
            temp_vecX[3] = vecAx[i,j,2]

            Lpsi = zeros(nval)
            Rpsi = zeros(nval)

            for l in 1:nval
                Lpsi[l] = Qcon[i-1,j,l] / Qcon[i-1,j,1]
                Rpsi[l] = Qcon[i,j,l]   / Qcon[i,j,1]
            end
            Lpsi[4] = (Qcon[i-1,j,4] + Qbase[i-1,j,4]) / Qcon[i-1,j,1]
            Rpsi[4] = (Qcon[i,j,4]   + Qbase[i,j,4])   / Qcon[i,j,1]

            if mdot > 0
                for l in 1:nval
                    E_adv_hat[i,j,l] = mdot * Lpsi[l] + ph * temp_vecX[l]
                end
            else
                for l in 1:nval
                    E_adv_hat[i,j,l] = mdot * Rpsi[l] + ph * temp_vecX[l]
                end
            end            
        end
    end

    for i in 2:cellxmax -1
        for j in 2:cellymax+1 -1
        
            # j-1セル
            cell_vecAy = zeros(2)
            cell_vecAy[1] = 0.5*(vecAy[i,j-1,1]+vecAy[i,j,1])
            cell_vecAy[2] = 0.5*(vecAy[i,j-1,2]+vecAy[i,j,2])
            
            rhoL = Qbase[i,j-1,1]
            VL = Qbase[i,j-1,2]*cell_vecAy[1] + Qbase[i,j-1,3]*cell_vecAy[2]
            pL = Qbase[i,j-1,4]
            
            # jセル
            cell_vecAy = zeros(2)
            cell_vecAy[1] = 0.5*(vecAy[i,j,1]+vecAy[i,j+1,1])
            cell_vecAy[2] = 0.5*(vecAy[i,j,2]+vecAy[i,j+1,2])
            
            rhoR = Qbase[i,j,1]
            VR = Qbase[i,j,2]*cell_vecAy[1] + Qbase[i,j,3]*cell_vecAy[2]
            pR = Qbase[i,j,4]

            # AUSM+up
            mdot, ph = AUSM_plusup_half(rhoL, rhoR, VL, VR, pL, pR, Minf, g)
            #=
            if i == 100 && j == 50
                println(" y ")
                println(mdot)
                println(ph)
            end
            =#
            
            # flux half
            temp_vecY = zeros(nval)
            temp_vecY[2] = vecAy[i,j,1]
            temp_vecY[3] = vecAy[i,j,2]

            Lpsi = zeros(nval)
            Rpsi = zeros(nval)

            for l in 1:nval
                Lpsi[l] = Qcon[i,j-1,l] / Qcon[i,j-1,1]
                Rpsi[l] = Qcon[i,j,l]   / Qcon[i,j,1]
            end
            Lpsi[4] = (Qcon[i,j-1,4] + Qbase[i,j-1,4]) / Qcon[i,j-1,1]
            Rpsi[4] = (Qcon[i,j,4]   + Qbase[i,j,4])   / Qcon[i,j,1]
            
            if mdot > 0
                for l in 1:nval
                    F_adv_hat[i,j,l] = mdot * Lpsi[l] + ph * temp_vecY[l]
                end
            else
                for l in 1:nval
                    F_adv_hat[i,j,l] = mdot * Rpsi[l] + ph * temp_vecY[l]
                end
            end
            
            #=
            if i == 35 && j ==2
                println(" y ")
                println(mdot)
                println(temp_vecY)
                println(Lpsi)
                println(Rpsi)
            end
            =#

        end
    end
    
    return E_adv_hat, F_adv_hat
end

function AUSM_plusup_half(rhoL, rhoR, UL, UR, pL, pR, Minf, g)
    # param
    beta  = 1/8
    kp    = 0.25
    sigma = 1.0
    ku    = 0.75
    
    # L, R
    aL = (g * pL / rhoL)^0.5
    aR = (g * pR / rhoR)^0.5

    # half
    ah   = 0.5*( aL + aR )
    rhoh = 0.5*( rhoL + rhoR )

    Mbar = (( UL^2 + UR^2 ) / ( 2 * ah^2 ))^0.5
    Mo   = (min(1,max( Mbar^2, Minf^2 )))^0.5
    fa   = Mo * (2-Mo)

    alpha = 3/16 * ( -4 + 5*fa )

    ML = UL/ah
    MR = UR/ah

    M_p4 = 0
    M_m4 = 0
    p_p5 = 0
    p_m5 = 0
    if abs(ML) >= 1
        M_p4 = 0.5*(ML + abs(ML))
        p_p5 = 0.5*(ML + abs(ML)) / ML
    else
        Mtp  = 0.25*(ML + 1)^2
        Mtm  = -0.25*(ML - 1)^2
        #M2m1 = (ML^2 - 1)^2
        M_p4 = Mtp * (1 - 16*beta*Mtm)
        #M_p4 = (1-ML)* (Mtp + beta*M2m1) + ML*0.5*(ML + abs(ML))
        p_p5 = Mtp * ((2-ML) - 16*alpha*ML*Mtm)
    end
    

    if abs(MR) >= 1
        M_m4 = 0.5*(MR - abs(MR))
        p_m5 = 0.5*(MR - abs(MR)) / MR
    else
        Mtp  = 0.25*(MR + 1)^2
        Mtm  = -0.25*(MR - 1)^2
        M_m4 = Mtm * (1 + 16*beta*Mtp)
        p_m5 = Mtm * ((-2-MR) + 16*alpha*MR*Mtp)
    end
    
    #=
    Mtp  = 0.25*(ML + 1)^2
    Mtm  = -0.25*(ML - 1)^2
    M2m1 = (ML^2 - 1)^2
    M_p4 = (1-ML)* (Mtp + beta*M2m1) + ML*0.5*(ML + abs(ML))
    
    Mtp  = 0.25*(MR + 1)^2
    Mtm  = -0.25*(MR - 1)^2
    M2m1 = (MR^2 - 1)^2
    M_m4 = (1-MR)* (Mtm - beta*M2m1) + MR*0.5*(MR - abs(MR))
    =#

    # M half
    Mp = -kp * max( 1 - sigma*Mbar^2, 0) * (pR - pL)/(rhoh*ah^2)

    # Mh = M_p4 + M_m4 + Mp/fa
    Mh = M_p4 + M_m4
    
    # mdot half
    mdot = ah * Mh
    if Mh > 0
        mdot = mdot * rhoL
    else
        mdot = mdot * rhoR
    end
    
    # p half
    pu = -ku * p_p5 * p_m5 * (rhoL + rhoR) * ah *(UR - UL)

    ph = p_p5*pL + p_m5*pR + fa * pu
    return mdot, ph, ah, Mh
end

function AUSM_plus(Qbase, Qcon, cellxmax, cellymax, vecAx, vecAy, specific_heat_ratio, volume, nval)
    E_adv_hat = zeros(cellxmax+1,   cellymax, nval)
    F_adv_hat = zeros(  cellxmax, cellymax+1, nval)
    g         = specific_heat_ratio

    for i in 2:cellxmax+1 -1
        for j in 2:cellymax -1
            # i-1セル
            cell_vecAx = zeros(2)
            cell_vecAx[1] = 0.5*(vecAx[i-1,j,1]+vecAx[i,j,1])
            cell_vecAx[2] = 0.5*(vecAx[i-1,j,2]+vecAx[i,j,2])
            
            rhoL = Qbase[i-1,j,1]
            UL = Qcon[i-1,j,2] / Qcon[i-1,j,1]*cell_vecAx[1] + Qcon[i-1,j,3] / Qcon[i-1,j,1]*cell_vecAx[2]
            pL = Qbase[i-1,j,4]
            
            # iセル
            cell_vecAx = zeros(2)
            cell_vecAx[1] = 0.5*(vecAx[i,j,1]+vecAx[i+1,j,1])
            cell_vecAx[2] = 0.5*(vecAx[i,j,2]+vecAx[i+1,j,2])
            
            rhoR = Qbase[i,j,1]
            #UR = Qbase[i,j,2]*cell_vecAx[1] + Qbase[i,j,3]*cell_vecAx[2]
            UR = Qcon[i,j,2] / Qcon[i,j,1]*cell_vecAx[1] + Qcon[i,j,3] / Qcon[i,j,1]*cell_vecAx[2]
            pR = Qbase[i,j,4]
            
            # AUSM+
            mdot, ph = AUSM_plus_half(rhoL, rhoR, UL, UR, pL, pR, g)

            # flux half
            temp_vecX = zeros(nval)
            temp_vecX[2] = vecAx[i,j,1]
            temp_vecX[3] = vecAx[i,j,2]

            Lpsi = zeros(nval)
            Rpsi = zeros(nval)

            for l in 1:nval
                Lpsi[l] = Qcon[i-1,j,l] / Qcon[i-1,j,1]
                Rpsi[l] = Qcon[i,j,l]   / Qcon[i,j,1]
            end
            Lpsi[4] = (Qcon[i-1,j,4] + Qbase[i-1,j,4]) / Qcon[i-1,j,1]
            Rpsi[4] = (Qcon[i,j,4]   + Qbase[i,j,4])   / Qcon[i,j,1]

            if mdot > 0
                for l in 1:nval
                    E_adv_hat[i,j,l] = mdot * Lpsi[l] + ph * temp_vecX[l]
                end
            else
                for l in 1:nval
                    E_adv_hat[i,j,l] = mdot * Rpsi[l] + ph * temp_vecX[l]
                end
            end            
        end
    end

    for i in 2:cellxmax -1
        for j in 2:cellymax+1 -1
        
            # j-1セル
            cell_vecAy = zeros(2)
            cell_vecAy[1] = 0.5*(vecAy[i,j-1,1]+vecAy[i,j,1])
            cell_vecAy[2] = 0.5*(vecAy[i,j-1,2]+vecAy[i,j,2])
            
            rhoL = Qbase[i,j-1,1]
            VL = Qbase[i,j-1,2]*cell_vecAy[1] + Qbase[i,j-1,3]*cell_vecAy[2]
            pL = Qbase[i,j-1,4]
            
            # jセル
            cell_vecAy = zeros(2)
            cell_vecAy[1] = 0.5*(vecAy[i,j,1]+vecAy[i,j+1,1])
            cell_vecAy[2] = 0.5*(vecAy[i,j,2]+vecAy[i,j+1,2])
            
            rhoR = Qbase[i,j,1]
            VR = Qbase[i,j,2]*cell_vecAy[1] + Qbase[i,j,3]*cell_vecAy[2]
            pR = Qbase[i,j,4]

            # AUSM+up
            mdot, ph = AUSM_plus_half(rhoL, rhoR, VL, VR, pL, pR, g)
            
            # flux half
            temp_vecY = zeros(nval)
            temp_vecY[2] = vecAy[i,j,1]
            temp_vecY[3] = vecAy[i,j,2]

            Lpsi = zeros(nval)
            Rpsi = zeros(nval)

            for l in 1:nval
                Lpsi[l] = Qcon[i,j-1,l] / Qcon[i,j-1,1]
                Rpsi[l] = Qcon[i,j,l]   / Qcon[i,j,1]
            end
            Lpsi[4] = (Qcon[i,j-1,4] + Qbase[i,j-1,4]) / Qcon[i,j-1,1]
            Rpsi[4] = (Qcon[i,j,4]   + Qbase[i,j,4])   / Qcon[i,j,1]
            
            if mdot > 0
                for l in 1:nval
                    F_adv_hat[i,j,l] = mdot * Lpsi[l] + ph * temp_vecY[l]
                end
            else
                for l in 1:nval
                    F_adv_hat[i,j,l] = mdot * Rpsi[l] + ph * temp_vecY[l]
                end
            end
        end
    end
    
    return E_adv_hat, F_adv_hat
end

function AUSM_plus_half(rhoL, rhoR, UL, UR, pL, pR, g)
    # param
    beta  = 1/8
    alpha = 3/16
    
    # L, R
    aL = (g * pL / rhoL)^0.5
    aR = (g * pR / rhoR)^0.5

    # half
    ah   = 0.5*( aL + aR )
    rhoh = 0.5*( rhoL + rhoR )

    ML = UL/ah
    MR = UR/ah

    M_p4 = 0
    M_m4 = 0
    p_p5 = 0
    p_m5 = 0
    if abs(ML) >= 1
        M_p4 = 0.5*(ML + abs(ML))
        p_p5 = 0.5*(ML + abs(ML)) / ML
    else
        M_p4 = 0.5*(ML+1)^2 + beta*(ML^2-1)^2
        p_p5 = 0.25*(ML+1)^2 * (2-ML) + alpha*ML*(ML^2-1)^2
    end
    
    if abs(MR) >= 1
        M_m4 = 0.5*(MR - abs(MR))
        p_m5 = 0.5*(MR - abs(MR)) / MR
    else
        M_m4 = -0.5*(MR-1)^2 - beta*(MR^2-1)^2
        p_m5 = 0.25*(MR-1)^2 * (2+MR) - alpha*MR*(MR^2-1)^2
    end
    
    # Mh = M_p4 + M_m4 + Mp/fa
    Mh = M_p4 + M_m4
    
    # mdot half
    mdot = ah * Mh
    if Mh > 0
        mdot = mdot * rhoL
    else
        mdot = mdot * rhoR
    end
    
    # p half
    ph = p_p5*pL + p_m5*pR

    return mdot, ph
end