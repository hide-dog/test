function base_to_conservative(Qbase, cellxmax, cellymax, cellzmax, specific_heat_ratio)
    """
    Qbase=[rho,u,v,w,p]
    Qcon=[rho,rhou,rhov,rhow,e]
    """
    Qcon = zeros(cellxmax, cellymax, cellzmax, 5)
    Threads.@threads for i in 1:cellxmax
        for j in 1:cellymax
            for k in 1:cellzmax
                Qcon[i,j,k,1] = Qbase[i,j,k,1]
                Qcon[i,j,k,2] = Qbase[i,j,k,1]*Qbase[i,j,k,2]
                Qcon[i,j,k,3] = Qbase[i,j,k,1]*Qbase[i,j,k,3]
                Qcon[i,j,k,4] = Qbase[i,j,k,1]*Qbase[i,j,k,4]
                Qcon[i,j,k,5] = Qbase[i,j,k,5]/(specific_heat_ratio-1) + 
                                Qbase[i,j,k,1]*(Qbase[i,j,k,2]^2+Qbase[i,j,k,3]^2+Qbase[i,j,k,4]^2)/2
            end
        end
    end
    return Qcon
end

function conservative_to_base(Qcon, cellxmax, cellymax, cellzmax, specific_heat_ratio)
    Qbase = zeros(cellxmax, cellymax, cellzmax, 5)
    Threads.@threads for i in 1:cellxmax
        for j in 1:cellymax
            for k in 1:cellzmax
                Qbase[i,j,k,1] = Qcon[i,j,k,1]
                Qbase[i,j,k,2] = Qcon[i,j,k,2]/Qcon[i,j,k,1]
                Qbase[i,j,k,3] = Qcon[i,j,k,3]/Qcon[i,j,k,1]
                Qbase[i,j,k,4] = Qcon[i,j,k,4]/Qcon[i,j,k,1]
                Qbase[i,j,k,5] = (Qcon[i,j,k,5]-Qcon[i,j,k,1]*(Qbase[i,j,k,2]^2+Qbase[i,j,k,3]^2+Qbase[i,j,k,4]^2)/2)*
                                    (specific_heat_ratio-1)
            end
        end
    end
    return Qbase
end

function Q_to_flux(Qcon, Qbase)
    throw(UndefVarError(:x))
    E=ones(size(Qcon)[1],4)
    for i in 1:size(E)[1]
        E[i,1]=Qcon[i,2]
        E[i,2]=Qcon[i,2]*Qbase[i,2]+Qbase[i,4]
        E[i,3]=Qcon[i,2]*Qbase[i,3]
        E[i,4]=(Qcon[i,4]+Qbase[i,4])*Qbase[i,2]
    end

    F=ones(size(Qcon)[1],4)
    for i in 1:size(F)[1]
        F[i,1]=Qcon[i,3]
        F[i,2]=Qcon[i,3]*Qbase[i,2]
        F[i,3]=Qcon[i,3]*Qbase[i,3]+Qbase[i,4]
        F[i,4]=(Qcon[i,4]+Qbase[i,4])*Qbase[i,3]
    end
    return E,F
end

function setup_Qcon_hat(Qcon, cellxmax, cellymax, cellzmax, volume)
    Qcon_hat = zeros(cellxmax, cellymax, cellzmax, 5)
    Threads.@threads for i in 1:cellxmax
        for j in 1:cellymax
            for k in 1:cellzmax
                for l in 1:5
                    Qcon_hat[i,j,k,l] = Qcon[i,j,k,l]*volume[i,j,k]
                end
            end
        end
    end
    return Qcon_hat
end

function Qhat_to_Q(Qcon_hat, cellxmax, cellymax, cellzmax, volume)
    Qcon = zeros(cellxmax, cellymax, cellzmax, 5)
    Threads.@threads for i in 1:cellxmax
        for j in 1:cellymax
            for k in 1:cellzmax
                for l in 1:5
                    Qcon[i,j,k,l] = Qcon_hat[i,j,k,l]/volume[i,j,k]
                end
            end
        end
    end
    return Qcon
end
