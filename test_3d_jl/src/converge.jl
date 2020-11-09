function set_res(Delta_Qcon_hat, Delta_Qcon_hat_temp, cellxmax, cellymax, cellzmax)
    res = zeros(cellxmax, cellymax, cellzmax, 5)
    
    for i in 2:cellxmax-1
        for j in 2:cellymax-1
            for k in 2:cellzmax-1
                for l in 1:5
                    res[i,j,k,l] = abs(Delta_Qcon_hat[i,j,k,l]-Delta_Qcon_hat_temp[i,j,k,l])
                end
            end
        end
    end 

    # println(res[3,3,:])
    # println(dt*RHS[3,3,:])
    return res
end

function check_converge(res, RHS, cellxmax, cellymax,cellzmax, init_small)
    norm2 = zeros(5)

    tempAxb = zeros(5)
    tempb = zeros(5)
    for i in 2:cellxmax-1
        for j in 2:cellymax-1
            for k in 2:cellzmax-1
                for l in 1:5
                    tempAxb[l] = tempAxb[l] + res[i,j,k,l]^2
                    tempb[l] = tempb[l] + RHS[i,j,k,l]^2
                end
            end
        end
    end
    
    
    #println(tempb)
    #println(tempAxb)
    #throw(UndefVarError(:x))    

    for l in 1:5
        norm2[l] = (tempAxb[l]/(tempb[l]+init_small)) ^ 0.5
    end
    return norm2
end