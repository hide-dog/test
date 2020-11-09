function setup_RHS(cellxmax, cellymax, cellzmax, E_adv_hat, F_adv_hat, G_adv_hat, E_vis_hat, F_vis_hat, G_vis_hat)
    RHS = zeros(cellxmax, cellymax, cellzmax, 5)
    for i in 2:cellxmax-1
        for j in 2:cellymax-1
            for k in 2:cellzmax-1
                for l in 1:5
                    RHS[i,j,k,l] = -((E_adv_hat[i+1,j,k,l]-E_vis_hat[i+1,j,k,l]) - (E_adv_hat[i,j,k,l]-E_vis_hat[i,j,k,l]) +
                                     (F_adv_hat[i,j+1,k,l]-F_vis_hat[i,j+1,k,l]) - (F_adv_hat[i,j,k,l]-F_vis_hat[i,j,k,l]) +
                                     (G_adv_hat[i,j,k+1,l]-G_vis_hat[i,j,k+1,l]) - (G_adv_hat[i,j,k,l]-G_vis_hat[i,j,k,l]))
                end
            end
        end
    end
    return RHS
end