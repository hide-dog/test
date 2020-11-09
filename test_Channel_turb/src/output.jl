using Printf

function output_result(stepnum, Qbase, cellxmax, cellymax, cellzmax, out_file_front, out_ext, out_dir, Rd, yplus)
    
    stepnum = string(stepnum)
    while length(stepnum) < 6
        stepnum = "0"*stepnum
    end
    
    fff = out_dir*"/" * out_file_front*stepnum*out_ext
    open(fff,"w") do f
        write(f,"result:rho[kg/m^3], u[m/s], v[m/s], w[m/s], p[Pa], T[K], yplus[-]\n")
        for i in 2:cellxmax-1
            for j in 2:cellymax-1
                for k in 2:cellzmax-1
                    a1 = @sprintf("%8.8f", Qbase[i,j,k,1])
                    a2 = @sprintf("%8.8f", Qbase[i,j,k,2])
                    a3 = @sprintf("%8.8f", Qbase[i,j,k,3])
                    a4 = @sprintf("%8.8f", Qbase[i,j,k,4])
                    a5 = @sprintf("%8.8f", Qbase[i,j,k,5])
                    T = Qbase[i,j,k,5]/(Qbase[i,j,k,1]*Rd)
                    a6 = @sprintf("%8.8f", T)
                    a7 = @sprintf("%8.8f", yplus[i,j,k])
                    write(f, a1*" "*a2*" "*a3*" "*a4*" "*a5*" "*a6*" "*a7*"\n")
                end
            end
        end
    end
    println("\nwrite "*fff)
end