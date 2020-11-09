using Printf

function output_result(stepnum, Qbase, cellxmax, cellymax, cellzmax, out_file_front, out_ext, out_dir, Rd)
    
    stepnum = string(stepnum)
    while length(stepnum) < 6
        stepnum = "0"*stepnum
    end
    
    fff = out_dir*"/" * out_file_front*stepnum*out_ext
    open(fff,"w") do f
        write(f,"result:rho[kg/m^3], u[m/s], v[m/s], w[m/s], p[Pa], T[K]\n")
        for i in 2:cellxmax-1
            for j in 2:cellymax-1
                for k in 2:cellzmax-1
                    a1 = @sprintf("%8.8e", Qbase[i,j,k,1])
                    a2 = @sprintf("%8.8e", Qbase[i,j,k,2])
                    a3 = @sprintf("%8.8e", Qbase[i,j,k,3])
                    a4 = @sprintf("%8.8e", Qbase[i,j,k,4])
                    a5 = @sprintf("%8.8e", Qbase[i,j,k,5])
                    T = Qbase[i,j,k,5]/(Qbase[i,j,k,1]*Rd)
                    a6 = @sprintf("%8.8e", T)
                    write(f, a1*" "*a2*" "*a3*" "*a4*" "*a5*" "*a6*"\n")
                end
            end
        end
    end
    println("\nwrite "*fff)
end

function reset_write(fwrite)
    open( fwrite, "w" ) do f 
    end
end

function output_physicaltime(fwrite, t, dt)
    a = string(t)
    b = @sprintf("%8.8f", t*dt)
    open(fwrite, "a") do f
        write(f, "\n")
        write(f, "step = "*a)
        write(f, "     ")
        write(f, "time[s] = "*b)
        write(f, "\n")
        write(f, "\n")
    end
end

function output_innertime(fwrite, tau, norm2)
    a = string(tau)
    open(fwrite, "a") do f
        write(f, a)
        write(f, " ")
        for l in 1:5
            a = @sprintf("%8.8e", norm2[l])
            write(f, a)
            write(f, " ")
        end
        write(f, "\n")
    end
end