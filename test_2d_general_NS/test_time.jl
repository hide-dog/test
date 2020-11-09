using Dates
using Printf

t1 = now()

sleep(3)

t2 = now()

etime = t2 - t1
outtime = ["temp"]
outtime[1] = string(etime)

et = replace(outtime[1], " milliseconds"=>"")
et = parse(Float64, et)

# et_s = string( Int(floor(et / 10^3)))
# et_m = string( Int(floor(et / 10^3 / 60)))
# et_h = string( Int(floor(et / 10^3 / 60 /60)))

et_h = floor(et / 10^3 / 60 /60)
et = et - et_h * (10^3 * 60 * 60)
et_m = floor(et / 10^3 / 60 )
et = et - et_m * (10^3 * 60)
et_s = floor(et / 10^3)
et = et - et_s * (10^3)

et_h = @sprintf( "%02d", et_h)
et_m = @sprintf( "%02d", et_m)
et_s = @sprintf( "%02d", et_s)
et = et_h*":"*et_m*":"*et_s

println("\n")
println("\n")
println("---------------------------- \n")
println("  fin program  ")
println("\n")
println("  elapse time (hh:mm:ss) = " * et)
println("\n")
println("---------------------------- \n")
