function set_wallpoint(nodes, bdcon, cellxmax, cellymax, cellzmax)
	# cal number of ponits
	nop = 0
	if Int(bdcon[1][1]) == 2 || Int(bdcon[1][1]) == 3
		nop = nop + (cellymax-1)*(cellzmax-1)
	end
	if Int(bdcon[2][1]) == 2 || Int(bdcon[2][1]) == 3
		nop = nop + (cellymax-1)*(cellzmax-1)
	end
	if Int(bdcon[3][1]) == 2 || Int(bdcon[3][1]) == 3
		nop = nop + (cellzmax-1)*(cellxmax-1)
	end
	if Int(bdcon[4][1]) == 2 || Int(bdcon[4][1]) == 3
		nop = nop + (cellzmax-1)*(cellxmax-1)
	end
	if Int(bdcon[5][1]) == 2 || Int(bdcon[5][1]) == 3
		nop = nop + (cellxmax-1)*(cellymax-1)
	end
	if Int(bdcon[6][1]) == 2 || Int(bdcon[6][1]) == 3
		nop = nop + (cellxmax-1)*(cellymax-1)
	end


	wallpoint = zeros(nop,3)
	ite = 1

	# x-
	if Int(bdcon[1][1]) == 2 || Int(bdcon[1][1]) == 3
		for j in 2:cellymax
			for k in 2:cellzmax
				wallpoint[ite,1] = nodes[2,j,k,1]
				wallpoint[ite,2] = nodes[2,j,k,2]
				wallpoint[ite,3] = nodes[2,j,k,3]
				ite = ite + 1
			end
		end
	end
	# x+
	if Int(bdcon[2][1]) == 2 || Int(bdcon[2][1]) == 3
		for j in 2:cellymax
			for k in 2:cellzmax
				wallpoint[ite,1] = nodes[cellzmax,j,k,1]
				wallpoint[ite,2] = nodes[cellzmax,j,k,2]
				wallpoint[ite,3] = nodes[cellzmax,j,k,3]
				ite = ite + 1
			end
		end
	end
	# y-
	if Int(bdcon[3][1]) == 2 || Int(bdcon[3][1]) == 3
		for k in 2:cellzmax
			for i in 2:cellxmax
				wallpoint[ite,1] = nodes[i,2,k,1]
				wallpoint[ite,2] = nodes[i,2,k,2]
				wallpoint[ite,3] = nodes[i,2,k,3]
				ite = ite + 1
			end
		end
	end
	# y+
	if Int(bdcon[4][1]) == 2 || Int(bdcon[4][1]) == 3
		for k in 2:cellzmax
			for i in 2:cellxmax
				wallpoint[ite,1] = nodes[i,cellymax,k,1]
				wallpoint[ite,2] = nodes[i,cellymax,k,2]
				wallpoint[ite,3] = nodes[i,cellymax,k,3]
				ite = ite + 1
			end
		end
	end
	# z-
	if Int(bdcon[5][1]) == 2 || Int(bdcon[5][1]) == 3
		for i in 2:cellxmax
			for j in 2:cellymax
				wallpoint[ite,1] = nodes[i,j,2,1]
				wallpoint[ite,2] = nodes[i,j,2,2]
				wallpoint[ite,3] = nodes[i,j,2,3]
				ite = ite + 1
			end
		end
	end
	# z+
	if Int(bdcon[6][1]) == 2 || Int(bdcon[6][1]) == 3
		for i in 2:cellxmax
			for j in 2:cellymax
				wallpoint[ite,1] = nodes[i,j,cellzmax,1]
				wallpoint[ite,2] = nodes[i,j,cellzmax,2]
				wallpoint[ite,3] = nodes[i,j,cellzmax,3]
				ite = ite + 1
			end
		end
	end
	return wallpoint, nop
end

function set_wally(cellcenter, wallpoint, nop, cellxmax, cellymax, cellzmax)
	wally = zeros(cellxmax, cellymax, cellzmax)

	distance_list = zeros(nop,2)
	distance = 0

	pickup_point = 5

	# 最大値を5つ抽出
	for i in 2:cellxmax-1
		for j in 2:cellymax-1
			for k in 2:cellzmax-1
				for n in 1:nop
						distance = ((cellcenter[i,j,k,1]-wallpoint[n,1])^2 +
									(cellcenter[i,j,k,2]-wallpoint[n,2])^2 +
									(cellcenter[i,j,k,3]-wallpoint[n,3])^2)^0.5
						distance_list[n,1] = n
						distance_list[n,2] = distance
				end

				qdistance_list = copy(distance_list)
				qdistance_list = quicksort(qdistance_list,1,nop,2,2)

				#println(distance_list[1,2])

				# 5wall_point
				wp = zeros(Int, pickup_point)
				for m in 1:pickup_point
					wp[m] = Int(qdistance_list[m,1])
				end
				
				# max(y)
				# 近傍5点を用いて一番遠いyが正解
				# 例1 1点はy+側、2点はy-側のとき、y-とy+で作れる三角形の高さが一番短いが間違い
				# 逆に一番遠いものが正解
				ytemp = zeros(10)
				ite1 = [1,1,1,1,1,1,2,2,2,3]
				ite2 = [2,2,2,3,3,4,3,3,4,4]
				ite3 = [3,4,5,4,5,5,4,5,5,5]

				# println(qdistance_list[5])
				for m in 1:10
					wp1 = wp[ite1[m]]
					wp2 = wp[ite2[m]]
					wp3 = wp[ite3[m]]
					
					dis1 = distance_list[wp1,2]
					dis2 = distance_list[wp2,2]
					dis3 = distance_list[wp3,2]
					
					wallpoint_dis1 = ((wallpoint[wp1,1]-wallpoint[wp2,1])^2 + (wallpoint[wp1,2]-wallpoint[wp2,2])^2 + (wallpoint[wp1,3]-wallpoint[wp2,3])^2)^0.5
					wallpoint_dis2 = ((wallpoint[wp2,1]-wallpoint[wp3,1])^2 + (wallpoint[wp2,2]-wallpoint[wp3,2])^2 + (wallpoint[wp2,3]-wallpoint[wp3,3])^2)^0.5
					wallpoint_dis3 = ((wallpoint[wp3,1]-wallpoint[wp1,1])^2 + (wallpoint[wp3,2]-wallpoint[wp1,2])^2 + (wallpoint[wp3,3]-wallpoint[wp1,3])^2)^0.5

					V = cal_tetrahedron_volume(dis1, dis2, dis3, wallpoint_dis1, wallpoint_dis2, wallpoint_dis3)

					S = cal_area(wallpoint_dis1, wallpoint_dis2, wallpoint_dis3)
					
					ytemp[m] = V/S
				end
				wally[i,j,k] = maximum(ytemp)
			end
		end
	end

	#=
	・壁面に下した垂線の計算について
	面積S = (s(s-a)(s-b)(s-c))^0.5 = hc/2  , h=垂線,c=底面(壁面の2点間距離)
	これで三辺の長さがわかれば垂線が算出できる
	=#
	return wally
end

function cal_tetrahedron_volume(a1,a2,a3,a4,a5,a6)
	temp1 = a1^2 * a5^2 * (a2^2 + a3^2 + a4^2 + a6^2 - a1^2 - a5^2)
	temp2 = a2^2 * a6^2 * (a1^2 + a3^2 + a4^2 + a5^2 - a2^2 - a6^2)
	temp3 = a3^2 * a4^2 * (a1^2 + a2^2 + a5^2 + a6^2 - a3^2 - a4^2)
	temp4 = -(a1*a2*a4)^2 - (a2*a3*a5)^2 - (a1*a3*a6)^2 - (a4*a5*a6)^2
	V2 = 1/144 * (temp1 + temp2 + temp3 + temp4)
	
	# ---------------------------------
	# 5つの近傍点で計算するため，構造格子だと，xが同じ値の場合がありエラーになるから
	# ---------------------------------
	volume = 0.0
	try volume = (V2)^0.5
	catch
		
	end
	
	return volume
end

function cal_area(a1,a2,a3)
	s = 0.5*(a1 + a2 + a3)

	# ---------------------------------
	# 5つの近傍点で計算するため，構造格子だと，xが同じ値の場合がありエラーになるから
	# ---------------------------------
	surface = s * (s-a1) * (s-a2) * (s-a3)
	
	if surface <= 0.0
		surface = 1.0e50
	end

	surface = surface^0.5

	return surface
end

function cal_yplus(Qbase, wally, mu, cellxmax, cellymax, cellzmax)
	yplus = zeros(cellxmax, cellymax, cellzmax)

	for i in 2:cellxmax-1
		for j in 2:cellymax-1
			for k in 2:cellzmax-1
				nu = mu[i,j,k] / Qbase[i,j,k,1]
				yplus[i,j,k], uplus = cal_wall_val_spalding(Qbase[i,j,k,2], wally[i,j,k], nu)
				
				if isequal(yplus[i,j,k], NaN) == true
					println("      ")
					println(Qbase[i,j,k,2])
					println(wally[i,j,k])
					println(nu)
					throw(UndefVarError(:x))
				end
			end
		end
	end

	return yplus
end

function cal_wall_val_spalding(u,y,nu)
	k = 0.4
	B = 5.5

	u_tau = 1.0
	for i in 1:20
		#=
		println("   zzz   ")
		println(u)
		println(y)
		println(nu)
		println(spalding(u_tau,u,y,nu,k,B))
		println(spalding_dash(u_tau,u,y,nu,k,B))
		println(u_tau)
		=#
			
		old = u_tau
		u_tau = u_tau-spalding(u_tau,u,y,nu,k,B)/spalding_dash(u_tau,u,y,nu,k,B) + 1.0e-50
		
		if abs(u_tau-old)/u_tau < 10^(-4)
			break
		end
	end

	yplus = u_tau*y/nu
	uplus = u/u_tau
	return yplus,uplus
end

function spalding(u_tau,u,y,nu,k,B)
	t = k*u/u_tau
	F = u/u_tau-u_tau*y/nu+exp(-k*B)*(exp(t)-1-t-t^2/2-t^3/6)
	return F
end

function spalding_dash(u_tau,u,y,nu,k,B)
	t = k*u/u_tau
	Fdash = -u/u_tau^2-y/nu+exp(-k*B)*(t/u_tau)*(-exp(t)+1+t+t^2/2)
	return Fdash
end


# test
#=
uu = 0.0
yy = 0.00014
nunu = 1.5e-5
println(cal_wall_val_spalding(uu,yy,nunu))
=#


function SGSstress()
	# 境界におけるtauが必要
end

function cal_Strain_rate_tensor()
	
end