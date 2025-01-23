function binary!(;networks)
    for n in 1:length(networks)

       @views net=networks[n]

       for j in 1:size(net,2)
            for i in 1:size(net,1)
                if net[i,j]>0.0
                    net[i,j]=1.0
                end
            end
        end
    end
end

function adjacency(;incidence_matrices::AbstractArray)

    adj=deepcopy(incidence_matrices)
    
    for n in 1:length(incidence_matrices)

        incidence_matrix=incidence_matrices[n]

        n_p=size(incidence_matrix, 1)
        n_a=size(incidence_matrix, 2)

        A=vcat(hcat(zeros(n_p, n_p), incidence_matrix), hcat(incidence_matrix', zeros(n_a, n_a)))

        adj[n]=A

    end

    return adj

end

function Q_matrix(A::Array{Float64})

    Q=zeros(size(A,1), size(A,2))
    A_sum=dropdims(sum(A, dims=2), dims=2)

    @inbounds for j in 1:size(Q,2)
        @inbounds for i in 1:size(Q,1)

            Q[i,j]=A[i,j]/A_sum[i]

        end
    end

    return(Q)

end

function networks_description(;networks, adj, q_list, net_name, matrix_type)

    results=[DataFrame() for _ in 1:length(networks)]

    for i in 1:length(networks)

        network=networks[i]
        A=adj[i]
        Q=q_list[i]

        if matrix_type=="incidence"

            n_p=size(network, 1)
            n_a=size(network, 2)
            n_sp=n_p+n_a        
                    
        elseif matrix_type=="adjacency"
        
            n_sp=size(network,1)
            n_p=0
            n_a=0            

        end      
                
        if matrix_type=="incidence"

            p_degree=dropdims(sum(network, dims=2), dims=2)
            a_degree=dropdims(sum(network, dims=1), dims=1)

            degree=vcat(p_degree, a_degree)

            cp_plant=ifelse.(p_degree.>=mean(p_degree), "core", "periphery")
            cp_animal=ifelse.(a_degree.>=mean(a_degree), "core", "periphery")

            cp=vcat(cp_plant, cp_animal)
        
        elseif matrix_type=="adjacency"
                
            cp=["unknown" for _ in 1:n_sp]
            degree=dropdims(sum(network, dims=2), dims=2)

        end

        results[i]=DataFrame(A=repeat([A], n_sp), Q=repeat([Q], n_sp), network=net_name[i], n_p=n_p, n_a=n_a, n_sp=n_sp, degree=degree, cp=cp)

        
    end

    return results
        
end

function QAB!(;dp::Array{Float64}, Q::Array{Float64}, γA::Array{Float64}, γB::Array{Float64}, QA::Array{Float64}, QB::Array{Float64}, time::Int64)

    @inbounds for j in 1:size(Q,2)
        @inbounds for i in 1:size(Q,1)
            QA[i,j]=(γA[i]*Q[i,j]*(dp[time,j]))
            QB[i,j]=(γB[i]*Q[i,j]*(1.0-dp[time,j]))
        end
    end

end

function coevo!(;dp::Array{Float64}, θ::Array{Float64}, mean_w::Array{Float64}, wA::Array{Float64}, wB::Array{Float64}, QA::Array{Float64}, QB::Array{Float64}, QA_sum::Array{Float64}, QB_sum::Array{Float64}, pdif::Array{Float64}, c::Float64, time::Int64)

    sum!(QA_sum, QA)
    sum!(QB_sum, QB)

    @inbounds for i in 1:size(θ,1)

        wA[i]=c+QA_sum[i]-θ[i,1]*dp[time,i]
        wB[i]=c+QB_sum[i]-θ[i,2]*(1.0-dp[time,i])
        mean_w[i]=dp[time,i]*wA[i] + (1.0-dp[time,i])*wB[i]

        dp[time+1,i]=dp[time,i]+((dp[time,i]*(1.0-dp[time,i])*(wA[i]-wB[i]))/mean_w[i])
        pdif[i]=abs(dp[time+1, i] - dp[time, i])

    end

end

function T_matrix(;Q::Array{Float64}, θ::Array{Float64}, γA::Array{Float64}, γB::Array{Float64})
    
    QT=zeros(size(Q,1), size(Q,2))

    for j in 1:size(Q,2)
        for i in 1:size(Q,1)
            QT[i,j]=(Q[i,j]*(γA[i]+γB[i]))/(θ[i,1]+θ[i,2])
        end
    end
    
    Θ=Diagonal(1.0.-(dropdims(sum(QT, dims=2), dims=2)))
    T=inv(I-QT)*Θ

    return T

end