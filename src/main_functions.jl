function net_coevolution(;net_info::DataFrame, γ, tmax::Int64, ϵ::Float64, sim::Int64)

    A=net_info.A[1]
    Q=net_info.Q[1]    
    net_name=net_info.network[1]
    n_a=net_info.n_a[1]
    n_p=net_info.n_p[1]
    n_sp=net_info.n_sp[1]
    degree=net_info.degree
    cp=net_info.cp
    teq=0

    dp=zeros(tmax, n_sp)
    pdif=zeros(n_sp)
    
    QA=zeros(n_sp, n_sp)
    QB=zeros(n_sp, n_sp)
    QA_sum=zeros(n_sp)
    QB_sum=zeros(n_sp)
    wA=zeros(n_sp)
    wB=zeros(n_sp)
    mean_w=zeros(n_sp)
    
    p0=rand(0.3:0.001:0.7, n_sp)
    θ=ones(n_sp, 2)

    if γ <= 0.5
        γB=rand(0.0:0.001:((2.0*γ)), n_sp)
        γA=(2.0*γ).-γB
    end

    if γ > 0.5
        γB=rand(((2.0*γ)-1.0):0.001:1.0, n_sp)
        γA=(2.0*γ).-γB
    end

    @views dp[1,:].=p0

   @inbounds for t in 1:(tmax-1)

        QAB!(dp=dp, Q=Q, γA=γA, γB=γB, QA=QA, QB=QB, time=t)
        coevo!(dp=dp, θ=θ, mean_w=mean_w, wA=wA, wB=wB, QA=QA, QB=QB, QA_sum=QA_sum, QB_sum=QB_sum, pdif=pdif, c=1.0, time=t)        

        if (mean(pdif) < ϵ) & (mean(pdif) > 0.0)
            teq=t+1
            break
        end

        #teq=tmax

    end

    @views peq=dp[teq,:]
    qeq=(1.0.-peq)
    pvar=peq.*qeq

    if n_p>0 && n_a>0
        type=vcat(repeat(["plant"], n_p), repeat(["animal"], n_a))
    elseif n_p==0 && n_a==0
        type=repeat(["unknown"], n_sp)
    end

    sp_id=["SP$i" for i in 1:n_sp]

    Tmat=T_matrix(Q=Q, θ=θ, γA=γA, γB=γB)

    T_diag=deepcopy(Tmat)
    T_diag[diagind(T_diag)].=0.0
    T_i=dropdims(sum(T_diag, dims=2), dims=2)
    F_i=T_i.-(γ./(1.0.+γ))

    df=DataFrame(network=net_name, sp_id=sp_id, type=type, position=cp, degree=degree, peq=peq, qeq=qeq, pvar=pvar, y=γ, yA=γA, yB=γB, T_i=T_i, F_i=F_i, sim=sim)

    return df

end

function net_multisim(p)
    
    @unpack γ,tmax,ϵ,nsim,net_info = p
    
    net_name=net_info.network[1]

    r=[DataFrame() for _ in 1:nsim]

    for n in 1:nsim
        r[n]=net_coevolution(net_info=net_info, γ=γ, tmax=tmax, ϵ=ϵ, sim=n)
    end

    CSV.write(datadir("sims", "$(net_name)_m$(γ).csv"), vcat(r...))

end