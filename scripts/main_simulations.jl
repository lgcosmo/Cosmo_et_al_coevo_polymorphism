#This script runs the numerical simulations that return the results used to build Figure 3 of the manuscript.

using Distributed

addprocs(12) #Parallel computing

@everywhere using DrWatson
@everywhere @quickactivate "coevo_polymorphism"

@everywhere begin

    using DataFrames
    using Distributions
    using Statistics
    using LinearAlgebra
    using DelimitedFiles
    using CSV

    BLAS.set_num_threads(1)

    include(srcdir("main_functions.jl"))
    include(srcdir("aux_functions.jl"))

    mi_settings=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]

    networks_paths=readdir(srcdir("empirical_networks", "weboflife"); join=true)

    networks_names=readdir(srcdir("empirical_networks", "weboflife"))
    networks_names=replace.(networks_names, r".csv.*" => "")

    networks=readdlm.(networks_paths, ',')

    binary!(networks=networks)
    adj=adjacency(incidence_matrices=networks)
    q_list=Q_matrix.(adj)

    descriptors=networks_description(networks=networks, adj=adj, q_list=q_list, net_name=networks_names, matrix_type="incidence")

    c=Dict(:net_info => descriptors, :γ => mi_settings, :tmax => 5000, :ϵ => 1e-6, :nsim => 1000)
    c_list=dict_list(c)

end

pmap(net_multisim, c_list)