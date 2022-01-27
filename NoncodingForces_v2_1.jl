# added the possibility of inferring fields
module NoncodingForces_v2_1

using Octavian
using LinearAlgebra
using FiniteDiff
using Statistics

#export 

const dna_alphabet = ['A', 'C', 'G', 'T']
const rna_alphabet = ['A', 'C', 'G', 'U']
const len_alphabet = 4


"""
Return the matrices used for the computation of the
partition function for each motif in motifs. In particular, M_{ij} is the number of times 
the given motif appears in n(i) + n(j), with n(i) the i-th nucleotide.
"""
function generate_motsM(motifs::Vector{String}, using_rna::Bool=false)
    if using_rna
        [[count(m, n1*n2, overlap=true) for n1 in rna_alphabet, n2 in rna_alphabet] for m in motifs]
    else
        [[count(m, n1*n2, overlap=true) for n1 in dna_alphabet, n2 in dna_alphabet] for m in motifs]
    end
end


"""
Return a matrix with the informations about the nt frequences, so that
when it is multiplied elementwise with that obtained through 
generate_motsM (after implementing the exponentiation with the force, see eval_log_Z),
the transfer matrices used to compute Z are finally obtained.
"""
function generate_nucsM(fields::Vector{Float64}, last::Bool)
    if !last
        [exp(f) for f in fields, i in 1:len_alphabet]
    else
        [exp(f1 + f2) for f1 in fields, f2 in fields]
    end
end


"""Compute the transfer matrix through generate_motsM, generate_freqsM and the forces given,
then takes the correct power to obtian Z."""
function eval_log_Z(fields, forces, motifs, L)
    v = ones(len_alphabet)
    motsM = generate_motsM(motifs)
    M = ones((len_alphabet, len_alphabet))
    for j in 1:length(forces)
        M .*= exp.(forces[j] * motsM[j])
    end
    TM = M .* generate_nucsM(fields, false) # this is the transfer matrix
    last_mat = M .* generate_nucsM(fields, true) # last matrix is special
    v = last_mat * v
    tTM = Matrix{Float64}(I, len_alphabet, len_alphabet)
    log_factors = 0
    for i in 1: L-2
        if i%10 == 0 # each 10 steps normalize v and save log of norm in log_factors, to avoid overflow for long sequences
            f = norm(tTM)
            log_factors += log(f)
            matmul!(tTM, tTM, TM, 1/f)
        else
            matmul!(tTM, tTM, TM)
        end
    end
    v = tTM * v
    return log(sum(v)) + log_factors
end


"""
This function returns the nucleotides or motifs that have to be inferred after
as many variables as possible are fixed to 0 through gauge transformations.
"""
function GaugeAwayVariables(motifs::Vector{String}, fields::Bool; using_rna=false)
    curr_alphabet = using_rna ? rna_alphabet : dna_alphabet
    n_motifs = length(motifs)
    if fields
        # check which nts I must include (no gauge transf to set them to 0)
        nt_gauge_mask = []
        gauge_dof_available = true
        for (i, nt) in enumerate(curr_alphabet)
            if sum(getindex.(motifs, 1) .== nt) == 4
                push!(nt_gauge_mask, 0)        
            elseif gauge_dof_available # use K_1 to fix to 0 when cannot use K_2
                push!(nt_gauge_mask, 0)
                gauge_dof_available = false
            else
                push!(nt_gauge_mask, 1)
            end
        end
        nt_gauge_mask = BitArray(nt_gauge_mask)
    end
    # check which motifs I must include (no gauge transf to set them to 0)
    mot_gauge_mask = ones(n_motifs)
    dof_used = 0
    for (i, nt) in enumerate(curr_alphabet)
        mask1 = (getindex.(motifs, 1) .== nt) .& (getindex.(motifs, 1) .!= getindex.(motifs, 2))
        mask2 = (getindex.(motifs, 2) .== nt) .& (getindex.(motifs, 1) .!= getindex.(motifs, 2))
        if (sum(mask1) + sum(mask2) == 6) & (dof_used < 3) # gauge can be fixed
            target = (getindex.(motifs, 1) .== nt) .& (getindex.(motifs, 2) .== 'T')
            target_pos = collect(1:n_motifs)[target][1]
            mot_gauge_mask[target_pos] = 0
            dof_used += 1
        end  
    end
    # if also fields are inferred, I might have one remaining dof
    if fields && gauge_dof_available
        target = (getindex.(motifs, 1) .== 'T') .& (getindex.(motifs, 2) .== 'T')
        target_pos = collect(1:n_motifs)[target][1]
        mot_gauge_mask[target_pos] = 0
    elseif !fields & (n_motifs == length(curr_alphabet)^2)# the fact that the sum of motifs is L-1 is not used if fields are not inferred
        target = (getindex.(motifs, 1) .== 'T') .& (getindex.(motifs, 2) .== 'T')
        target_pos = collect(1:n_motifs)[target][1]
        mot_gauge_mask[target_pos] = 0
    end
    mot_gauge_mask = BitArray(mot_gauge_mask)
    #fields && println(nt_gauge_mask)
    #println(mot_gauge_mask)
    if fields
        return curr_alphabet[nt_gauge_mask], motifs[mot_gauge_mask]
    else
        return motifs[mot_gauge_mask]
    end
end


"""
Returns a new dictionary with all the fields/couplings inferred,
such that the exp of the fields sum to 1 and as much information
as possible is included in the field, by putting other motifs to 0
(if possible, those of the form AA, CC, GG, TT and those ending 
with T).
"""
function FixFinalGauge!(dict_vars::Dict{String, Float64}; using_rna=false)
    curr_alphabet = using_rna ? string.(rna_alphabet) : string.(dna_alphabet) 
    all_vars = collect(keys(dict_vars))
    nts = [x for x in all_vars if length(x)==1]
    mots = [x for x in all_vars if length(x)==2]
    for (i, nt) in enumerate(curr_alphabet)
        if !(nt in nts)
            dict_vars[nt] = 0
        end
    end    
    # fix gauge (lattice-gas T gauge)
    # step1
    for (i, nt) in enumerate(curr_alphabet)
        mask = (string.(getindex.(mots, 1)) .== nt)
        if sum(mask) == 4 # gauge can be fixed
            to_be_modified = mots[mask]
            k_2 = - dict_vars[nt * nt]
            dict_vars[nt] -= k_2
            [dict_vars[mot] += k_2 for mot in to_be_modified]
        end
    end
    # step2
    k_1 = - log(sum(exp.([dict_vars[n] for n in curr_alphabet])))
    [dict_vars[n] += k_1 for n in curr_alphabet]
    # step3
    for (i, nt) in enumerate(curr_alphabet[1:3])
        mask1 = (string.(getindex.(mots, 1)) .== nt) .& (getindex.(mots, 1) .!= getindex.(mots, 2))
        mask2 = (string.(getindex.(mots, 2)) .== nt) .& (getindex.(mots, 1) .!= getindex.(mots, 2))
        if sum(mask1) + sum(mask2) == 6 # gauge can be fixed
            to_be_modified1 = mots[mask1]
            to_be_modified2 = mots[mask2]
            k_3 = - dict_vars[nt * "T"]
            [dict_vars[mot] += k_3 for mot in to_be_modified1]
            [dict_vars[mot] -= k_3 for mot in to_be_modified2]
        end 
    end
    return dict_vars
end


function DimerForce_withFields(seq::Union{String,Vector{String}}, nucleotides::Vector{Char}, motifs::Vector{String}; 
                     tolerance::Float64=0.01, max_iter::Int=100, add_pseudocount=true, using_rna=false)
    curr_alphabet = using_rna ? rna_alphabet : dna_alphabet
    n_nucleotides = length(nucleotides)
    n_motifs = length(motifs)
    if typeof(seq) == String
        seqs = [seq]
    else
        # all sequences used should have the same length
        seqs = seq
        @assert all(length.(seqs) .== length(seqs[1])) "Sequences provided must all have the same length."
    end
    L = length(seqs[1]) 
    n_obs_nucs = [mean(count.(string(m), seqs)) for m in nucleotides]
    n_obs_mots = [mean(count.(m, seqs, overlap=true)) for m in motifs]   
    if add_pseudocount
        n_obs_nucs = [x+1 for x in n_obs_nucs]
        n_obs_mots = [x+1 for x in n_obs_mots]
    end    
    n_obs = [n_obs_nucs; n_obs_mots]
    # I need a function of a unique vector, so I use a "closure", a function defined within function to have access to freqs, motifs, L
    function closed_eval_log_Z(x)
        fields = zeros(4)
        k = 1
        for (i, nt) in enumerate(curr_alphabet)
            if nt in nucleotides
                fields[i] = x[k]
                k+=1
            end
        end
        forces = x[k:end]
        return eval_log_Z(fields, forces, motifs, L)
    end
    ###########################
    #vars = vcat(fill(log(1/n_nucleotides), n_nucleotides), zeros(n_motifs))
    vars =  zeros(n_nucleotides + n_motifs)
    ###########################
    for l in 1:max_iter
        ns = FiniteDiff.finite_difference_gradient(closed_eval_log_Z, vars)        
        dn = FiniteDiff.finite_difference_hessian(closed_eval_log_Z, vars)
        delta = inv(dn) * (n_obs .- ns)
        if maximum(abs.(n_obs .- ns)) <= tolerance
            #println("Done! Exiting...")
            break
        end        
        vars .+= delta
        ############################
        ## restore 1-sum gauge for fields
        #k_1 = - log(sum(exp.(vars[1:4])))
        #vars[1:4] .+= k_1
        ############################
    end

    # format result
    res = Dict(zip(string.(nucleotides), vars[1:n_nucleotides]))
    for (i, m) in enumerate(motifs)
        res[m] = vars[n_nucleotides+i]
    end
    return res
end


function DimerForce_onlyForces(seq::Union{String,Vector{String}}, motifs::Vector{String}, fields::Vector{Float64}; 
                               tolerance::Float64=0.01, max_iter::Int=100, add_pseudocount=true, using_rna=false)
    n_motifs = length(motifs)
    if typeof(seq) == String
        seqs = [seq]
    else
        # all sequences used should have the same length
        seqs = seq
        @assert all(length.(seqs) .== length(seqs[1])) "Sequences provided must all have the same length."
    end
    L = length(seqs[1]) 
    n_obs = [mean(count.(m, seqs, overlap=true)) for m in motifs]           
    if add_pseudocount
        n_obs = [x+1 for x in n_obs]
    end
    # I need a function of a unique vector, so I use a "closure", a function defined within function to have access to freqs, motifs, L
    function closed_eval_log_Z(fs)
        return eval_log_Z(fields, fs, motifs, L)
    end
    forces = zeros(n_motifs)       
    for l in 1:max_iter
        ns = FiniteDiff.finite_difference_gradient(closed_eval_log_Z, forces)        
        dn = FiniteDiff.finite_difference_hessian(closed_eval_log_Z, forces)
        df = inv(dn) * (n_obs .- ns)
        if maximum(abs.(n_obs .- ns)) <= tolerance
            break
        end        
        forces .+= df
    end
    
    # format result
    return Dict(zip(motifs, forces))
end


"""
    DimerForce(seq::Union{String,Vector{String}}, motifs::Vector{String}; freqs=missing, tolerance::Float64=0.01, 
                     max_iter::Int=100, add_pseudocount=false, using_rna=false)
If frequencies are given, return the forces on the motifs 'motifs' computed for sequence 'seq' with 
the frequency bias given.
If frequencies are not given, the local fields (NOT frequencies!) are inferred and returned (in a
gauge such that the exp of them sum to 1, so the fields can be interpreted as log of frequencies).
Notice that seq is used only to compute the number of observed motifs (and nucleotides if the frequencies
are not provided by user). 
Also notice that if motifs are dinucleotides, only up to 12 of them are 
independent. 
Finally, tolerance and max_iter are parameters for the Newton-Raphson algorithm 
used to solve the system of equations.
If add_pseudocount, a single pseudocount is added for each observed number of
nucleotides and dinucleotides.
seq can be a single sequence or a vector of sequences (in this second case, the sequences
must have all the same length).
"""
function DimerForce(seq::Union{String,Vector{String}}, motifs::Vector{String}; freqs=missing, tolerance::Float64=0.01, 
                     max_iter::Int=100, add_pseudocount=false, using_rna=false)
    if ismissing(freqs) # infer also fields
        new_nts, new_mots = GaugeAwayVariables(motifs, true; using_rna)
        discarded_mots = [mot for mot in motifs if !(mot in new_mots)]
        t_res = DimerForce_withFields(seq, new_nts, new_mots; tolerance=tolerance, max_iter=max_iter, 
                                     add_pseudocount=add_pseudocount, using_rna=using_rna)
        [t_res[mot] = 0 for mot in discarded_mots] # include in final output motifs put to 0 through gauge
        return FixFinalGauge!(t_res; using_rna)
    else # do not infer fields
        @assert isapprox(sum(freqs), 1) "When frequencies are provided, they must sum to 1."
        new_mots = GaugeAwayVariables(motifs, false; using_rna)
        discarded_mots = [mot for mot in motifs if !(mot in new_mots)]
        # from freqs to fields
        fields = log.(freqs)
        t_res = DimerForce_onlyForces(seq, new_mots, fields; tolerance=tolerance, max_iter=max_iter, 
                                     add_pseudocount=add_pseudocount, using_rna=using_rna)
        [t_res[mot] = 0 for mot in discarded_mots] # include in final output motifs put to 0 through gauge
        return t_res
    end
end


"""
Given a sequence seq, the single-nucleotide fields and the dinucleotide forces (as a single dict), 
compute the energy of this sequence.
"""
function compute_energy(seq::String, fields_forces::Dict{String, Float64})
    e = 0
    for k in keys(fields_forces)
        e += count(k, seq) * fields_forces[k]
    end
    return e
end


"""
Given a sequence seq, the single-nucleotide fields and the dinucleotide forces (as a single dict),
compute the log-likelihood (energy minus log of Z) of this sequence.
logZ can be passed directly if pre-computed, otherwise is it computed each time this function is called.
"""
function compute_loglikelihood(seq::String, fields_forces::Dict{String, Float64}; logZ=missing, using_rna=false)
    e = compute_energy(seq, fields_forces)
    if ismissing(logZ)
        ks = keys(fields_forces)
        if !using_rna
            k1 = string.(dna_alphabet)
        else
            k1 = string.(rna_alphabet)
        end
        k2 = [k for k in ks if length(k)==2]
        fields = [fields_forces[k] for k in k1]
        forces = [fields_forces[k] for k in k2]
        L = length(seq)
        logz = eval_log_Z(fields, forces, k2, L)
    else
        logz = logZ
    end
    return e - logz
end


"""
    marginal_2points(fields_forces::Dict{String, Float64}, L::Int, pos::Int; using_rna=false)
Given a model specified by fields and forces and sequence length, compute the 2-point marginal in position
(i-1, i). The result is returned as a dictionary "s_(i-1) s_i" => probability.
"""
function marginal_2points(fields_forces::Dict{String, Float64}, L::Int, pos::Int; using_rna=false)
    # collect fields and forces
    ks = keys(fields_forces)
    if !using_rna
        k1 = string.(dna_alphabet)
    else
        k1 = string.(rna_alphabet)
    end
    k2 = [k for k in ks if length(k)==2]
    fields = [fields_forces[k] for k in k1]
    forces = [fields_forces[k] for k in k2]
    
    # build transfer matrices
    motsM = generate_motsM(k2)
    M = ones((len_alphabet, len_alphabet))
    for j in 1:length(forces)
        M .*= exp.(forces[j] * motsM[j])
    end
    TM = M .* generate_nucsM(fields, false) # this is the transfer matrix
    last_mat = M .* generate_nucsM(fields, true) # last matrix is special

    # compute marginals
    @assert (pos >= 2) & (pos <= L) "Please provide a position 'pos' between 2 and L (included)."
    pre_res = 0
    if pos == 2
        v = last_mat * ones(len_alphabet)
        tTM = Matrix{Float64}(I, len_alphabet, len_alphabet)
        log_factors = 0
        for i in 1:(L-3)
            if i%10 == 0 # each 10 steps normalize v and save log of norm in log_factors, to avoid overflow for long sequences
                f = norm(tTM)
                log_factors += log(f)
                matmul!(tTM, tTM, TM, 1/f)
            else
                matmul!(tTM, tTM, TM)
            end
        end
        v = tTM * v
        pre_res = transpose(transpose(TM) .* v)
        # take the log, reinsert the log_factors
        pre_res = log.(pre_res) .+ log_factors
        # normalize
        pre_res = pre_res .- sign(pre_res[1]) * maximum(abs.(pre_res)) # to have logprobs closer to 0 
        pre_res = exp.(pre_res) ./ sum(exp.(pre_res))
    elseif pos == L
        v = transpose(ones(len_alphabet)) * TM         
        tTM = Matrix{Float64}(I, len_alphabet, len_alphabet)
        log_factors = 0
        for i in 1:(L-3)
            if i%10 == 0 # each 10 steps normalize v and save log of norm in log_factors, to avoid overflow for long sequences
                f = norm(tTM)
                log_factors += log(f)
                matmul!(tTM, tTM, TM, 1/f)
            else
                matmul!(tTM, tTM, TM)
            end
        end
        v = v * tTM
        pre_res = transpose(transpose(last_mat) .* v)
        # take the log, reinsert the log_factors
        pre_res = log.(pre_res) .+ log_factors
        # normalize
        pre_res = pre_res .- sign(pre_res[1]) * maximum(abs.(pre_res)) # to have logprobs closer to 0
        pre_res = exp.(pre_res) ./ sum(exp.(pre_res))        
    else
        v2 = last_mat * ones(len_alphabet)
        log_factors2 = 0
        tTM = Matrix{Float64}(I, len_alphabet, len_alphabet)
        for i in 1:(L-1-pos)
            if i%10 == 0 # each 10 steps normalize v and save log of norm in log_factors, to avoid overflow for long sequences
                f = norm(tTM)
                log_factors2 += log(f)
                matmul!(tTM, tTM, TM, 1/f)
            else
                matmul!(tTM, tTM, TM)
            end
        end
        v2 = tTM * v2
        v1 = transpose(ones(len_alphabet)) * TM 
        log_factors1 = 0
        tTM = Matrix{Float64}(I, len_alphabet, len_alphabet)
        for i in 1:(pos-3)
            if i%10 == 0 # each 10 steps normalize v and save log of norm in log_factors, to avoid overflow for long sequences
                f = norm(tTM)
                log_factors1 += log(f)
                matmul!(tTM, tTM, TM, 1/f)
            else
                matmul!(tTM, tTM, TM)
            end
        end
        v1 = v1 * tTM
        log_factors = log_factors1 + log_factors2
        pre_res = transpose(transpose(TM) .* v1)
        pre_res = transpose(transpose(pre_res) .* v2)
        # take the log, reinsert the log_factors
        pre_res = log.(pre_res) .+ log_factors
        # normalize
        pre_res = pre_res .- sign(pre_res[1]) * maximum(abs.(pre_res)) # to have logprobs closer to 0
        pre_res = exp.(pre_res) ./ sum(exp.(pre_res))        
    end
    pre_res
    # prepare output
    res = Dict{String, Float64}()
    for (i,a) in enumerate(k1)
        for (j,b) in enumerate(k1)
            res[a * b] = pre_res[i,j]
        end
    end
    return res
end


"""
    marginal_1point(marginal_2points; first::Bool=true)
Given the 2point marginal, sum over the second symbol (if first=true, otherwise 
over the first) to obtain the 1-point marginal.
"""
function marginal_1point(marginal_2points::Dict{String, Float64}; first::Bool=true)
    alphabet = string.(unique(collect(prod(keys(marginal_2points)))))
    mat = Array{Float64}(undef, 4,4)
    for (i,a) in enumerate(alphabet)
        for (j,b) in enumerate(alphabet)
            mat[i,j] = marginal_2points[a*b]
        end
    end
    if first
        return Dict(zip(alphabet, vec(sum(mat, dims=2))))
    else
        return Dict(zip(alphabet, vec(sum(mat, dims=1))))
    end
end


"""
    marginal_1given_previous(marginal_2points::Dict{String, Float64}; first::Bool=true)
Given the 2point marginal, fix the first symbol (if `first`=true, otherwise 
fix the second) to obtain the probability of the second given that the first is equal
to `fixed` (or viceversa if `first`=false).
"""
function marginal_1given_previous(marginal_2points::Dict{String, Float64}, fixed::String; first::Bool=true)
    alphabet = string.(unique(collect(prod(keys(marginal_2points)))))
    mat = Array{Float64}(undef, 4,4)
    for (i,a) in enumerate(alphabet)
        for (j,b) in enumerate(alphabet)
            mat[i,j] = marginal_2points[a*b]
        end
    end
    @assert fixed in alphabet "Please fix the symbol as one from: $(alphabet)."
    if first
        row = findfirst(fixed .== alphabet)
        t_res = mat[row,:]
        res = t_res ./ sum(t_res)
    else
        col = findfirst(fixed .== alphabet)
        t_res = mat[:,col]
        res = t_res ./ sum(t_res)
    end
    return Dict(zip(alphabet, res))
end




end # end of the package
