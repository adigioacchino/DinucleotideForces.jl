using DinucleotideForces
using Random
using Test

########################################################
# useful functions
########################################################
# approximate, fast way of computing CG forces (fitted params for human C,G freqs, 
# suppl mat of see https://academic.oup.com/mbe/article/38/6/2428/6130826)
function ApproxCGForceHfields(seq::AbstractString)
    L = length(seq)
    f_motif = count("CG", seq) / L
    return (log(f_motif / (0.201)^2) - 0.03) / 0.95
end

# remove motifs from the sequence (so that there are approx 1 - frac motifs)
function RemoveMotifs(seq::AbstractString, motif::AbstractString, 
                      frac::Float64=0.1, rng::Xoshiro=Xoshiro(1))
    seq_vec = string.(collect(seq))
    motif_poss = findall(motif, seq, overlap=true)
    nts = ["A","C","G","T"]
    for mp in motif_poss
        r = rand(rng)
        if r < frac
            pos = rand(rng, [1,2])
            good_nts = nts[nts .!= seq_vec[mp[pos]]]
            seq_vec[mp[pos]] = rand(rng, good_nts)
        end
    end
    return join(seq_vec)
end

# add motifs from the sequence (so that there are 1 + frac motifs)
function AddMotifs(seq::AbstractString, motif::AbstractString, 
                   frac::Float64=0.1, rng::Xoshiro=Xoshiro(1))
    seq_vec = string.(collect(seq))
    toadd = round(Int, count(motif, seq) * frac)
    added = 0
    while added < toadd
        i = rand(rng, 1:length(seq_vec))
        if join(seq_vec[i]*seq_vec[i+1]) != motif
            seq_vec[i] = string(motif[1])
            seq_vec[i+1] = string(motif[2])
        else
            continue
        end
        added += 1
    end
    return join(seq_vec)    
end
########################################################
########################################################


########################################################
# useful variables
########################################################
# this is an Influenza H5N1 PB2 segment, strain used: A/Anhui/1/2005
testseq = "ATGGAGAGAATAAAAGAATTAAGGGATCTAATGTCACAGTCCCGCACTCGCGAGATACTAACAAAAACCACTGTGGACCATATGGCCATAATCAAGAAGTACACATCAGGAAGACAAGAGAAGAACCCTGCTCTCAGAATGAAATGGATGATGGCAATGAAATATCCAATCACAGCGGACAAGAGAATAACAGAGATGATTCCTGAAAGGAATGAACAAGGGCAGACGCTCTGGAGCAAGACAAATGATGCCGGATCGGACAGGTTGATGGTGTCTCCCTTAGCTGTAACTTGGTGGAATAGGAATGGGCCGACGACAAGTGCAGTCCATTATCCAAAGGTTTACAAAACATACTTTGAGAAGGCTGAAAGGCTAAAACATGGAACCTTCGGTCCCGTCCATTTTCGAAACCAAGTTAAAATACGCCGCCGAGTTGATATAAATCCTGGCCATGCAGATCTCAGTGCTAAAGAAGCACAAGATGTCATCATGGAGGTCGTTTTCCCAAATGAAGTGGGAGCTAGAATATTGACATCAGAGTCACAATTGACAATAACGAAAGAGAAGAAAGAAGAGCTCCAAGATTGTAAGATTGCTCCCTTAATGGTTGCATACATGTTGGAAAGGGAACTGGTCCGCAAAACCAGATTCCTACCGGTAGCAAGCGGAACAAGCAGTGTGTACATTGAGGTATTGCATTTGACTCAAGGGACCTGCTGGGAACAGATGTACACTCCAGGCGGAGAAGTGAGAAACGACGATGTTGACCAGAGTTTGATCATCGCTGCCAGAAACATTGTTAGGAGAGCAACGGTATCAGCGGATCCACTGGCATCACTGCTGGAGATGTGTCACAGCACACAAATTGGTGGGATAAGGATGGTGGACATCCTTAGGCAAAACCCAACTGAGGAACAAGCTGTGGGTATATGCAAAGCAGCAATGGGTCTGAGGATCAGTTCATCCTTTAGCTTTGGAGGCTTCACTTTCAAAAGAACAAGTGGATCATCCGTCACGAAGGAAGAGGAAGTGCTTACAGGCAACCTCCAAACATTGAAAATAAGAGTACATGAGGGGTATGAAGAGTTCACAATGGTTGGACGGAGGGCAACAGCTATCCTGAGGAAAGCAACTAGAAGGCTGATTCAGTTGATAGTAAGTGGAAGAGACGAACAATCAATCGCTGAGGCAATCATTGTAGCAATGGTGTTCTCACAGGAGGATTGCATGATAAAGGCAGTCCGGGGCGATTTGAATTTCGTAAACAGAGCAAACCAAAGATTAAACCCCATGCATCAACTCCTGAGACATTTTCAAAAGGACGCAAAAGTGCTATTTCAGAATTGGGGAATTGAACCCATTGATAATGTCATGGGGATGATCGGAATATTACCTGACCTGACTCCCAGCACAGAAATGTCACTGAGAAGAGTAAGAGTTAGTAAAGTGGGAGTGGATGAATATTCCAGCACTGAGAGAGTAATTGTAAGTATTGACCGTTTCTTAAGGGTTCGAGATCAGCGGGGGAACGTACTCTTATCTCCCGAAGAGGTCAGCGAAACCCAGGGAACAGAGAAATTGACAATAACATATTCATCATCAATGATGTGGGAAATCAACGGTCCTGAGTCAGTGCTTGTTAACACCTATCAATGGATCATCAGAAACTGGGAAACTGTGAAGATTCAATGGTCTCAAGACCCCACGATGCTGTACAATAAGATGGAGTTTGAACCGTTCCAATCCTTGGTACCTAAGGCTGCCAGAGGTCAATACAGTGGATTTGTGAGAACACTATTCCAACAAATGCGTGACGTACTGGGGACATTTGATACTGTCCAGATAATAAAGCTGCTACCATTTGCAGCAGCCCCACCAGAGCAGAGCAGAATGCAGTTTTCTTCTCTAACTGTGAATGTGAGAGGCTCAGGAATGAGAATACTCGTAAGGGGCAATTCCCCTGTGTTCAACTACAATAAGGCAACCAAAAGGCTTACCGTTCTTGGAAAGGACGCAGGTGCATTAACAGAGGATCCAGATGAGGGGACAACCGGAGTGGAGTCTGCAGTACTGAGGGAATTCCTAATTCTAGGCAAGGAGGACAAAAGATATGGACCAGCATTGAGTATCAATGAACTGAGCAACCTTGCGAAAGGGGAGAAAGCTAATGTGCTGATAGGACAAGGAGACGTGGTGTTGGTAATGAAACGGAAACGGGACTCTAGCATACTTACTGACAGCCAGACAGCGACCAAAAGAATTCGGATGGCCATCAATTAG"
nts = ["A","C","G","T"]
all_dinucleotides = [x*y for x in nts for y in nts]
test_model = DimerForce(testseq, all_dinucleotides)
########################################################
########################################################


########################################################
# TESTS
########################################################
# test forces of random sequences
@testset "random uniform sequences" begin
    randseq = join(rand(Xoshiro(0), nts, 1_000_000))
    for dn in all_dinucleotides
        result = DimerForce(randseq, [dn], 1000)
        @test result[dn] ≈ 0. atol=0.05
        for n in nts
            @test exp(result[n]) ≈ 0.25 atol=0.05
        end
    end
    res_all = DimerForce(randseq, all_dinucleotides, 1000)
    for dn in all_dinucleotides
        @test res_all[dn] ≈ 0. atol=0.05
    end
end

# test sequences with 1 dinucleotide under or overexpressed
@testset "single-dinucleotide mutated sequences" begin
    randseq = join(rand(Xoshiro(0), nts, 1_000_000))
    # special treatment for CG for which I fitted ApproxCGForceHfields
    fs = [29.9, 20.1, 20.1, 29.9]/100
    for a in 0.5:0.1:0.9
        mut0 = RemoveMotifs(randseq, "CG", a)
        @test DimerForce(mut0, ["CG"], 1_000, freqs=fs)["CG"] ≈ ApproxCGForceHfields(mut0) atol=0.02
    end
    # for all the others, just check the sign
    for dn in all_dinucleotides
        mut1 = AddMotifs(randseq, dn)
        res1 = DimerForce(mut1, [dn], 1_000)
        @test res1[dn] > 0 
        mut2 = RemoveMotifs(randseq, dn)
        res2 = DimerForce(mut2, [dn], 1_000)
        @test res2[dn] < 0
    end
end

# take a sequence, infer a model, do variation of all parameters, check that loglik is lower
@testset "check on loglikelihoods" begin
    randseq = join(rand(Xoshiro(1), nts, 5_000))
    max_loglik_model = DimerForce(randseq, all_dinucleotides)
    loglik = DinucleotideForces.ComputeLoglikelihood(randseq, max_loglik_model)
    epsilon = 0.05
    for n in nts
        alt_model = copy(max_loglik_model)
        alt_model[n] = alt_model[n] + epsilon
        alt_loglik = DinucleotideForces.ComputeLoglikelihood(randseq, alt_model)
        @test loglik > alt_loglik
        alt_model = copy(max_loglik_model)
        alt_model[n] = alt_model[n] - epsilon
        alt_loglik = DinucleotideForces.ComputeLoglikelihood(randseq, alt_model)
        @test loglik > alt_loglik
    end
    for dn in all_dinucleotides
        alt_model = copy(max_loglik_model)
        alt_model[dn] = alt_model[dn] + epsilon
        alt_loglik = DinucleotideForces.ComputeLoglikelihood(randseq, alt_model)
        @test loglik > alt_loglik
        alt_model = copy(max_loglik_model)
        alt_model[dn] = alt_model[dn] - epsilon
        alt_loglik = DinucleotideForces.ComputeLoglikelihood(randseq, alt_model)
        @test loglik > alt_loglik
    end
end

# test multiseq: take a sequence, split and check that the result is consistent
@testset "model trained on multiple sequences" begin
    l = 760
    testseq_chunks = [testseq[i*l+1:(i+1)*l] for i in 0:2]
    chunks_model = DimerForce(testseq_chunks, all_dinucleotides)
    for k in keys(test_model)
        @test test_model[k] ≈ chunks_model[k] atol=0.02
    end
end

# test sampling & inference
@testset "testing sampling and inference" begin
    sampled_seqs = [DinucleotideForces.SampleSequence(test_model, 1_000, Xoshiro(i)) for i in 1:500]
    ss_model = DimerForce(sampled_seqs, all_dinucleotides)
    for k in keys(test_model)
        @test test_model[k] ≈ ss_model[k] atol=0.03
    end
end