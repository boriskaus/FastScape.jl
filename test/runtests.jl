using Test, FastScape, Statistics, DelimitedFiles

@testset verbose = true "FastScape.jl" begin

    @testset "Mountain" begin
        include("Mountain.jl")
        @test maximum(h) ≈ 2528.0934077776587 rtol = 1e-6
        @test mean(h) ≈ 970.238605311561 rtol = 1e-6
    end

    @testset "Fan" begin
        include("Fan.jl")
        @test maximum(h) ≈ 664.0063768713197 rtol = 5e-2
        @test mean(h) ≈ 113.71736488098453 rtol = 5e-2
    end

    @testset "Margin" begin
        include("Margin.jl")
        @test mean(h) ≈ -498.10691586070374 rtol = 1e-6
    end

    @testset "DippingDyke" begin
        include("DippingDyke.jl")
        @test maximum(h) ≈ 2042.8661078613902 rtol = 1e-6
        @test mean(h) ≈ 697.1218846524354 rtol = 1e-6
    end
end