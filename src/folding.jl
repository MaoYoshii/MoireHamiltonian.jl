using UnPack
using LinearAlgebra


export  Moire,
        lattice_img

"""
    Moire(Gmax,nmax,lattice₁,lattice₂)
    Gmax    : the maximum norm of folding vector
    nmax    : the range to search folding vectors
    lattice : constructors of lattice
"""
struct Moire
    fold :: AbstractArray{Int}              
    G₁ :: AbstractArray{Real}               # reciprocal lattice of layer 1
    G₂ :: AbstractArray{Real}               # reciprocal lattice of layer 1
    numfold :: Int                          # fold k space of layer 2 for numfold times
    M :: Int                          # fold k space of layer 2 for numfold times
    #matrix :: AbstractArray{ComplexF64}     # matrix for moire hamiltonian

    function Moire(Gmax,nmax,lattice₁,lattice₂)
        function make_fold(Gmax,nmax,b)
            G = Array{Float64,2}(undef,(2nmax+1)^2,2)
            F = Array{Int,2}(undef,(2nmax+1)^2,2)
            index = 0
            @inbounds for g₁ in -nmax:nmax
                @inbounds for g₂ in -nmax:nmax
                    g = [g₁ g₂]*b
                    if norm(g) ≤ Gmax
                        index += 1
                        F[index,:] = [g₁ g₂]
                        G[index,:] = g
                    end
                end
            end
            return G[1:index,:] , F[1:index,:]
        end
        G₁ , fold=make_fold(Gmax,nmax,lattice₁.b)
        new(fold,
            G₁,
            similar(G₁),
            size(G₁,1),
            lattice₁.M+lattice₂.M*size(G₁,1)
        )
    end
end

"""
    lattice_img(Moire,lattice₁)
    lattice₂    : lattice constructor of layer1
"""
function lattice_img(Moire,lattice₁)
    @unpack a, τ, b, Ham!, M= lattice₁
    @unpack fold, G₁, G₂, numfold , M = Moire
    
    p1 = plot(title = "Lattice") ; p2 = plot(title="Reciprocal lattice");
    lat = fold*a' ; G = fold*b
    scatter!(p1,lat[:,1] , lat[:,2]) ; scatter!(p2,G[:,1],G[:,2])
    plot(p1,p2,aspect_ratio=:equal,msw = 0,ms=2,mc = :blue)
end