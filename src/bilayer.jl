using UnPack

export  Bilayer,
        bilayer!

exppi(x) = cospi(x) + im*sinpi(x)

function phase(G,τ)
    if size(τ,2) == 1
        return exppi(G'*τ)
    else
        return exppi.(G'*τ)
    end
end

function bilayer!(kx,ky,V!,lattice₁,lattice₂,Moire,Matrix)
    @unpack fold, G₁, G₂, numfold, M = Moire
    @. Matrix = 0.
    G₂ = ( G₁/lattice₂.b .|> x -> round(x,RoundNearest) ) |> x -> x*lattice₂.b # 1/2は近い偶数になる
    @inbounds Threads.@threads for g in 1:numfold
        from = g*lattice₂.M+1 ; to = g*lattice₂.M+lattice₁.M
        V!(kx+G₁[g,1],ky+G₁[g,2],view(Matrix,1:lattice₁.M,from:to))
        rmul!( view(Matrix ,1:lattice₁.M,from:to) , phase( view(-G₂,g,:) , lattice₂.τ ) )
        lmul!( phase( view(G₁,g,:) , lattice₁.τ ) , view(Matrix ,1:lattice₁.M,from:to) )
    end
    @. Matrix += Matrix'
    from  = 1; to = lattice₁.M
    lattice₁.Ham!(kx,ky,view(Matrix,from:to,from:to))
    @inbounds Threads.@threads for g in 1:numfold
        from = g*lattice₂.M+1 ; to = g*lattice₂.M+lattice₁.M
        lattice₂.Ham!(kx+G₁[g,1],ky+G₁[g,2],view(Matrix,from:to,from:to))
    end
end

struct Bilayer
    Ham!
    function Bilayer(V,lattice₁,lattice₂,Moire,inplace)
        @unpack fold, G₁, G₂, numfold, M = Moire
        if inplace
            function _bilayer!(kx,ky,Matrix)
                @. Matrix = 0.
                G₂ = ( G₁/lattice₂.b .|> x -> round(x,RoundNearest) ) |> x -> x*lattice₂.b # 1/2は近い偶数になる
                @inbounds for g in 1:numfold
                    from = lattice₁.M+(g-1)*lattice₂.M+1 ; to = lattice₁.M+g*lattice₂.M
                    V(kx+G₁[g,1],ky+G₁[g,2],view(Matrix,1:lattice₁.M,from:to))
                    rmul!( view(Matrix ,1:lattice₁.M,from:to) , phase( view(-G₂,g,:) , lattice₂.τ ) )
                    lmul!( phase( view(G₁,g,:) , lattice₁.τ ) , view(Matrix ,1:lattice₁.M,from:to) )
                end
                @. Matrix += Matrix'
                from  = 1; to = lattice₁.M
                lattice₁.Ham!(kx,ky,view(Matrix,from:to,from:to))
                @inbounds for g in 1:numfold
                    from = lattice₁.M+(g-1)*lattice₂.M+1 ; to = lattice₁.M+g*lattice₂.M
                    lattice₂.Ham!(kx+G₁[g,1],ky+G₁[g,2],view(Matrix,from:to,from:to))
                end
            end
            return new(_bilayer!)
        
        else
            function _bilayer(kx,ky)
                Matrix = zeros(ComplexF64,lattice₁.M+numfold*lattice₂.M,lattice₁.M+numfold*lattice₂.M)
                G₂ = ( G₁/lattice₂.b .|> x -> round(x,RoundNearest) ) |> x -> x*lattice₂.b # 1/2は近い偶数になる
                @inbounds for g in 1:numfold
                    from = lattice₁.M+(g-1)*lattice₂.M+1 ; to = lattice₁.M+g*lattice₂.M
                    rmul!( view(Matrix ,1:lattice₁.M,from:to) , phase( view(-G₂,g,:) , lattice₂.τ ) )
                    lmul!( phase( view(G₁,g,:) , lattice₁.τ ) , view(Matrix ,1:lattice₁.M,from:to) )
                end
                @. Matrix += Matrix'
                from  = 1; to = lattice₁.M
                Matrix[from:to,from:to] = lattice₁.Ham!(kx,ky)
                @inbounds for g in 1:numfold
                    from = lattice₁.M+(g-1)*lattice₂.M+1 ; to = lattice₁.M+g*lattice₂.M
                    Matrix[from:to,from:to] = lattice₂.Ham!(kx+G₁[g,1],ky+G₁[g,2])
                end
                Matrix
            end
            return new(_bilayer)
        end
    end
end
