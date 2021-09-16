using UnPack

export  Bilayer

exppi(x) = cospi(x) + im*sinpi(x)

function phase(G,τ)
    if size(τ,2) == 1
        return exppi(G'*τ)
    else
        return exppi.(G'*τ)
    end
end

function bilayer!(kx,ky,V!,lattice₁,lattice₂,Moire)
    @unpack fold, G₁, G₂, numfold, matrix = Moire
    from  = 1; to = lattice₁.M
    @. matrix = 0.
    G₂ = ( G₁/lattice₂.b .|> x -> round(x,RoundNearest) ) |> x -> x*lattice₂.b # 1/2は近い偶数になる

    @inbounds for g in 1:numfold
        from += lattice₂.M ; to += lattice₂.M
        V!(kx+G₁[g,1],ky+G₁[g,2],view(matrix,1:lattice₁.M,from:to))
        rmul!( view(matrix ,1:lattice₁.M,from:to) , phase( view(-G₂,g,:) , lattice₂.τ ) )
        lmul!( phase( view(G₁,g,:) , lattice₁.τ ) , view(matrix ,1:lattice₁.M,from:to) )
    end
    @. matrix += matrix'
    from  = 1; to = lattice₁.M
    lattice₁.Ham!(kx,ky,view(matrix,from:to,from:to))
    @inbounds for g in 1:numfold
        from += lattice₂.M ; to += lattice₂.M
        lattice₂.Ham!(kx+G₁[g,1],ky+G₁[g,2],view(matrix,from:to,from:to))
    end
end
struct Bilayer
    bilayer!
    function Bilayer(V!,lattice₁,lattice₂,Moire)
        function _bilayer!(kx,ky,V!,lattice₁,lattice₂,Moire,Matrix)
            @unpack fold, G₁, G₂, numfold, matrix = Moire
            from  = 1; to = lattice₁.M
            @. Matrix = 0.
            G₂ = ( G₁/lattice₂.b .|> x -> round(x,RoundNearest) ) |> x -> x*lattice₂.b # 1/2は近い偶数になる
            @inbounds for g in 1:numfold
                from += lattice₂.M ; to += lattice₂.M
                V!(kx+G₁[g,1],ky+G₁[g,2],view(Matrix,1:lattice₁.M,from:to))
                rmul!( view(Matrix ,1:lattice₁.M,from:to) , phase( view(-G₂,g,:) , lattice₂.τ ) )
                lmul!( phase( view(G₁,g,:) , lattice₁.τ ) , view(Matrix ,1:lattice₁.M,from:to) )
            end
            @. Matrix += Matrix'
            from  = 1; to = lattice₁.M
            lattice₁.Ham!(kx,ky,view(Matrix,from:to,from:to))
            @inbounds for g in 1:numfold
                from += lattice₂.M ; to += lattice₂.M
                lattice₂.Ham!(kx+G₁[g,1],ky+G₁[g,2],view(Matrix,from:to,from:to))
            end
        end
        new((kx,ky,Matrix) -> _bilayer!(kx,ky,V!,lattice₁,lattice₂,Moire,Matrix))
    end
end