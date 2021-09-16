export  lattice

struct lattice
    a :: AbstractArray{Real}
    τ :: AbstractArray{Real}
    b :: AbstractArray{Real}
    Ham!
    M :: Int
    function lattice(a,τ,Ham!,M) 
        new(a,τ,2inv(a),Ham!,M)
    end
end