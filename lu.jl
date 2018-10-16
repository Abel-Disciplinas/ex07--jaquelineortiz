"""`declu(A; ϵ = 1e-12)`

A função `declu` calcula a fatoração LU de uma matriz densa A sem pivoteamento. L e U
são armazenadas na própria A, com aᵢⱼ = uᵢⱼ se j ≥ i e aᵢⱼ = ℓᵢⱼ se i > j.
Se algum pivô for menor que ϵ, um erro ocorre.
"""
function declu(A, ϵ = 1e-12)
    n = size(A, 1)
    for j = 1:n
        if abs(A[j,j]) < ϵ
            error("Matriz muito próxima de ser singular")
        end
        I = j+1:n
        @views A[I,j] ./= A[j,j]
        @views A[I,I] .-= A[I,j] * A[j,I]'
    end
    return A
end

"""`reslu(A, b)`

Resolve o sistema Ax = b, **substituindo** o b por x.
"""
function reslu(A, b)
    n = length(b)
    for i = 1:n
        for j = 1: i-1
            b[i] -= A[i,j] * b[j]
        end
    end
    for i = n : -1 : 1
        for j = i+1 : n
            b[i] -= A[i,j] * b[j]
        end
        b[i] = b[i] / A[i,i]
    end
    return b
end

function reslu2(A, b)
    n = length(b)
    for j = 1 : n 
        for i = j+1 : n
            b[i] -= A[i,j] * b[j]
        end
    end           
    for i = n : -1 : 1
        s = b[i]
        for j = i+1 : n
            b[i] -= A[i,j] * b[j]
        end
        b[i] = b[i] / A[i,i]
    end
    return b
end

function reslu3(A, b)
    n = length(b)
    for i = 1 : n
        for j = 1 : i-1
            b[i] -= A[i,j] *b[j]
        end
    end
    for j = n : -1 : 1
        b[j] = b[j] / A[j,j]
        for i = j-1 : -1 : 1
            b[i] -= A[i,j] * b[j]
        end
    end
    return b
end

function reslu4(A, b)
    n = length(b)
    for j = 1 : n
        for i = j+1:n
            b[i] -= A[i,j] * b[j]
        end
    end
    for j = n : -1 : 1
        b[j] = b[j] / A[j,j]
        for i = j-1 : -1 : 1
            b[i] -= A[i,j] * b[j]
        end
    end
    return b
end
