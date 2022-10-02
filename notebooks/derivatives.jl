function truncated(c, k, x)
    if x ≥ c
        (x - c)^k
    else
        zero(x)
    end
end

function truncated(c, k, x, r)
    if x ≥ c && r ≤ k
        -prod((k-r+1):k) * (x - c)^(k-r)
    else
        zero(x)
    end
end


function der(::typeof(truncated), n)
    if n > 0
        (c,k,x) -> truncated(c, k, x, n)
    else
        truncated
    end
end

der(truncated, 3)(2, 2, 3.5)
