function truncated(c, k, x)
    if x ≥ c
        (x - c)^k
    else
        zero(x)
    end
end
