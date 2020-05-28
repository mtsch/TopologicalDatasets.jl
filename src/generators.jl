function sphere(n; dim=3, cutoff=1)
    points = Vector{NTuple{dim, Float64}}(undef, n)
    for i in 1:n
        vec = normalize!(randn(dim))
        while vec[dim] > cutoff
            vec = normalize!(randn(3))
        end
        points[i] = tuple(vec...)
    end

    return points
end
