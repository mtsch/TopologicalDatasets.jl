struct GeodesicMetric{T, D}
    dists::Matrix{T}
    graph::SimpleWeightedGraph{Int64, T}
    points::Vector{NTuple{D, T}}
end

function GeodesicMetric(points::Vector, r1::Real, r2::Real; metric=Euclidean())
    r1 > r2 && @warn "you probably want r1 ≤ r2"
    points = shuffle(points)
    tree = KDTree(SVector.(points), metric)

    # Create a subsample of points r1 apart.
    visited = falses(length(points))
    new_points = eltype(points)[]
    for i in eachindex(points)
        visited[i] && continue
        push!(new_points, points[i])
        visited[inrange(tree, SVector(points[i]), r1)] .= true
    end
    points = new_points

    I = Int[]
    J = Int[]
    V = typeof(metric(SVector(points[1]), SVector(points[2])))[]
    tree = KDTree(SVector.(points), metric)
    for i in eachindex(new_points)
        for j in inrange(tree, SVector(points[i]), 2r2)
            d = metric(SVector(points[i]), SVector(points[j]))
            append!(I, (i, j))
            append!(J, (j, i))
            append!(V, (d, d))
        end
    end
    graph = SimpleWeightedGraph(sparse(I, J, V, length(new_points), length(new_points)))
    dists = floyd_warshall_shortest_paths(graph).dists

    return GeodesicMetric(dists, graph, points)
end

function Base.show(io::IO, gm::GeodesicMetric{T}) where T
    print(io, "GeodesicMetric{$T}")
end
function Base.show(io::IO, ::MIME"text/plain", gm::GeodesicMetric)
    print(io, gm)
    print(io, " on $(length(gm.points)) points with $(ne(gm.graph)) edges")
end

distances(gm::GeodesicMetric) = gm.dists
points(gm::GeodesicMetric) = gm.points

@recipe function f(gm::GeodesicMetric)
    graph = get(plotattributes, :graph, true)
    @series begin
        seriestype := :path
        label --> "graph"

        xs = Float64[]
        ys = Float64[]
        zs = Float64[]
        for e in edges(gm.graph)
            p = gm.points[src(e)]
            q = gm.points[dst(e)]

            length(p) ≥ 1 && append!(xs, (p[1], q[1], NaN))
            length(p) < 2 && append!(ys, (0, 0, NaN))
            length(p) ≥ 2 && append!(ys, (p[2], q[2], NaN))
            length(p) ≥ 3 && append!(zs, (p[3], q[3], NaN))
        end

        if length(gm.points[1]) ≤ 2
            xs, ys
        else
            xs, ys, zs
        end
    end

    @series begin
        seriestype := :scatter
        markersize --> 1
        gm.points
    end
end
