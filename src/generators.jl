###
### Abstract
###

abstract type AbstractGenerator{I, O} end

Base.eltype(::AbstractGenerator{<:Any, O}) where O = NTuple{O, Float64}
in_dim(::AbstractGenerator{I}) where I = I
out_dim(::AbstractGenerator{<:Any, O}) where O = O
weight(::AbstractGenerator) = 1

Base.rand(g::AbstractGenerator{I}) where I = g(rand(I)...)

"""
    sample(generator, n[; origin=ntuple(zero, I), svec=false, random=true])

Sample `n` points from generator, translated by `origin`. If `svec` is set to true, return
`SVector`s insted of `Tuple`s. If `random` is true, sample randomly, otherwise sample
by splitting the input into equal parts. This will not yield completely equidistant data.
"""
function StatsBase.sample(
    gen::AbstractGenerator{<:Any, I}, n;
    origin=ntuple(zero, I), svec=false, random=true
) where I
    if random
        points = sample_rand(gen, n)
    else
        points = sample_regular(gen, n)
    end

    if svec
        return [SVector{I, T}.(p .+ origin) for p in points]
    else
        return [p .+ origin for p in points]
    end
end

function sample_rand(gen::AbstractGenerator{I}, n) where I
    points = Vector{eltype(gen)}(undef, n)
    for i in 1:n
        input = ntuple(_ -> rand(), I)
        n_iters = 0
        while (p = rand(gen)) === nothing
            n_iters += 1
            if n_iters ≥ 10_000
                error("no point found after $n_iters iterations")
            end
        end
        points[i] = p
    end
    return points
end

@generated function sample_regular(gen::AbstractGenerator{I}, n) where I
    quote
        k = round(Int, n^(1//$I))
        points = eltype(gen)[]
        i = 1
        @nloops $I j (_ -> range(0, 1, length=k+1)[1:end-1]) begin
            input = @ncall $I tuple j
            p = gen(input...)
            p !== nothing && push!(points, p)
        end
        return points
    end
end

###
### Basic spaces
###

"""
    Sphere{I, O} <: AbstractGenerator{I, O}

`D`-dimensional sphere with radius `r`.

# Constructor

`Sphere(D[; r=1])`
"""
struct Sphere{I, O} <: AbstractGenerator{I, O}
    r::Float64
end
Sphere(I; r=1) = Sphere{I, I+1}(r)
Sphere(I, r) = Sphere{I, I+1}(r)

function Base.rand(s::Sphere{I}) where I
    vec = normalize!(randn(I + 1))
    return s.r .* tuple(vec...)
end

function (s::Sphere{I, O})(args::Vararg{<:Any, I}) where {I, O}
    first = ntuple(I) do i
        x_i = cos(2π * args[i])
        for j in i-1:-1:1
            x_i *= s.r * sin(2π * args[j])
        end
        x_i
    end
    last = s.r * prod(t -> sin(2π * t), args)
    tuple(first..., last)
end

"""
    Torus <: AbstractGenerator{2, 3}

Torus with outer radius `r1` and inner radius `r2`.

# Constructor

`Torus([; r1=2, r2=1])`
"""
@kwdef struct Torus <: AbstractGenerator{2, 3}
    r1::Float64 = 2
    r2::Float64 = 1
end

function (t::Torus)(θ, φ)
    r1 = t.r1
    r2 = t.r2
    θ *= 2π
    φ *= 2π
    return ((r1 + r2*cos(θ)) * cos(φ), (r1 + r2*cos(θ)) * sin(φ), r2*sin(θ))
end

"""
    Knot <: AbstractGenerator{1, 3}

Torus knot parametrized as:

```math
  t \\mapsto (m \\cos(p t) + n \\cos((2 - q) t),
             m \\sin(p t) + n \\sin((2 - q) t),
             h \\sin(q t))
```

# Constructor

`Knot([; p=2, q=3, n=1, m=1.5, h=1])`
"""
@kwdef struct Knot <: AbstractGenerator{1, 3}
    p::Float64 = 2
    q::Float64 = 3
    n::Float64 = 1
    m::Float64 = 1.5
    h::Float64 = 1
end

function (k::Knot)(t)
    t *= 2π
    p, q, n, m, h = k.p, k.q, k.n, k.m, k.h
    return (m * cos(p * t) + n * cos((2 - q) * t),
            m * sin(p * t) + n * sin((2 - q) * t),
            h * sin(q * t))
end

"""
    Klein <: AbstractGenerator{2, 4}

Klein bottle embedded in ``\\mathbb{R}^4``. Parametrized as:

```math
(\\vartheta, \\varphi) \\mapsto ((r1 + r2 \\cos(\\vartheta)) \\cos(\\varphi),
                              (r1 + r2 \\cos(\\vartheta)) \\sin(\\varphi),
                              r2 \\sin(\\vartheta) \\cos(\\frac{\\varphi}{2}),
                              r2 \\sin(\\vartheta) \\sin(\\frac{\\varphi}{2}))
```

# Constrcutor

`Klein([; r1=1, r2=2])`
"""
@kwdef struct Klein <: AbstractGenerator{2, 4}
    r1::Float64 = 1
    r2::Float64 = 1
end

function (k::Klein)(θ, φ)
    r1, r2 = k.r1, k.r2
    θ *= 2π
    φ *= 2π
    return ((r1 + r2 * cos(θ)) * cos(φ),
            (r1 + r2 * cos(θ)) * sin(φ),
            r2 * sin(θ) * cos(φ/2),
            r2 * sin(θ) * sin(φ/2))
end

###
### Subtractable
###

"""
    Cube{I} <: AbstractGenerator{I, I}

A `D`-dimensional hypercube. Dimensions determined by arguments. Can be subtracted from
other generators.

# Constructor

`Cube(D[, dimensions...])`
"""
@kwdef struct Cube{I} <: AbstractGenerator{I, I}
    rs::NTuple{I, Float64}
end
Cube(I, args...) = Cube{I}(args)
Cube(I) = Cube{I}(ntuple(_ -> 1.0, I))

function (c::Cube{I})(args::Vararg{<:Any, I}) where I
    return args
end

function Base.in(v, c::Cube{I}) where I
    if length(v) > I
        error("dimnension to high to be in Ball{$I}.")
    else
        return all((v .> 0) .& (v .< c.rs))
    end
end

"""
    Ball{I} <: AbstractGenerator{I, I}

`D`-dimensional ball with radius `r`. Can be subtracted from other generators.

# Constructor

`Ball(D[; r=1])`
"""
struct Ball{I} <: AbstractGenerator{I, I}
    r::Float64
end
Ball(I; r=1) = Ball{I}(r)
Ball(I, r) = Ball{I}(r)

function Base.in(v, b::Ball{I}) where I
    if length(v) > I
        error("dimnension to high to be in Ball{$I}.")
    else
        return norm(SVector(v)) ≤ b.r
    end
end

function Base.rand(b::Ball{I}) where I
    while true
        p = ntuple(_ -> b.r * (2rand() - 1), I)
        if p in b
            return p
        end
    end
end

###
### Modifiers
###

"""
    Noisy{I, O, A<:AbstractGenerator{I, O}} <: AbstractGenerator{I, O}

Add some noise to a generator.

# Constructor

`Noisy(generator[; amount=0.1, gauss=false])`
"""
struct Noisy{I, O, A<:AbstractGenerator{I, O}} <: AbstractGenerator{I, O}
    generator::A
    amount::Float64
    gauss::Bool
end
Noisy(gen; amount=0.1, gauss=false) = Noisy(gen, amount, gauss)

function (n::Noisy{I, O})(args::Vararg{<:Any, I}) where {I, O}
    p = n.generator(args...)
    if n.gauss
        return p .+ n.amount .* tuple(randn(O)...)
    else
        return p .+ n.amount .* tuple(rand(O)...)
    end
end

###
### Combinations
###

"""
    Product{I, O, A<:AbstractGenerator, B<:AbstractGenerator} <: AbstractGenerator{I, O}

Cartesian product of two generators.

# Example

```jldoctest
gen = Sphere(1) × Sphere(1)
gen(0, 0.2)

# output

(1.0, 0.0, 0.30901699437494745, 0.9510565162951535)
```
"""
struct Product{I, O, A<:AbstractGenerator, B<:AbstractGenerator} <: AbstractGenerator{I, O}
    left::A
    right::B
end
function Product(
    left::AbstractGenerator{I1, O1}, right::AbstractGenerator{I2, O2}
) where {I1, I2, O1, O2}
    return Product{I1 + I2, O1 + O2, typeof(left), typeof(right)}(left, right)
end
LinearAlgebra.:×(left::AbstractGenerator, right::AbstractGenerator) = Product(left, right)

Base.show(io::IO, p::Product) = print(io, "($(p.left) × $(p.right))")

function (p::Product{I, O})(args::Vararg{<:Any, I}) where {I, O}
    left_i = in_dim(p.left)
    l = p.left(args[1:left_i]...)
    r = p.right(args[left_i+1:end]...)
    return tuple(l..., r...)
end

"""
    Difference{I, O, A<:AbstractGenerator, B<:AbstractGenerator} <: AbstractGenerator{I, O}

Difference between two generators. Point in the first generator are dropped if they are `in`
the second generator. Constructed with `-`.

# Limitations

* Currently, only `Ball` and `Cube` can be subtracted.
* `Ball` or `Cube` dimension must match output dimension of generator.
* When sampling non-randomly, fewer than `n` points may be sampled.
"""
struct Difference{I, O, A<:AbstractGenerator, B<:AbstractGenerator} <: AbstractGenerator{I, O}
    left::A
    right::B
end
function Difference(
    left::AbstractGenerator{I, O}, right::AbstractGenerator{O, O}
) where {I, O}
    return Difference{I, O, typeof(left), typeof(right)}(left, right)
end
Base.:-(left::AbstractGenerator, right::AbstractGenerator) = Difference(left, right)

Base.show(io::IO, d::Difference) = print(io, "($(d.left) - $(d.right))")

function (m::Difference{I})(args::Vararg{<:Any, I}) where I
    p = m.left(args...)
    if p in m.right
        return nothing
    else
        return p
    end
end

function Base.rand(m::Difference)
    p = rand(m.left)
    if p in m.right
        return nothing
    else
        return p
    end
end

"""
    Sum{I, O, T} <: AbstractGenerator{I, O}

Sample points from both generators. Constructed with `+`.

# Limitations

* Sampling will draw an equal number of points from all spaces added together even if spaces
  are of different "sizes".
* Generators must be translated manually or spaces will overlap.
"""
struct Sum{I, O, T} <: AbstractGenerator{I, O}
    gens::T
end
function Sum(gens)
    I = maximum(in_dim.(gens)) + 1
    O = maximum(out_dim.(gens))
    Sum{I, O, typeof(gens)}(gens)
end
Base.:+(left::AbstractGenerator, right::AbstractGenerator) = Sum((left, right))
Base.:+(u::Sum, right::AbstractGenerator) = Sum((u.gens..., right))
Base.:+(left::AbstractGenerator, u::Sum) = Sum((left, u.gens...))

Base.show(io::IO, s::Sum) = print(io, "(", join(string.(s.gens), " + "), ")")

function (u::Sum{I, O})(args::Vararg{Any, I}) where {I, O}
    weights = weight.(u.gens)
    selector = args[1] * sum(weights)
    gen = u.gens[searchsortedlast(cumsum([weights...]), selector) + 1]
    return tuple(gen(args[2:in_dim(gen)+1]...)..., ntuple(_ -> 0.0, O - out_dim(gen))...)
end

function Base.rand(u::Sum{<:Any, O}) where O
    weights = weight.(u.gens)
    selector = rand() * sum(weights)
    gen = u.gens[searchsortedlast(cumsum([weights...]), selector) + 1]
    return tuple(rand(gen)..., ntuple(_ -> 0.0, O - out_dim(gen))...)
end

"""
    Translated{I, O, G} <: AbstractGenerator{I, O}

Move space around by adding (with `+`) a `Vector` or `Tuple` of numbers to it.
"""
struct Translated{I, O, G} <: AbstractGenerator{I, O}
    gen::G
    offset::NTuple{O, Float64}
end
function Translated(gen::AbstractGenerator{I}, offset) where {I}
    O = max(out_dim(gen), length(offset))
    off = (offset..., ntuple(_ -> 0.0, Val(O - length(offset))))
    return Translated{I, O, typeof(gen)}(gen, offset)
end
Base.:+(gen::AbstractGenerator, off::Union{AbstractVector, Tuple}) = Translated(gen, off)
Base.:+(off::Union{AbstractVector, Tuple}, gen::AbstractGenerator) = Translated(gen, off)

Base.show(io::IO, t::Translated) = print(io, "($(t.gen) + $(t.offset))")

function (t::Translated{I, O})(args::Vararg{Any, I}) where {I, O}
    return tuple(t.gen(args...)..., ntuple(_ -> 0.0, O - out_dim(t.gen))...) .+ t.offset
end

function Base.rand(t::Translated{I, O}) where {I, O}
    return tuple(rand(t.gen)..., ntuple(_ -> 0.0, O - out_dim(t.gen))...) .+ t.offset
end

#TODO: rotate, warp
# weight, + : select space with weight probability.
# a la 2*Sphere(1, 2) + Sphere(1) -- select twice as many points from first sphere.

function projection(points, to=3)
    T = eltype(first(points))
    n = length(first(points))
    if n ≤ to
        return points
    else
        points_matrix = reshape(reinterpret(T, points), (n, length(points)))
        proj = transform(fit(PCA, points_matrix), points_matrix)
        return [ntuple(i -> proj[i, j], to) for j in 1:length(points)]
    end
end

"""
   projplot

Plot collection of `Tuple`s or `SVector`s as points in euclidean space with sane
defaults. If dimension is higer than `dim` keyword argument (defaults to 3), project points
to `dim` dimensions using PCA and color points according to first dimension not shown.
"""
projplot

@userplot ProjPlot
@recipe function f(arg::ProjPlot)
    pts = arg.args[1]

    @series begin
        dim = length(first(pts))
        max_dim = get(plotattributes, :dim, 3)
        if dim ≥ max_dim + 1
            proj = projection(pts, max_dim + 1)
            pts = [p[1:max_dim] for p in proj]
            marker_z --> [p[max_dim + 1] for p in proj]
            markerstrokecolor --> :auto
            markersize --> 2
            marker --> :heat
        else
            markersize --> (dim ≤ 2 ? 2 : 1)
        end
        seriestype := :scatter

        pts
    end
end
