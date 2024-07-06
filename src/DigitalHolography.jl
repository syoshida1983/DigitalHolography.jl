module DigitalHolography

export PSDH2
export PSDH3
export PSDH4
export ParallelPSDH
export GeneralizedPSDH

"""
    PSDH2(I₁, I₂, Iᵣ, δ)

return object wave extracted by the two-step phase-shifting method (see Ref. 1).

# Arguments
- `I₁`, `I₂`: Interferograms corresponding to phases ``\\phi`` and ``\\phi - \\delta``, respectively.
- `Iᵣ`: intensity of reference wave.
- `δ`: phase difference ``\\delta``.

> 1. [X. F. Meng, L. Z. Cai, X. F. Xu, X. L. Yang, X. X. Shen, G. Y. Dong, and Y. R. Wang, "Two-step phase-shifting interferometry and its application in image encryption," Opt. Lett. **31**, 1414-1416 (2006)](https://doi.org/10.1364/OL.31.001414)
"""
function PSDH2(I₁, I₂, Iᵣ, δ)
    u = 2(1 - cos(δ))
    v = @. 2(1 - cos(δ))*(I₁ + I₂) + 4Iᵣ*sin(δ)^2
    w = @. I₁^2 + I₂^2 - 2*I₁*I₂*cos(δ) + 4Iᵣ^2*sin(δ)^2
    a = @. (v - √(v^2 - 4u*w + 0im))/(2u)
    return @. (I₁ - a)/(2√Iᵣ) + im*(I₂ - I₁*cos(δ) - (1 - cos(δ))a)/(2sin(δ)√Iᵣ)
end

"""
    PSDH3(I₁, I₂, I₃, δ)

return object wave extracted by the three-step phase-shifting method (see Ref. 1, 2).
Phases of `I₁`, `I₂`, and `I₃` correspond to ``\\phi - \\delta``, ``\\phi``, and ``\\phi + \\delta``, respectively.

> 1. [Katherine Creath, "V Phase-Measurement Interferometry Techniques," Progress in Optics **26**, 349-393 (1988)](https://doi.org/10.1016/S0079-6638(08)70178-1)
> 2. [Horst Schreiber, and John H. Bruning, "Phase Shifting Interferometry," in Daniel Malacara (ed.), *Optical Shop Testing*, John Wiley & Sons, Ltd, pp. 547-666 (2006)](https://doi.org/10.1002/9780470135976.ch14)
"""
function PSDH3(I₁, I₂, I₃, δ)
    ϕ = @. atan((I₁ - I₃)/sin(δ), (2I₂ - I₁ - I₃)/(1 - cos(δ)))
    v = @. (I₁ + I₃ - 2I₂*cos(δ))/(2(1 - cos(δ)))
    w = @. √(((1 - cos(δ))*(I₁ - I₃))^2 + (sin(δ)*(2I₂ - I₁ - I₃))^2)/(2sin(δ)*(1 - cos(δ)))
    return @. √((v - √(v^2 - w^2 + 0im))/2)*exp(im*ϕ)
end

"""
    PSDH4(I₁, I₂, I₃, I₄)

return object wave extracted by the four-step phase-shifting method. See (Creath, 1988), (Schreiber and Bruning, 2006).
Phase difference `δ` corresponding to `I₁`, `I₂`, `I₃`, and `I₄` are assumed to be ``\\delta = 0, \\dfrac{\\pi}{2}, \\pi``, and ``\\dfrac{3\\pi}{2}``, respectively.
"""
function PSDH4(I₁, I₂, I₃, I₄)
    ϕ = @. atan(I₂ - I₄, I₁ - I₃)
    v = @. (I₁ + I₂ + I₃ + I₄)/4
    w = @. √((I₁ - I₃)^2 + (I₂ - I₄)^2)/2
    return @. √((v - √(v^2 - w^2 + 0im))/2)*exp(im*ϕ)
end

"""
    ParallelPSDH(I)

return object wave extracted by the parallel four-step phase-shifting method (see Ref. 1, 2).
Using 2x2 pixels with phase differences of ``\\dfrac{\\pi}{2}`` each as units, the object wave is extracted through the four-step phase-shifting method.

> 1. [Y. Awatsuji, M. Sasada, and T. Kubota, "Parallel quasi-phase-shifting digital holography," Appl. Phys. Lett., **85**, 1069-1071 (2004).](https://doi.org/10.1063/1.1777796)
> 2. [Yasuhiro Awatsuji, "Parallel Phase-Shifting Digital Holography," in Bahram Javidi, Enrique Tajahuerce, Pedro Andrés (eds.), *Multi-Dimensional Imaging*, John Wiley & Sons, Ltd, pp. 1-23](https://doi.org/10.1002/9781118705766.ch1)
"""
function ParallelPSDH(I)
    Ny, Nx = size(I)
    u = Array{ComplexF64}(undef, Ny, Nx)

    @fastmath @inbounds for j ∈ 1:Nx - 1, i ∈ 1:Ny - 1
        I₁ = I[i+(i+1)%2, j+j%2]        # +cos
        I₂ = I[i+(i+1)%2, j+(j+1)%2]    # +sin
        I₃ = I[i+i%2, j+(j+1)%2]        # -cos
        I₄ = I[i+i%2, j+j%2]            # -sin
        ϕ = atan(I₂ - I₄, I₁ - I₃)
        v = (I₁ + I₂ + I₃ + I₄)/4
        w = √((I₁ - I₃)^2 + (I₂ - I₄)^2)/2
        u[i, j] = √((v - √(v^2 - w^2 + 0im))/2)*exp(im*ϕ)
    end

    u[:,end] = u[:,end-1]
    u[end,:] = u[end-1,:]
    u[end,end] = u[end-1,end-1]

    return u
end

"""
    GeneralizedPSDH(I)

return object wave extracted by the generalized phase-shifting method. See (Creath, 1988), (Schreiber and Bruning, 2006).
Store N interferograms in a three-dimensional array as `I[:,:,1]`, `I[:,:,2]`, ..., `I[:,:,N]`.
Phase difference ``\\delta_{n}`` corresponding to the n-th interferogram is assumed to be ``\\delta_{n} = \\dfrac{2\\pi}{N}n``.
"""
function GeneralizedPSDH(I)
    N = size(I, 3)
    x = sum(I[:,:,i].*cos(i*2π/N) for i in 1:N)
    y = sum(I[:,:,i].*sin(i*2π/N) for i in 1:N)
    ϕ = @. atan(-y, x)
    v = sum(I[:,:,i] for i in 1:N)/N
    w = @. 2√(x^2 + y^2)/N
    return @. √((v - √(v^2 - w^2 + 0im))/2)*exp(im*ϕ)
end

end
