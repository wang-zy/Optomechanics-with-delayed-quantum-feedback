module TwoPhotonSweep

export g0Sweep

using Grid

function Conversion(Theta)
    # Convert a general tensor into 2-D matrix by conbining its first two indexes.
    TensorSize = size(Theta)
    N = length(TensorSize) - 2
    DimX = TensorSize[1] * TensorSize[2]
    DimY = div(prod(TensorSize), DimX)
    reshape(Theta, DimX, DimY), TensorSize, N
end

function MPS_svd(Theta)
    # V is a matrix after tensor contraction
    U, S, V = svd(Theta)
    SvdThreshold = 0.0001
    Dmax = 150

    bond = length(S)
    s = 0.0
    for i in bond:-1:1
        if S[i] < SvdThreshold || bond > Dmax
            bond -= 1
            s += S[i] * S[i]
        else
            break
        end
    end
    Norm = 1.0 / sqrt(1.0-s)
    U[:,1:bond], V'[1:bond,:], S[1:bond]*Norm
end

function Nsvd(Theta)
    # Theta[alpha, i_1, ... , i_n, gamma]
    # return all matrices and vectors directly
    Tmp, TensorSize, N = Conversion(Theta)
    GammaNew = Vector{Array{Complex{Float64},2}}(N)
    LambdaNew = Vector{Array{Float64,1}}(N-1)
    for dim in N:-1:2
        if dim == 2
            SR = MPS_svd(Tmp) # SvdResult
            GammaNew[N] = SR[2]
            GammaNew[N-1] = SR[1]
            LambdaNew[N-1] = SR[3]
        else
            SR = MPS_svd(Tmp)
            GammaNew[N-dim+1] = SR[1]
            LambdaNew[N-dim+1] = SR[3]
            Beta, y = size(SR[2])
            i = TensorSize[3+N-dim]
            DimX = Beta * i
            DimY = div(y, i)
            Tmp = reshape(diagm(SR[3])*SR[2], DimX, DimY)
        end
    end
    GammaNew, LambdaNew, N
end

function LeftMatrixToTensor(X, Lambda)
    # recast the left matrix X into a tensor form
    DimX, DimBeta = size(X)
    DimAlpha = length(Lambda)
    #println("$DimAlpha")
    DimI = div(DimX, DimAlpha)
    for i in eachindex(X)
        X[i] /= Lambda[mod1(i,DimAlpha)]
    end
    reshape(X, DimAlpha, DimI, DimBeta)
end

function RightMatrixToTensor(Y, Lambda)
    # recast the right matrix Y into a tensor form
    DimBeta, DimY = size(Y)
    DimGamma = length(Lambda)
    DimJ = div(DimY, DimGamma)
    for i in eachindex(Y)
        Y[i] /= Lambda[cld(i,DimBeta*DimJ)]
    end
    reshape(Y, DimBeta, DimJ, DimGamma)
end

function Nrecast(Theta, LambdaLeft, LabmbdaRight)
    # Recast a N qubit tensor in MPS form
    # Left: smaller label in the MPS chain
    GammaTmp, LambdaNew, N = Nsvd(Theta)
    GammaNew = Vector{Array{Complex{Float64},3}}(N)
    for i in 1:N
        if i == N
            GammaNew[N] = RightMatrixToTensor(GammaTmp[N], LabmbdaRight)
        elseif i == 1
            GammaNew[1] = LeftMatrixToTensor(GammaTmp[1], LambdaLeft)
        else
            GammaNew[i] = LeftMatrixToTensor(GammaTmp[i], LambdaNew[i-1])
        end
    end
    GammaNew, LambdaNew
end

function GateContraction(U, Gamma, Lambda)
    # Gamma, Lambda: MPS chain -> L G ... L G L
    # < i_1, ..., i_n | U | j_1, ..., j_n > = U[j_1,...,j_n, i_1,...,i_n]
    Size = collect(size(U))
    N = div(length(Size), 2)
    Dim = Size[1:N]
    Beta = zeros(Int, N+1)
    for i in 1:N+1
        Beta[i] = length(Lambda[i])
    end
    tmp::Array{Int32, 1} = [Beta[1]; Dim; Beta[N+1]]
    Theta = Array{Complex{Float64}}(tmp...)
    #Theta = Array{Complex{Float64}}(Beta[1], Dim..., Beta[N+1])

    for gamma in 1:Beta[N+1], I in Counter(Dim), alpha in 1:Beta[1]
        for J in Counter(Dim), B in Counter(collect(Beta[2:N]))
            t::Complex{Float64} = U[J...,I...] * Gamma[1][alpha,J[1],B[1]] * Gamma[N][B[N-1],J[N],gamma] * Lambda[2][B[1]]
            for k = 2:N-1
                t *= Gamma[k][B[k-1], J[k], B[k]] * Lambda[k+1][B[k]]
            end
            Theta[alpha,I...,gamma] += t
        end
        Theta[alpha,I...,gamma] *= Lambda[1][alpha] * Lambda[N+1][gamma]
    end
    Theta
end

function TwoGateContraction(U, Gamma, Lambda)
    # Specific version of GateContraction with N = 2
    # To improve performance by ensuring type stability
    # May be replaced by the general version later
    Size = size(U)
    Dim1 = Size[1]
    Dim2 = Size[2]
    DimAlpha = length(Lambda[1])
    DimGamma = length(Lambda[3])
    Beta = length(Lambda[2])
    Theta = zeros(Complex{Float64}, DimAlpha, Dim1, Dim2, DimGamma)

    for gamma in 1:DimGamma, I2 in 1:Dim2, I1 in 1:Dim1, alpha in 1:DimAlpha
        for J2 in 1:Dim2, J1 in 1:Dim1, B in 1:Beta
            Theta[alpha,I1,I2,gamma] += U[J1,J2,I1,I2] * Gamma[1][alpha,J1,B] * Gamma[2][B,J2,gamma]* Lambda[2][B]
        end
        Theta[alpha,I1,I2,gamma] *= Lambda[1][alpha] * Lambda[3][gamma]
    end
    Theta
end

function NQubitGate(U, Gamma, Lambda)
    #Theta = GateContraction(U, Gamma, Lambda)
    Theta = TwoGateContraction(U, Gamma, Lambda)
    Gamma, Lambda = Nrecast(Theta, Lambda[1], Lambda[length(Lambda)])
    Gamma, Lambda[1]
end

function SwapGate(Gamma, Lambda)
    DimAlpha, DimI, DimBeta = size(Gamma[1])
    DimJ = size(Gamma[2])[2]
    DimGamma = size(Gamma[2])[3]
    Theta = zeros(Complex{Float64}, DimAlpha * DimJ, DimI * DimGamma)
    for gamma in 1:DimGamma, i in 1:DimI, j in 1:DimJ, alpha in 1:DimAlpha
        s = zero(Complex{Float64})
        for beta in 1:DimBeta
            s += Gamma[1][alpha,i,beta] * Lambda[2][beta] * Gamma[2][beta,j,gamma]
        end
        Theta[alpha+(j-1)*DimAlpha, i+(gamma-1)*DimI] = s * Lambda[1][alpha] * Lambda[3][gamma]
    end
    X, Y, NewLambda = MPS_svd(Theta)
    tmp = Vector{Array{Complex{Float64},3}}(2)
    tmp[1] = LeftMatrixToTensor(X, Lambda[1])
    tmp[2] = RightMatrixToTensor(Y, Lambda[3])
    tmp, NewLambda
end

# Convert the physical parameters into MPS Simulation parameters.
type MPS_parameter
    rdt::Float64 # gamma * dt
    wdt::Float64 # omega_m * dt
    gdt::Float64 # g_0 * dt

    length::Int # related to total time of the simulation

    DIM_A::Int # dimension of the optical cavity
    DIM_S::Int # dimension of the mechanical system
    DIM_B::Int # dimension of the time bins

    width::Int # photon linewidth in time domain
    Tm::Int
    TT::Float64

    # infomation of the input photons
    vacuum::Float64 # probablity of |0> componenet
    Bounce::Int # whether or not to include the feedback process

    function MPS_parameter(;gamma=100, omega_m=10, g_0=20, t_total=1, N=10, dim_a=3, dim_s=15, dim_b=3, vacuum=0.5, widthfactor = 2.0, NumBounce=1)
        dt = 1 / N
        G = max(gamma, omega_m, g_0)
        Rdt = dt * gamma / G
        Wdt = dt * omega_m / G
        Gdt = dt * g_0 / G

        T_m = round(Int, 2 * pi / (omega_m / G) * N)
        Length = round(Int, T_m * t_total)
        Width = round(Int, 1 / gamma * G * N * widthfactor)

        new(Rdt, Wdt, Gdt, Length, dim_a, dim_s, dim_b, Width, T_m, t_total, vacuum, NumBounce)
    end
end

function Annihilation(D)
    # return the Annihilation operator with dimension D
    tmp = zeros(Complex{Float64}, D, D)
    for i in 1:D-1
        tmp[i, i+1] = Complex{Float64}(sqrt(i), 0)
    end
    tmp
end

function OperatorToTensor(U, Dim...)
    # Convert operator (kron product) to tensor
    # U[(i1-1)D2..Dn+(i2-1)D3..Dn+..+in, (j1-1)D2..Dn+(j2-1)D3..Dn+..+jn] =
    # < i_1, ..., i_n | U | j_1, ..., j_n > = U[j_1,...,j_n, i_1,...,i_n]
    # Dim is array of dimension of each site
    Utensor = zeros(Complex{Float64}, Dim..., Dim...)
    N = length(Dim)
    for I in Counter(Dim), J in Counter(Dim)
        X::Int32 = I[N]
        Y::Int32 = J[N]
        for k in 1:N-1
            X += (I[k]-1) * prod(Dim[k+1:N])
            Y += (J[k]-1) * prod(Dim[k+1:N])
        end
        Utensor[J...,I...] = U[X, Y]
    end
    Utensor
end

function TwoOperator(p::MPS_parameter)
    # U = exp[ -i*ω_m*dt*b^†*b + i*g_0*dt*a^†*a*(b+b^†) + sqrt(γdt/2)*(e^{iϕ}B(k-l)^†*a-e^{-iϕ}B(k-l)*a^†)
    #          + sqrt(γdt/2)*(B(k)^†*a-B(k)*a^†)]
    # Order: left(k-l), middle(system(opto before mech)), right(k)
    DIM_B = p.DIM_B
    DIM_A = p.DIM_A
    DIM_S = p.DIM_S
    a = kron(Annihilation(DIM_A), eye(DIM_S))
    b = kron(eye(DIM_A), Annihilation(DIM_S))
    HS = -im*p.wdt*b'*b - im*p.gdt*a'*a*(b'+b)
    B = Annihilation(DIM_B)
    H = kron(HS, eye(DIM_B)) + sqrt(p.rdt) * (kron(a, B') - kron(a', B))
    OperatorToTensor(expm(H), DIM_A*DIM_S, DIM_B)
end

function MPS_chain(p::MPS_parameter)
    # generate the MPS chain for simulation
    system = zeros(Complex{Float64}, 1, p.DIM_A*p.DIM_S, 1)
    system[1, 1, 1] = Complex{Float64}(1, 0)
    bins = zeros(Complex{Float64}, 1, p.DIM_B, 1)
    bins[1,1,1] = Complex{Float64}(1, 0)

    Gamma = fill(bins, p.length + 1)
    Gamma[1] = system
    Lambda = fill([1.0], p.length + 2)
    Gamma, Lambda
end

function TotalDM(Gamma, LambdaLeft, LambdaRight, DIM_A, DIM_S)
    # return the density matrix of the optomechanical system
    DimSystem = DIM_A * DIM_S
    Left = length(LambdaLeft)
    Right = length(LambdaRight)
    DM = zeros(Complex{Float64}, DimSystem, DimSystem)
    for x in 1:Left, j in 1:DimSystem, i in 1:DimSystem, y in 1:Right
        DM[i,j] += LambdaLeft[x]^2 * LambdaRight[y]^2 * Gamma[x,i,y] * Gamma[x,j,y]'
    end
    DM
end

function DensityMatrix(Gamma, LambdaLeft, LambdaRight, DIM_A, DIM_S)
    # return the density matrix of mechanical and optical system respectly
    T_DM = TotalDM(Gamma, LambdaLeft, LambdaRight, DIM_A, DIM_S)
    M_DM = zeros(Complex{Float64}, DIM_S, DIM_S)
    O_DM = zeros(Complex{Float64}, DIM_A, DIM_A)
    for k in 1:DIM_A, j in 1:DIM_S, i in 1:DIM_S
        M_DM[i,j] = M_DM[i,j] + T_DM[(k-1)*DIM_S+i, (k-1)*DIM_S+j]
    end
    for j in 1:DIM_A, i in 1:DIM_A, k in 1:DIM_S
        O_DM[i,j] = O_DM[i,j] + T_DM[(i-1)*DIM_S+k, (j-1)*DIM_S+k]
    end
    M_DM, O_DM
end

function extract(Gamma, Lambda, L, p::MPS_parameter)
    # extract information from the MPS chain
    # L is the label of the system
    DensityMatrix(Gamma[L], Lambda[L], Lambda[L+1], p.DIM_A, p.DIM_S)
end

function OneSiteAve(Gamma, Lambda, C)
    # Note that here Lambda is still an array (of length 2), while Gamma is not
    DimAlpha, Dim, DimBeta = size(Gamma)
    s = zero(Complex{Float64})
    for beta in 1:DimBeta, i in 1:Dim, i_ in 1:Dim, alpha in 1:DimAlpha
        s += C[i_, i] * (Lambda[1][alpha] * Gamma[alpha, i_, beta] * Lambda[2][beta])' * Lambda[1][alpha] * Gamma[alpha, i, beta] * Lambda[2][beta]
    end
    real(s)
end

function AveNum(dm)
    # return the average number from input density matrix
    Dim = size(dm)[1]
    s = 0.0
    for i in 2:Dim
        s += real(dm[i,i]) * (i-1)
    end
    s
end

function LineShape(width)
    # generate the line shape of a light pulse
    line = zeros(width)
    for i in 1:width
        #line[i] = exp(-(i-(width+1)/2)^2 / width^2 * 16)
        line[i] = exp(-(i-(width+1)/2)^2 / width^2 * 16)
    end
    s = 0.0
    for i in 1:width
        s += line[i]^2
    end
    for i in 1:width
        line[i] = line[i] / sqrt(s)
    end
    line
end

function DMGeneration(lineshape, vacuum, k)
    # return the density matrix for part of a photon, called by the OnePhoton function
    lineshape = lineshape * sqrt(1 - vacuum)
    lineshape = [sqrt(vacuum); lineshape]
    RhoK = zeros(k+1, k+1)
    tmp = 0.0
    for i in 1:k+1
        tmp += lineshape[i] * lineshape[i]
        for j in 1:k+1
            RhoK[i, j] = lineshape[i] * lineshape[j]
        end
    end
    RhoK[1, 1] = vacuum + 1 - tmp
    RhoK
end

function OnePhoton(lineshape, vacuum, DimB)
    # generate the MPS chain for single photon state with a given lineshape
    N = length(lineshape) # total length of the chain
    Gamma = Vector{Array{Complex{Float64},3}}(N)
    Lambda = Vector{Array{Float64,1}}(N-1)
    EigThreshold = 1e-7
    StartDict = Dict{Float64, Array{Float64,1}}()
    EndDict = Dict{Float64, Array{Float64,1}}()
    for k in 1:N
        if k == 1
            RhoK = DMGeneration(lineshape, vacuum, k)
            E, V = eig(RhoK)
            for i in 1:length(E)
                StartDict[E[i]] = V[:, i]
            end
            for key in keys(StartDict)
                if key < EigThreshold
                    delete!(StartDict, key)
                end
            end
            Lambda[1] = sort(collect(keys(StartDict)), rev = true)
            gamma = zeros(Complex{Float64}, 1, DimB, length(Lambda[1]))
            for alpha in 1:length(Lambda[1])
                gamma[1, 1, alpha] = StartDict[Lambda[1][alpha]][1]
                gamma[1, 2, alpha] = StartDict[Lambda[1][alpha]][2]
            end
            Gamma[1] = gamma
        elseif k == N
            gamma = zeros(Complex{Float64}, length(Lambda[N-1]), DimB, 1)
            lineshape = lineshape * sqrt(1 - vacuum)
            lineshape = [sqrt(vacuum); lineshape]
            for alpha in 1:length(Lambda[N-1])
                tmp = StartDict[Lambda[N-1][alpha]]
                gamma[alpha, 2, 1] = lineshape[N+1] * tmp[1] / sqrt(Lambda[N-1][alpha])
                tt = lineshape[1:N]
                gamma[alpha, 1, 1] = dot(tt, tmp) / sqrt(Lambda[N-1][alpha])
            end
            Gamma[N] = gamma
        else
            RhoK = DMGeneration(lineshape, vacuum, k)
            E, V = eig(RhoK)
            for i in 1:length(E)
                EndDict[E[i]] = V[:, i]
            end
            for key in keys(EndDict)
                if key < EigThreshold
                    delete!(EndDict, key)
                end
            end
            Lambda[k] = sort(collect(keys(EndDict)), rev = true)
            gamma = zeros(Complex{Float64}, length(Lambda[k-1]), DimB, length(Lambda[k]))
            for alpha in 1:length(Lambda[k]), beta in 1:length(Lambda[k-1])
                tmp1 = StartDict[Lambda[k-1][beta]]
                tmp2 = EndDict[Lambda[k][alpha]]
                gamma[beta, 2, alpha] = tmp1[1] * tmp2[length(tmp2)] / sqrt(Lambda[k-1][beta])
                tt = tmp2[1:length(tmp2)-1]
                gamma[beta, 1, alpha] = dot(tmp1, tt) / sqrt(Lambda[k-1][beta])
            end
            Gamma[k] = gamma
            StartDict = deepcopy(EndDict)
            EndDict = Dict{Float64, Array{Float64,1}}()
        end
    end
    for i in 1:N-1
        Lambda[i] = sqrt(Lambda[i])
    end
    Gamma, Lambda
end

function EECal(Lambda)
    # calculate the EE1 between two photons
    # and also EE2 between two photons and optomechanical system
    pos1 = round(Int, 0.5 * length(Lambda))
    tmp1 = Lambda[pos1]
    pos2 = round(Int, 0.9 * length(Lambda))
    tmp2 = Lambda[pos2]
    EE1 = 0.0
    for j in 1:length(tmp1)
        EE1 += -log2(tmp1[j]^2) * tmp1[j]^2
    end
    EE2 = 0.0
    for j in 1:length(tmp2)
        EE2 += -log2(tmp2[j]^2) * tmp2[j]^2
    end
    EE1, EE2
end

function InnerProduct(Gamma1, Lambda1, Gamma2, Lambda2)
    # Take the inner product of two MPS chain
    N = length(Lambda1) - 1
    Dim = zeros(Int, N)
    Beta1 = zeros(Int, N+1)
    Beta2 = zeros(Int, N+1)
    for i in 1:N+1
        Beta1[i] = length(Lambda1[i])
        Beta2[i] = length(Lambda2[i])
    end
    for j in 1:N
        Dim[j] = size(Gamma1[j])[2]
    end

    TMP = zeros(Complex{Float64}, Beta2[2], Beta1[2])
    for B in 1:Beta1[2], B_ in Beta2[2], i in 1:Dim[1], A1 in 1:Beta1[1], A2 in 1:Beta2[1]
        TMP[B_, B] += (Lambda2[1][A2] * Gamma2[1][A2, i, B_])' * Lambda1[1][A1] * Gamma1[1][A1, i, B]
    end
    for n in 2:N-1
        tmp = zeros(Complex{Float64}, Beta2[n+1], Beta1[n+1])
        for B2 in 1:Beta1[n+1], B_2 in 1:Beta2[n+1], i in 1:Dim[n], B in 1:Beta1[n], B_ in 1:Beta2[n]
                tmp[B_2, B2] += TMP[B_, B] * (Lambda2[n][B_] * Gamma2[n][B_, i, B_2])' * Lambda1[n][B] * Gamma1[n][B, i, B2]
        end
        TMP = tmp
    end
    s = zero(Complex{Float64})
    for G_ in 1:Beta2[N+1], G in 1:Beta1[N+1], i in 1:Dim[N], B in 1:Beta1[N], B_ in 1:Beta2[N]
        s += TMP[B_, B] * (Lambda2[N][B_] * Gamma2[N][B_, i, G_] * Lambda2[N+1][G_])' * Lambda1[N][B] * Gamma1[N][B, i, G] * Lambda1[N+1][G]
    end
    s
end

function Fidelity(Gamma, Lambda, Gamma1, Lambda1, Gamma2, Lambda2, Gamma3, Lambda3)
    # calculate the fidelity of photons
    s = zeros(Complex{Float64}, 3)
    s[1] = InnerProduct(Gamma, Lambda, Gamma1, Lambda1)
    s[2] = InnerProduct(Gamma, Lambda, Gamma2, Lambda2)
    s[3] = InnerProduct(Gamma, Lambda, Gamma3, Lambda3)
    s
end

function MPS_Simulation(p::MPS_parameter, g0)
    # Perform the MPS simulation
    Gamma, Lambda = MPS_chain(p)
    U = TwoOperator(p)
	pAux = deepcopy(p)
    pAux.gdt = 0
    pAux.vacuum = 0
    Gamma2Aux, Lambda2Aux = MPS_chain(pAux) # |10>
    Gamma3Aux, Lambda3Aux = MPS_chain(pAux) # |01>
    Gamma4Aux, Lambda4Aux = MPS_chain(pAux) # |11>
    UAux = TwoOperator(pAux)

	filename = pwd() * "/data/" * string(g0) * ".txt"
    outfile = open(filename, "a")

    EE1 = Float64[]
    EE2 = Float64[]

	# prepare the initial state of the MPS chain
    Line = LineShape(p.width)
    x1 = round(Int, p.Tm * 0.1)
    x2 = round(Int, p.Tm * 0.25) + x1
    Gamma[x1+1:x1+p.width], Lambda[x1+2:x1+p.width] = OnePhoton(Line, p.vacuum, p.DIM_B)
    Gamma[x2+1:x2+p.width] = deepcopy(Gamma[x1+1:x1+p.width])
    Lambda[x2+2:x2+p.width] = deepcopy(Lambda[x1+2:x1+p.width])
    OneGamma, OneLambda = OnePhoton(Line, 0.0, p.DIM_B)
    Gamma2Aux[x1+1:x1+p.width] = deepcopy(OneGamma)
    Lambda2Aux[x1+2:x1+p.width] = deepcopy(OneLambda)
    Gamma3Aux[x2+1:x2+p.width] = deepcopy(OneGamma)
    Lambda3Aux[x2+2:x2+p.width] = deepcopy(OneLambda)
    Gamma4Aux[x1+1:x1+p.width] = deepcopy(OneGamma)
    Lambda4Aux[x1+2:x1+p.width] = deepcopy(OneLambda)
    Gamma4Aux[x2+1:x2+p.width] = deepcopy(OneGamma)
    Lambda4Aux[x2+2:x2+p.width] = deepcopy(OneLambda)

    push!(EE1, 0.0)
    push!(EE2, 0.0)

    for i in 1:p.length
        Gamma[i:i+1], Lambda[i+1] = NQubitGate(U, Gamma[i:i+1], Lambda[i:i+2])
        Gamma[i:i+1], Lambda[i+1] = SwapGate(Gamma[i:i+1], Lambda[i:i+2])
        Gamma2Aux[i:i+1], Lambda2Aux[i+1] = NQubitGate(UAux, Gamma2Aux[i:i+1], Lambda2Aux[i:i+2])
        Gamma2Aux[i:i+1], Lambda2Aux[i+1] = SwapGate(Gamma2Aux[i:i+1], Lambda2Aux[i:i+2])
        Gamma3Aux[i:i+1], Lambda3Aux[i+1] = NQubitGate(UAux, Gamma3Aux[i:i+1], Lambda3Aux[i:i+2])
        Gamma3Aux[i:i+1], Lambda3Aux[i+1] = SwapGate(Gamma3Aux[i:i+1], Lambda3Aux[i:i+2])
        Gamma4Aux[i:i+1], Lambda4Aux[i+1] = NQubitGate(UAux, Gamma4Aux[i:i+1], Lambda4Aux[i:i+2])
        Gamma4Aux[i:i+1], Lambda4Aux[i+1] = SwapGate(Gamma4Aux[i:i+1], Lambda4Aux[i:i+2])
    end
    e1, e2 = EECal(Lambda)
    push!(EE1, e1)
    push!(EE2, e2)
    for j in 2:p.Bounce
        for i in p.length:-1:1
            Gamma[i:i+1], Lambda[i+1] = SwapGate(Gamma[i:i+1], Lambda[i:i+2])
			Gamma2Aux[i:i+1], Lambda2Aux[i+1] = SwapGate(Gamma2Aux[i:i+1], Lambda2Aux[i:i+2])
            Gamma3Aux[i:i+1], Lambda3Aux[i+1] = SwapGate(Gamma3Aux[i:i+1], Lambda3Aux[i:i+2])
            Gamma4Aux[i:i+1], Lambda4Aux[i+1] = SwapGate(Gamma4Aux[i:i+1], Lambda4Aux[i:i+2])
        end
        for i in 1:p.length
            Gamma[i:i+1], Lambda[i+1] = NQubitGate(U, Gamma[i:i+1], Lambda[i:i+2])
            Gamma[i:i+1], Lambda[i+1] = SwapGate(Gamma[i:i+1], Lambda[i:i+2])
			Gamma2Aux[i:i+1], Lambda2Aux[i+1] = NQubitGate(UAux, Gamma2Aux[i:i+1], Lambda2Aux[i:i+2])
            Gamma2Aux[i:i+1], Lambda2Aux[i+1] = SwapGate(Gamma2Aux[i:i+1], Lambda2Aux[i:i+2])
            Gamma3Aux[i:i+1], Lambda3Aux[i+1] = NQubitGate(UAux, Gamma3Aux[i:i+1], Lambda3Aux[i:i+2])
            Gamma3Aux[i:i+1], Lambda3Aux[i+1] = SwapGate(Gamma3Aux[i:i+1], Lambda3Aux[i:i+2])
            Gamma4Aux[i:i+1], Lambda4Aux[i+1] = NQubitGate(UAux, Gamma4Aux[i:i+1], Lambda4Aux[i:i+2])
            Gamma4Aux[i:i+1], Lambda4Aux[i+1] = SwapGate(Gamma4Aux[i:i+1], Lambda4Aux[i:i+2])
        end
		if mod(j, 2) == 0
            s = Fidelity(Gamma, Lambda, Gamma2Aux, Lambda2Aux, Gamma3Aux, Lambda3Aux, Gamma4Aux, Lambda4Aux)
            writedlm(outfile, s')
        end
        e1, e2 = EECal(Lambda)
        push!(EE1, e1)
        push!(EE2, e2)
    end
    writedlm(outfile, EE1')
    writedlm(outfile, EE2')
    close(outfile)
end

function g0Sweep(g0)
    println("$g0")
    GM = 200
    OM = 0.03
    NB = 20
    p = MPS_parameter(gamma=GM, omega_m=OM, g_0=g0, t_total=0.5, N=2, dim_a=2, dim_s=12, dim_b=2, vacuum=0.5, widthfactor=1000, NumBounce=NB)
    MPS_Simulation(p, g0)
end

end
