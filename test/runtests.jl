using Test, HANK_MNS


@testset "Testing some basic functionalities" begin

   p = set_parameters()
   
   #sanity check grid
   @test maximum(p.k_grid) == p.a_max
   @test minimum(p.k_grid) == p.a_max

   #sanity test income grid
   @test minimum(p.z) > 0

   #test reshape_c
   @test reshape_c(ones(p.nb,p.nz),p) == ones(p.nb*p.nz)
   
   #check income state transition process
   @test sum(p.Πz,dims=2) == ones(p.nz,1)
   @test sum(p.Γ) ≈ 1.0

   #if next period policy function is always one, unconstrained agents in highest income states should also
   #1 as policy function in current period if R = 1/β.
   cs = HANK_MNS.EGM(ones(p.nb,p.nz),p.β,ones(2)./p.β,ones(2),zeros(2),0.1*ones(2),p)
   @test cs[:,3] == ones(p.nb)

   #test transition matrix - rows should sum up to one
   Pi = forwardmat(cs,1/p.β,1.0,0.0,0.1,p)
   @test sum(Pi,dims = 2) ≈ ones(p.nk*p.nz,1)
   #Pi should be square matrix with size nk*nz
   @test size(Pi) == (p.nk*p.nz,p.nk*p.nz) 

   #invariant distribution induced by Pi should sum up to 1
   @test sum(inv_dist(Pi)) ≈ 1.0

end;


@testset "Compare computation of steady state and transition with MNS baseline" begin

    p0 = set_parameters()

 @test typeof(p0) == params

    b,y,SS = get_steady_state(p0,0.6,[0.95,0.99])

 #comparison with MNS results
 @test abs(b - 0.986) < 1e-4   #compare β
 @test abs(y -  1.0368) < 1e-4 #compare steady state output level

 p1 = set_parameters(β=b)

 tp = get_transition_full(20,200,p1,SS)

 #compare initial output response - ca. 0.1 percentage points for output
 @test abs((tp.Y[1]/SS.Y) - 1.001) <1e-4

 #compare initial inflation response - ca. 0.3 percentage points for inflation
 @test abs(tp.pΠ[1]-1.003) <1e-4

end;

