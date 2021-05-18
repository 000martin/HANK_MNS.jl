using Test, HANK_MNS

@testset "Compare computation of steady state and transition with MNS baseline" begin

    p0 = set_parameters()

 @test typeof(p0) == params

    b,y,SS = get_steady_state(p0,0.6,[0.95,0.99])

 #comparison with MNS results
 @test abs(b - 0.986) < 1e-4   #compare β
 @test abs(y -  1.0368) < 1e-4 #compare output

 p1 = set_parameters(β=b)

 tp = get_transition_full(20,200,p1,SS)

 #compare initial output response - ca. 0.1 percentage points for output
 @test abs((tr.Y[1]/SS.Y) - 1.001) <1e-4

 #compare initial inflation response - ca. 0.3 percentage points for inflation
 @test abs(tr.pΠ[1]-1.003) <1e-4

end;
