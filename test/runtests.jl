using Test, HANK_MNS

@test hello("Julia") == "Hello, Julia"
@test domath(2.0) ≈ 7.0
@test domath2(2.0) ≈ 8.0
