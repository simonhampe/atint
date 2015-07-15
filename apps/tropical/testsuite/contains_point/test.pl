my $c1 = uniform_linear_space<Max>(3,2);
my $c2 = cartesian_product<Min>( (uniform_linear_space<Min>(2,1)), (projective_torus<Min>(1)));
my $c3 = cross_variety<Max>(2,1,0);
my $c4 = empty_cycle<Min>(3);

compare_values("1",1,contains_point($c1, [1,0,1,-1,1]))
	and
compare_values("2",1,contains_point($c2, [1,0,0,0,-1]))
	and
compare_values("3",0,contains_point($c3, [1,0,1,0]))
	and
compare_values("4",0,contains_point($c4, [1,0,0,0,0]));

