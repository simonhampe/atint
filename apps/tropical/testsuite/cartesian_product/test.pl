my $l1 = uniform_linear_space<Max>(2,1);
my $locall1 = local_vertex($l1,0);
my $line = new Cycle<Max>(VERTICES=>[[1,0,0],[0,1,0],[0,0,1]],MAXIMAL_POLYTOPES=>[[0,1],[0,2]],WEIGHTS=>[2,2]);
my $minline = uniform_linear_space<Min>(2,1);
my $minprod = cartesian_product($minline, projective_torus<Min>(1));
my $minotherprod = cartesian_product( (projective_torus<Min>(1)), $minline);
my $points = point_collection<Max>([[0,1,2],[0,2,3]],[1,1]);

compare_object("1", cartesian_product($l1,$l1));

compare_object("2", cartesian_product($locall1, $line));

compare_object("3", cartesian_product($minprod, $minotherprod ));

compare_object("4", cartesian_product($points, empty_cycle<Max>(1)));
