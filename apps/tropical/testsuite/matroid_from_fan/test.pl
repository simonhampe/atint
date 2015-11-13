my $a = point_collection<Max>([[0,0,0,0]],[1]);
my $b = affine_linear_space<Min>([[0,1,1]]);
my $c = cartesian_product( (uniform_linear_space<Max>(2,1)), (projective_torus<Max>(1)));
my $d = m0n<Min>(5);

compare_object("1", matroid_from_fan($a));

compare_object("2", matroid_from_fan($b));

compare_object("3", matroid_from_fan($c));

compare_object("4", matroid_from_fan($d));
