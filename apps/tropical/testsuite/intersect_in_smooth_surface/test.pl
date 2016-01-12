#Create test objects

my $bisector = affine_linear_space<Max>([[0,1,1,0]]);
my $stdplane = uniform_linear_space<Max>(3,2);

my $lineinspace = uniform_linear_space<Min>(3,1);
my $linetimesline = cartesian_product( (uniform_linear_space<Min>(2,1)), projective_torus<Min>(1));

my $actualline = affine_linear_space<Min>([[0,0,0,1]]);

my $plane = projective_torus<Max>(2);
my $curveinplane = new Cycle<Max>(VERTICES=>[[1,0,0,0],[1,0,1,0],[0,0,0,1],[0,0,-1,-1],[0,0,1,-1]],MAXIMAL_POLYTOPES=>[[0,1],[0,2],[0,3],[1,2],[1,4]],WEIGHTS=>[1,1,1,1,1]);
my $stdline = uniform_linear_space<Max>(2,1);
my $lineinplane = affine_linear_space<Max>([[0,1,1,]]);

compare_object("1", intersect_in_smooth_surface( $stdplane, $bisector,$bisector));

compare_object("2", intersect_in_smooth_surface( $linetimesline, $lineinspace, $lineinspace));

compare_object("3", intersect_in_smooth_surface( $linetimesline, $actualline, $actualline));

compare_object("4", intersect_in_smooth_surface( $linetimesline, $lineinspace, $actualline));

compare_object("5", intersect_in_smooth_surface( $plane, $curveinplane, $stdline));

compare_object("6", intersect_in_smooth_surface( $plane, $lineinplane, $lineinplane));
