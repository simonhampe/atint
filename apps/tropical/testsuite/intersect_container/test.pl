my $a1 = projective_torus<Max>(2);
my $c1 = rational_fct_from_affine_numerator(toTropicalPolynomial("max(0,x,y,x+y, 2x-1, 2y-1)"))->DOMAIN;

my $a2 = new Cycle<Min>(VERTICES=>[[1,0,0,0],[1,0,1,0],[0,0,0,1],[0,0,-1,-1],[0,0,1,-1]],MAXIMAL_POLYTOPES=>[[0,1],[0,2],[0,3],[1,2],[1,4]], WEIGHTS=>[1,1,1,1,1]); 
my $c2 = orthant_subdivision<Min>([1,0,1/2,1]);

my $a3 = point_collection<Max>([[0,0,0]],[3]);
my $c3 = cross_variety<Max>(2,1,0);

my $a4 = affine_linear_space<Min>([[0,1,0,1],[0,0,1,1]]);
my $c4 = uniform_linear_space<Min>(3,3);

my $a5 = empty_cycle<Max>(5);
my $c5 = projective_torus<Max>(5);

compare_object("1", intersect_container($a1,$c1));

compare_object("2", intersect_container($a2,$c2));

compare_object("3", intersect_container($a3,$c3));

compare_object("4", intersect_container($a4,$c4));

compare_object("5", intersect_container($a5,$c5));
