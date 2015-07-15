my $c1 = uniform_linear_space<Min>(3,2);
my $c2 = new Cycle<Min>(VERTICES=>[[1,0,0,0],[1,0,1,0],[0,0,0,1],[0,0,-1,-1],[0,0,1,-1]],MAXIMAL_POLYTOPES=>[[0,1],[0,2],[0,3],[1,2],[1,4]], WEIGHTS=>[1,1,1,1,1]);
my $c3 = affine_linear_space<Max>([[0,1,0,1],[0,0,1,1]]);

compare_object("1", skeleton_complex($c1,1))
	and
compare_object("2", skeleton_complex($c2,0))
	and
compare_object("3", skeleton_complex($c2,0,1))
	and
compare_object("4", skeleton_complex($c3,1));
