my $cycle1 = point_collection<Min>([[0,1,2],[0,2,3],[0,-1,-1]],[1,1,1]);
my $cycle2 = uniform_linear_space<Max>(3,2);
my $cycle3 = affine_linear_space<Min>([[0,1,2,0],[0,0,2,3]]);
my $morph3 = new Morphism<Min>(MATRIX=>[[1,0,0,0],[0,0,1,0],[0,1,0,0],[0,0,0,1]]);


compare_object("1", affine_transform($cycle1, [[1,0,0],[0,0,1],[0,1,0]],[2,3,4]))
	and
compare_object("2", shift_cycle($cycle2, [1,2,3,4]))
	and
compare_object("3", affine_transform($cycle3,$morph3))
	and
compare_object("4", affine_transform( (empty_cycle<Max>(3)), [[1,0,0,0],[0,1,0,0],[0,0,1,0]]));

