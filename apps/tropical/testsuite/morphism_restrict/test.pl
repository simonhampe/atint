my $m1 = new Morphism<Max>(MATRIX=>[[1,0,1],[0,2,0],[3,0,-1]],TRANSLATE=>[1,1,-4]);
my $m2 = new Morphism<Min>(TRANSLATE=>[1,-10]);
my $c2 = new Cycle<Min>(PROJECTIVE_VERTICES=>[[1,0,0],[1,1,0],[1,2,3]],MAXIMAL_POLYTOPES=>[[0],[1],[2]],WEIGHTS=>[1,1,1]);
my $m3 = new Morphism<Max>(DOMAIN=>(affine_linear_space<Max>([[0,1,1]])),MATRIX=>[[1,0,0],[0,3,-2],[-1,0,2]]);

compare_object( '1', $m1->restrict( uniform_linear_space<Max>(2,1) ));

compare_object( '2', $m2->restrict($c2));

compare_object( '3', $m3->restrict(uniform_linear_space<Max>(2,1))); 
