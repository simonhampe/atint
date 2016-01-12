my $fct1 = new RationalFunction<Min>(NUMERATOR=>toTropicalPolynomial("min(2x,y+z)"), DENOMINATOR=>toTropicalPolynomial("min(2x,2y,2z)"));
my $cycle1 = uniform_linear_space<Min>(2,1);

my $dom2 = uniform_linear_space<Max>(3,2);
my $fct2 = new RationalFunction<Max>(DOMAIN=>$dom2,VERTEX_VALUES=>[1,2,3,4,5],LINEALITY_VALUES=>[]);
my $cycle2 = new Cycle<Max>(VERTICES=>[[1,0,0,0,0],[0,0,1,3,3],[0,0,1,-1/3,1],[0,0,-1,0,-6/5],[0,0,0,-1,-3/5]],MAXIMAL_POLYTOPES=>[[0,1],[0,2],[0,3],[0,4]],WEIGHTS=>[2,1,1,1]);

compare_object("1", $fct1->restrict($cycle1));

compare_object("2", $fct2->restrict($cycle2));
