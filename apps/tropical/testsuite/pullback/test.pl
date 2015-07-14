my $m1 = morphism_from_affine<Max>([[0,1],[1,0]],[1,1]);
my $r1 = new RationalFunction<Max>( NUMERATOR=>toTropicalPolynomial("max(2x,2y,2z)"), DENOMINATOR=>toTropicalPolynomial("max(2y,4x-2z)"));

my $d2 = new Cycle<Min>(VERTICES=>[[1,0,0,0],[1,0,1,0],[0,0,0,1],[0,0,-1,-1],[0,0,1,-1]],MAXIMAL_POLYTOPES=>[[0,1],[0,2],[0,3],[1,2],[1,4]], WEIGHTS=>[1,1,1,1,1]);
my $m2 = new Morphism<Min>(DOMAIN=>$d2, VERTEX_VALUES=>[[0,0,0,0],[0,1,1,0],[0,-1,0,0],[0,0,0,-1],[0,0,-1,0],[0,1,1,1]]);
my $r2 = new RationalFunction<Min>(NUMERATOR=>toTropicalPolynomial("min(x,z,y,w)"), DENOMINATOR=>toTropicalPolynomial("min(x)",qw(x y z w)));

my $m3 = morphism_from_affine<Max>([[-1,0],[0,-1]],[0,0]);
my $d3 = uniform_linear_space<Max>(2,1);
my $r3 = new RationalFunction<Max>(DOMAIN=>$d3,VERTEX_VALUES=>[0,1,2,3]);

my $d4 = new Cycle<Min>(VERTICES=>[[1,0,0,0,0],[1,0,1,1,0],[0,0,-1,0,0],[0,0,0,0,-1],[0,0,0,-1,0],[0,0,1,1,1]],MAXIMAL_POLYTOPES=>[[0,1],[0,2],[0,4],[1,3],[1,5]],WEIGHTS=>[1,1,1,1,1]);
my $r4 = new RationalFunction<Min>(DOMAIN=>$d4, VERTEX_VALUES=>[0,1,2,3,2,-1]);

my $d5 = point_collection<Max>([[0,1,2]],[1]);
my $m5 = new Morphism<Max>(DOMAIN=>$d5,VERTEX_VALUES=>[[0,3,5]]);
my $r5 = new RationalFunction<Max>(NUMERATOR=>toTropicalPolynomial("max(a,b,c)"), DENOMINATOR=>toTropicalPolynomial("max(a)",qw(a b c)));

my $m6 = morphism_from_affine<Min>([[1,0],[0,1]],[0,0]);
my $d6 = point_collection<Min>([[0,0,0]],[1]);
my $r6 = new RationalFunction<Min>(DOMAIN=>$d6,VERTEX_VALUES=>[3]);

compare_object("1", pullback($m1,$r1))
	and
compare_object("2", pullback($m2,$r2))
	and
compare_object("3", pullback($m3,$r3))
	and
compare_object("4", pullback($m2,$r4))
	and
compare_object("5", pullback($m5,$r5))
	and
compare_object("6", pullback($m6, $r6));

