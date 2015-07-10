my $fct1 = new RationalFunction<Min>(NUMERATOR=>toTropicalPolynomial("min(3x)",qw(x y z)), DENOMINATOR=>toTropicalPolynomial("min(3z)", qw(x y z)));

my $fct2 = new RationalFunction<Min>(NUMERATOR=>toTropicalPolynomial("min(2x+2y)", qw(x y z)), DENOMINATOR=>toTropicalPolynomial("min(3z)", qw(x y z)));

my $g1 = rational_fct_from_affine_numerator(toTropicalPolynomial("max(0,x,y)"));
my $g2 = new RationalFunction<Max>(DOMAIN=>$g1->DOMAIN, VERTEX_VALUES=>[0,1,1,1],LINEALITY_VALUES=>[]);

compare_object("1", 2*$fct1 + 3*$fct2)
	and
compare_object("2", $g1 + $g2);


