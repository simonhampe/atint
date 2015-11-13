my $v1=new Vector<Int>([1,1,1,1,-4]);
compare_object("1", (hurwitz_subdivision<Max>(1, $v1, Verbose=>0)));

my $v2=new Vector<Rational>([1]);
compare_object("2", (hurwitz_subdivision<Min>(1, $v1, $v2, Verbose=>0)));
