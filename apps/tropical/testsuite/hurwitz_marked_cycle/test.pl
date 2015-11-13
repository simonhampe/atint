my $v1=new Vector<Int>(1,1,1,1,-4);
compare_object("1", (hurwitz_marked_cycle<Max>(1, $v1)));

my $v2=new Vector<Rational>([1]);
compare_object("2", (hurwitz_marked_cycle<Min>(1, $v1, $v2)));
