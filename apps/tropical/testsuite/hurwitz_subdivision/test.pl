compare_object("1",hurwitz_subdivision<Max>(1, (new Vector<Int>([1,1,1,1,-4])), Verbose=>0) )
	and
compare_object("2",hurwitz_subdivision<Min>(1, (new Vector<Int>(1,1,1,1,-4)), (new Vector<Rational>([1])), Verbose=>0) );
