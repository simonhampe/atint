my $f1 = empty_cycle<Min>(3);
my $f2 = new Cycle<Max>(VERTICES=>[[1,0,0,0],[1,0,1,0],[0,0,0,1],[0,0,-1,-1],[0,0,1,-1]],MAXIMAL_POLYTOPES=>[[0,1],[0,2],[0,3],[1,2],[1,4]],WEIGHTS=>[1,1,1,1,1]);
my $f3 = cross_variety<Min>(2,1);

my @dec1 = fan_decomposition($f1);
my @dec2 = fan_decomposition($f2);
my @dec3 = fan_decomposition($f3);

compare_values("1n",0,scalar(@dec1));

compare_values("2n",2,scalar(@dec2));

compare_values("3n",4,scalar(@dec3));

compare_object("2_1",$dec2[0]);

compare_object("2_2",$dec2[1]);

compare_object("3_1",$dec3[0]);

compare_object("3_2",$dec3[1]);

compare_object("3_3",$dec3[2]);

compare_object("3_4",$dec3[3]);
