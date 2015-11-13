my $o1 = orthant_subdivision<Max>([1,0,0,0]);
my $o2 = empty_cycle<Min>(5);
my $o3 = matroid_fan_from_flats<Max>(matroid::uniform_matroid(3,4));

compare_object("1", coarsen($o1));

compare_object("2", coarsen($o3));
