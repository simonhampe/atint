my $m1 = matroid::uniform_matroid(3,4);
my $m2 = new matroid::Matroid(N_ELEMENTS=>4, BASES=>[[0,2,3],[1,2,3]]);

compare_object("1", matroid_fan_from_flats<Max>($m1))
	and
compare_object("2", matroid_fan_from_flats<Min>($m2));
