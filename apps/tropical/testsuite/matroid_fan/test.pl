compare_object( '1', matroid_fan<Max>( matroid::matroid_from_graph(graph::complete(4))));

compare_object( '2', matroid_fan<Min>([[1,1,0],[0,0,1]]));

compare_object( '3', matroid_fan<Max>(new matroid::Matroid(N_ELEMENTS=>2, BASES=>[[0]])));
