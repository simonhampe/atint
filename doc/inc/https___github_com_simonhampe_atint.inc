<div xmlns:html="http://www.w3.org/1999/xhtml" xmlns:pm="http://www.polymake.org/ns/docs#3" id="content">
  <h1 class="extension">extension: atint</h1>
  <p>  a-tint is an extension for tropical intersection theory.<br xmlns="http://www.w3.org/1999/xhtml"/>  (C) 2011 - 2015 Simon Hampe, simon.hampe@googlemail.com<br xmlns="http://www.w3.org/1999/xhtml"/>  See also: <br xmlns="http://www.w3.org/1999/xhtml"/>  https://github.com/simonhampe/atint</p>
  <h1>application: tropical</h1>
  <div class="level2">
    <h2>
      <b>Objects</b>
    </h2>
    <ul class="unfoldable">
      <li>
        <div class="li">
          <div class="icon" id="icon:tropical__Cycle" onclick="swap_content('span:tropical__Cycle'); return false;">
		 	
	</div>
          <div class="unfoldable">
            <div class="object">
              <h3>
                <a href="#" class="javalink" id="tropical__Cycle" onclick="swap_content('span:tropical__Cycle'); return false;">
                  <span class="object">Cycle</span>
                </a>
              </h3>
            </div>
            <div style="display: inline" id="span:tropical__Cycle" class="foldit">
              <div class="level3"/>
              <h4>Properties of Cycle</h4>
              <ul class="unfoldable">
                <li>
                  <div class="li">
                    <div class="icon" id="icon:tropical__Cycle__BALANCED_FACES" onclick="swap_content('span:tropical__Cycle__BALANCED_FACES'); return false;">
		 	
	</div>
                    <div class="unfoldable"><a href="#" class="javalink" id="tropical__Cycle__BALANCED_FACES" onclick="swap_content('span:tropical__Cycle__BALANCED_FACES'); return false;"><b>BALANCED_FACES</b>:
		</a>common::Vector&lt;Bool&gt;<br/><div style="display: inline" id="span:tropical__Cycle__BALANCED_FACES" class="foldit"><div class="descr"><p> A vector whose entries correspond to the rows of CODIMENSION_ONE_POLYTOPES. 
 The i-th entry is true, if and only if the complex is balanced at that face</p></div></div></div>
                  </div>
                </li>
                <li>
                  <div class="li">
                    <div class="icon" id="icon:tropical__Cycle__DEGREE" onclick="swap_content('span:tropical__Cycle__DEGREE'); return false;">
		 	
	</div>
                    <div class="unfoldable"><a href="#" class="javalink" id="tropical__Cycle__DEGREE" onclick="swap_content('span:tropical__Cycle__DEGREE'); return false;"><b>DEGREE</b>:
		</a>common::Integer<br/><div style="display: inline" id="span:tropical__Cycle__DEGREE" class="foldit"><div class="descr"><p> The degree of the tropical variety, i.e. the weight of the intersection product with a
 uniform tropical linear space of complementary dimension.</p></div></div></div>
                  </div>
                </li>
                <li>
                  <div class="li">
                    <div class="icon" id="icon:tropical__Cycle__IS_BALANCED" onclick="swap_content('span:tropical__Cycle__IS_BALANCED'); return false;">
		 	
	</div>
                    <div class="unfoldable"><a href="#" class="javalink" id="tropical__Cycle__IS_BALANCED" onclick="swap_content('span:tropical__Cycle__IS_BALANCED'); return false;"><b>IS_BALANCED</b>:
		</a>common::Bool<br/><div style="display: inline" id="span:tropical__Cycle__IS_BALANCED" class="foldit"><div class="descr"><p> Whether the cycle is balanced. As many functions in a-tint can 
 deal with non-balanced complexes, we include this as a property.</p></div></div></div>
                  </div>
                </li>
                <li>
                  <div class="li">
                    <div class="icon" id="icon:tropical__Cycle__IS_IRREDUCIBLE" onclick="swap_content('span:tropical__Cycle__IS_IRREDUCIBLE'); return false;">
		 	
	</div>
                    <div class="unfoldable"><a href="#" class="javalink" id="tropical__Cycle__IS_IRREDUCIBLE" onclick="swap_content('span:tropical__Cycle__IS_IRREDUCIBLE'); return false;"><b>IS_IRREDUCIBLE</b>:
		</a>common::Bool<br/><div style="display: inline" id="span:tropical__Cycle__IS_IRREDUCIBLE" class="foldit"><div class="descr"><p> Whether this complex is irreducible.</p></div></div></div>
                  </div>
                </li>
                <li>
                  <div class="li">
                    <div class="icon" id="icon:tropical__Cycle__LATTICE_BASES" onclick="swap_content('span:tropical__Cycle__LATTICE_BASES'); return false;">
		 	
	</div>
                    <div class="unfoldable"><a href="#" class="javalink" id="tropical__Cycle__LATTICE_BASES" onclick="swap_content('span:tropical__Cycle__LATTICE_BASES'); return false;"><b>LATTICE_BASES</b>:
		</a>common::IncidenceMatrix&lt;NonSymmetric&gt;<br/><div style="display: inline" id="span:tropical__Cycle__LATTICE_BASES" class="foldit"><div class="descr"><p> This incidence matrix gives a lattice basis for each maximal polytope. More precisely it 
 gives a lattice basis whose span contains the lattice of the maximal polytope. Row i 
 corresponds to cone i and gives lattice generator indices referring to <a href="&#10;&#9;&#9;&#9;#tropical__Cycle__LATTICE_GENERATORS" onclick="unfold('span:tropical__Cycle__LATTICE_GENERATORS');">LATTICE_GENERATORS</a>.
 If this property is computed via rules, it does indeed give a lattice basis for the
 cone lattice, but when it is computed during an operation like refinement or divisor
 it will in general be larger.
 If this property exists, lattice normals might be computed faster.</p></div></div></div>
                  </div>
                </li>
                <li>
                  <div class="li">
                    <div class="icon" id="icon:tropical__Cycle__LATTICE_GENERATORS" onclick="swap_content('span:tropical__Cycle__LATTICE_GENERATORS'); return false;">
		 	
	</div>
                    <div class="unfoldable"><a href="#" class="javalink" id="tropical__Cycle__LATTICE_GENERATORS" onclick="swap_content('span:tropical__Cycle__LATTICE_GENERATORS'); return false;"><b>LATTICE_GENERATORS</b>:
		</a>common::Matrix&lt;Integer, NonSymmetric&gt;<br/><div style="display: inline" id="span:tropical__Cycle__LATTICE_GENERATORS" class="foldit"><div class="descr"><p> This is an irredundant list of all lattice generators of all maximal polyhedra.
 If this property exists, lattice normals might be computed faster</p></div></div></div>
                  </div>
                </li>
                <li>
                  <div class="li">
                    <div class="icon" id="icon:tropical__Cycle__LATTICE_NORMAL_FCT_VECTOR" onclick="swap_content('span:tropical__Cycle__LATTICE_NORMAL_FCT_VECTOR'); return false;">
		 	
	</div>
                    <div class="unfoldable"><a href="#" class="javalink" id="tropical__Cycle__LATTICE_NORMAL_FCT_VECTOR" onclick="swap_content('span:tropical__Cycle__LATTICE_NORMAL_FCT_VECTOR'); return false;"><b>LATTICE_NORMAL_FCT_VECTOR</b>:
		</a>common::Map&lt;Pair&lt;Int, Int&gt;, Vector&lt;Rational&gt;&gt;<br/><div style="display: inline" id="span:tropical__Cycle__LATTICE_NORMAL_FCT_VECTOR" class="foldit"><div class="descr"><p> For each lattice normal vector, this gives a vector of length (number of rays) + (lineality dim.), such that
 if a rational function is given by values on the rays and lin space generators, 
 the value of the corresponding normal <a href="tropical.html#Cycle__LATTICE_NORMALS" onclick="unfold('span:Cycle__LATTICE_NORMALS');">LATTICE_NORMALS</a>-&gt;{i}-&gt;{j} can be computed by 
 multiplying the function value vector with the vector <a href="&#10;&#9;&#9;&#9;#tropical__Cycle__LATTICE_NORMAL_FCT_VECTOR" onclick="unfold('span:tropical__Cycle__LATTICE_NORMAL_FCT_VECTOR');">LATTICE_NORMAL_FCT_VECTOR</a>-&gt;{i}-&gt;{j}. 
 This is done in the following way:  
 We use the generating system (and indices refer to SEPARATED_VERTICES) &lt;(r_i-r_0)_i&gt;0, s_j, l_k&gt;, 
 where r_0 is the ray of the maximal cone with the lowest index in <a href="&#10;&#9;&#9;&#9;#tropical__Cycle__SEPARATED_VERTICES" onclick="unfold('span:tropical__Cycle__SEPARATED_VERTICES');">SEPARATED_VERTICES</a>,
 such that it fulfills x0 = 1, r_i are the remaining rays with x0 = 1, ordered according to their index 
 in <a href="&#10;&#9;&#9;&#9;#tropical__Cycle__SEPARATED_VERTICES" onclick="unfold('span:tropical__Cycle__SEPARATED_VERTICES');">SEPARATED_VERTICES</a>, s_j are the rays of the cone with x0 = 0 and l_k are the lineality space 
 generators. We will then store the coefficients a_i of (r_i - r_0) at the index of r_i, then - sum(a_i) at 
 the index of r_0 and the remaining coefficients at the appropriate places.
 In particular, the value of a lattice normal under a rational function can be computed
 simply by taking the scalar product of RAY_VALUES | LIN_VALUES with this FCT_VECTOR</p></div></div></div>
                  </div>
                </li>
                <li>
                  <div class="li">
                    <div class="icon" id="icon:tropical__Cycle__LATTICE_NORMAL_SUM" onclick="swap_content('span:tropical__Cycle__LATTICE_NORMAL_SUM'); return false;">
		 	
	</div>
                    <div class="unfoldable"><a href="#" class="javalink" id="tropical__Cycle__LATTICE_NORMAL_SUM" onclick="swap_content('span:tropical__Cycle__LATTICE_NORMAL_SUM'); return false;"><b>LATTICE_NORMAL_SUM</b>:
		</a>common::Matrix&lt;Rational, NonSymmetric&gt;<br/><div style="display: inline" id="span:tropical__Cycle__LATTICE_NORMAL_SUM" class="foldit"><div class="descr"><p> Rows of this matrix correspond to <a href="tropical.html#Cycle__CODIMENSION_ONE_POLYTOPES" onclick="unfold('span:Cycle__CODIMENSION_ONE_POLYTOPES');">CODIMENSION_ONE_POLYTOPES</a>, and each row contains the weighted sum:
 sum_{cone &gt; codim-1-face}( weight(cone) * <a href="tropical.html#Cycle__LATTICE_NORMALS" onclick="unfold('span:Cycle__LATTICE_NORMALS');">LATTICE_NORMALS</a>-&gt;{codim-1-face}-&gt;{cone})</p></div></div></div>
                  </div>
                </li>
                <li>
                  <div class="li">
                    <div class="icon" id="icon:tropical__Cycle__LATTICE_NORMAL_SUM_FCT_VECTOR" onclick="swap_content('span:tropical__Cycle__LATTICE_NORMAL_SUM_FCT_VECTOR'); return false;">
		 	
	</div>
                    <div class="unfoldable"><a href="#" class="javalink" id="tropical__Cycle__LATTICE_NORMAL_SUM_FCT_VECTOR" onclick="swap_content('span:tropical__Cycle__LATTICE_NORMAL_SUM_FCT_VECTOR'); return false;"><b>LATTICE_NORMAL_SUM_FCT_VECTOR</b>:
		</a>common::Matrix&lt;Rational, NonSymmetric&gt;<br/><div style="display: inline" id="span:tropical__Cycle__LATTICE_NORMAL_SUM_FCT_VECTOR" class="foldit"><div class="descr"><p> Rows of this matrix correspond to <a href="&#10;&#9;&#9;&#9;#tropical__Cycle__SEPARATED_CODIMENSION_ONE_POLYTOPES" onclick="unfold('span:tropical__Cycle__SEPARATED_CODIMENSION_ONE_POLYTOPES');">SEPARATED_CODIMENSION_ONE_POLYTOPES</a> and each row contains a function vector 
 for the corresponding row of <a href="&#10;&#9;&#9;&#9;#tropical__Cycle__LATTICE_NORMAL_SUM" onclick="unfold('span:tropical__Cycle__LATTICE_NORMAL_SUM');">LATTICE_NORMAL_SUM</a>. This function vector is computed 
 in the same way as described under <a href="&#10;&#9;&#9;&#9;#tropical__Cycle__LATTICE_NORMAL_FCT_VECTOR" onclick="unfold('span:tropical__Cycle__LATTICE_NORMAL_FCT_VECTOR');">LATTICE_NORMAL_FCT_VECTOR</a>.
 Note that for any codim-1-faces at which the complex is not balanced, the corresponding row is a zero row.
 If a face is balanced can be checked under <a href="&#10;&#9;&#9;&#9;#tropical__Cycle__BALANCED_FACES" onclick="unfold('span:tropical__Cycle__BALANCED_FACES');">BALANCED_FACES</a>.</p></div></div></div>
                  </div>
                </li>
                <li>
                  <div class="li">
                    <div class="icon" id="icon:tropical__Cycle__LOCAL_RESTRICTION" onclick="swap_content('span:tropical__Cycle__LOCAL_RESTRICTION'); return false;">
		 	
	</div>
                    <div class="unfoldable"><a href="#" class="javalink" id="tropical__Cycle__LOCAL_RESTRICTION" onclick="swap_content('span:tropical__Cycle__LOCAL_RESTRICTION'); return false;"><b>LOCAL_RESTRICTION</b>:
		</a>common::IncidenceMatrix&lt;NonSymmetric&gt;<br/><div style="display: inline" id="span:tropical__Cycle__LOCAL_RESTRICTION" class="foldit"><div class="descr"><p> This contains a list of sets of ray indices (referring to <a href="tropical.html#" onclick="unfold('span:');">VERTICES</a>). All of these sets should 
 describe polyhedra of the polyhedral complex (though not necessarily maximal ones).
 A polyhedron is now called compatible with this property, if it contains one of these polyhedra
 If this list is not empty, all computations will be done only on (or around) compatible cones. 
 The documentation of each property will explain in what way this restriction is enforced. 
 If this list is empty or not defined, there is no restriction.
 Careful: The implementation assumes that ALL maximal cones are compatible. If in doubt,
 you can create a complex with a local restriction from a given complex by using one of
 the "local_..." creation methods
 This list is assumed to be irredundant, i.e. there are no doubles (though this should not
 break anything, it is simply less efficient). It is, however, possible that one element
 is a subset of another.</p></div></div></div>
                  </div>
                </li>
                <li>
                  <div class="li">
                    <div class="icon" id="icon:tropical__Cycle__SEPARATED_CODIMENSION_ONE_POLYTOPES" onclick="swap_content('span:tropical__Cycle__SEPARATED_CODIMENSION_ONE_POLYTOPES'); return false;">
		 	
	</div>
                    <div class="unfoldable"><a href="#" class="javalink" id="tropical__Cycle__SEPARATED_CODIMENSION_ONE_POLYTOPES" onclick="swap_content('span:tropical__Cycle__SEPARATED_CODIMENSION_ONE_POLYTOPES'); return false;"><b>SEPARATED_CODIMENSION_ONE_POLYTOPES</b>:
		</a>common::IncidenceMatrix&lt;NonSymmetric&gt;<br/><div style="display: inline" id="span:tropical__Cycle__SEPARATED_CODIMENSION_ONE_POLYTOPES" class="foldit"><div class="descr"><p> An incidence matrix describing which codimension one polytope in the complex is generated by which vertices. 
 Each row corresponds to a codimension one polytope (More precisely, the i-th element represents the
 same codim 1 polytope as the i-th element of <a href="tropical.html#Cycle__CODIMENSION_ONE_POLYTOPES" onclick="unfold('span:Cycle__CODIMENSION_ONE_POLYTOPES');">CODIMENSION_ONE_POLYTOPES</a>). The indices in a row refer 
 to rows of <a href="&#10;&#9;&#9;&#9;#tropical__Cycle__SEPARATED_VERTICES" onclick="unfold('span:tropical__Cycle__SEPARATED_VERTICES');">SEPARATED_VERTICES</a>.</p></div></div></div>
                  </div>
                </li>
                <li>
                  <div class="li">
                    <div class="icon" id="icon:tropical__Cycle__SEPARATED_CONVERSION_VECTOR" onclick="swap_content('span:tropical__Cycle__SEPARATED_CONVERSION_VECTOR'); return false;">
		 	
	</div>
                    <div class="unfoldable"><a href="#" class="javalink" id="tropical__Cycle__SEPARATED_CONVERSION_VECTOR" onclick="swap_content('span:tropical__Cycle__SEPARATED_CONVERSION_VECTOR'); return false;"><b>SEPARATED_CONVERSION_VECTOR</b>:
		</a>common::Vector&lt;Int&gt;<br/><div style="display: inline" id="span:tropical__Cycle__SEPARATED_CONVERSION_VECTOR" class="foldit"><div class="descr"><p> A vector with an entry for each row in <a href="&#10;&#9;&#9;&#9;#tropical__Cycle__SEPARATED_VERTICES" onclick="unfold('span:tropical__Cycle__SEPARATED_VERTICES');">SEPARATED_VERTICES</a>. 
 More precisely, the i-th entry gives the row index of the ray in <a href="tropical.html#" onclick="unfold('span:');">VERTICES</a> 
 that is equal to the i-th row of <a href="&#10;&#9;&#9;&#9;#tropical__Cycle__SEPARATED_VERTICES" onclick="unfold('span:tropical__Cycle__SEPARATED_VERTICES');">SEPARATED_VERTICES</a>.</p></div></div></div>
                  </div>
                </li>
                <li>
                  <div class="li">
                    <div class="icon" id="icon:tropical__Cycle__SEPARATED_MAXIMAL_POLYTOPES" onclick="swap_content('span:tropical__Cycle__SEPARATED_MAXIMAL_POLYTOPES'); return false;">
		 	
	</div>
                    <div class="unfoldable"><a href="#" class="javalink" id="tropical__Cycle__SEPARATED_MAXIMAL_POLYTOPES" onclick="swap_content('span:tropical__Cycle__SEPARATED_MAXIMAL_POLYTOPES'); return false;"><b>SEPARATED_MAXIMAL_POLYTOPES</b>:
		</a>common::IncidenceMatrix&lt;NonSymmetric&gt;<br/><div style="display: inline" id="span:tropical__Cycle__SEPARATED_MAXIMAL_POLYTOPES" class="foldit"><div class="descr"><p> An incidence matrix describing which maximal polytope in the complex us generated by which rays. 
 Each row corresponds to a maximal polytope (More precisely, the i-th element represents the
 same maximal polytope as the i-th element of <a href="tropical.html#" onclick="unfold('span:');">MAXIMAL_POLYTOPES</a>). The indices in a row refer to rows 
 of <a href="&#10;&#9;&#9;&#9;#tropical__Cycle__SEPARATED_VERTICES" onclick="unfold('span:tropical__Cycle__SEPARATED_VERTICES');">SEPARATED_VERTICES</a>, i.e. the maximal polytope described by the i-th element is generated by 
 the vertices corresponding to these row indices.</p></div></div></div>
                  </div>
                </li>
                <li>
                  <div class="li">
                    <div class="icon" id="icon:tropical__Cycle__SEPARATED_VERTICES" onclick="swap_content('span:tropical__Cycle__SEPARATED_VERTICES'); return false;">
		 	
	</div>
                    <div class="unfoldable"><a href="#" class="javalink" id="tropical__Cycle__SEPARATED_VERTICES" onclick="swap_content('span:tropical__Cycle__SEPARATED_VERTICES'); return false;"><b>SEPARATED_VERTICES</b>:
		</a>common::Matrix&lt;Rational, NonSymmetric&gt;<br/><div style="display: inline" id="span:tropical__Cycle__SEPARATED_VERTICES" class="foldit"><div class="descr"><p> This is a matrix of vertices of the complex. More precisely, each ray r from <a href="tropical.html#" onclick="unfold('span:');">VERTICES</a> occurs  
 as a row in this matrix...
 - once, if r_0 = 1
 - k times, if r_0 = 0 and k is the number of equivalence classes of maximal cones 
 containing r with respect to the following relation:
 Two maximal cones m, m' containing r are equivalent, if they are equal or there exists a sequence 
 of maximal cones m = m_1,...m_r = m', such that r is contained in each m_i and each intersection 
 m_i cap m_i+1 contains at least one ray s with s_0 = 1.
 The reason for this is that, when for example specifying a piecewise affine linear function on a 
 polyhedral complex, the same far ray with x0 = 0 might be assigned two different values, if it is 
 contained in two "non-connected" maximal cones (where connectedness is to be understood as described above).
 If there is a <a href="&#10;&#9;&#9;&#9;#tropical__Cycle__LOCAL_RESTRICTION" onclick="unfold('span:tropical__Cycle__LOCAL_RESTRICTION');">LOCAL_RESTRICTION</a> the above equivalence relation is changed in such a 
 way that the affine ray s with s_0 = 1 that must be contained in the intersection of
 two subsequent cones must be a compatible ray</p></div></div></div>
                  </div>
                </li>
                <li>
                  <div class="li">
                    <div class="icon" id="icon:tropical__Cycle__WEIGHT_CONE" onclick="swap_content('span:tropical__Cycle__WEIGHT_CONE'); return false;">
		 	
	</div>
                    <div class="unfoldable"><a href="#" class="javalink" id="tropical__Cycle__WEIGHT_CONE" onclick="swap_content('span:tropical__Cycle__WEIGHT_CONE'); return false;"><b>WEIGHT_CONE</b>:
		</a>polytope::Cone&lt;Rational&gt;<br/><div style="display: inline" id="span:tropical__Cycle__WEIGHT_CONE" class="foldit"><div class="descr"><p> The intersection of <a href="&#10;&#9;&#9;&#9;#tropical__Cycle__WEIGHT_SPACE" onclick="unfold('span:tropical__Cycle__WEIGHT_SPACE');">WEIGHT_SPACE</a> with the positive orthant.</p></div></div></div>
                  </div>
                </li>
                <li>
                  <div class="li">
                    <div class="icon" id="icon:tropical__Cycle__WEIGHT_SPACE" onclick="swap_content('span:tropical__Cycle__WEIGHT_SPACE'); return false;">
		 	
	</div>
                    <div class="unfoldable"><a href="#" class="javalink" id="tropical__Cycle__WEIGHT_SPACE" onclick="swap_content('span:tropical__Cycle__WEIGHT_SPACE'); return false;"><b>WEIGHT_SPACE</b>:
		</a>common::Matrix&lt;Rational, NonSymmetric&gt;<br/><div style="display: inline" id="span:tropical__Cycle__WEIGHT_SPACE" class="foldit"><div class="descr"><p> A Z-basis (as rows) for the space of weight distributions on this 
 tropical cycle making it balanced (i.e. this cycle is irreducible, if and only if 
 <a href="&#10;&#9;&#9;&#9;#tropical__Cycle__WEIGHT_SPACE" onclick="unfold('span:tropical__Cycle__WEIGHT_SPACE');">WEIGHT_SPACE</a> has only one row and the gcd of <a href="tropical.html#Cycle__WEIGHTS" onclick="unfold('span:Cycle__WEIGHTS');">WEIGHTS</a> is 1.</p></div></div></div>
                  </div>
                </li>
                <li>
                  <div class="li">
                    <div class="icon" id="icon:tropical__Cycle__WEIGHT_SYSTEM" onclick="swap_content('span:tropical__Cycle__WEIGHT_SYSTEM'); return false;">
		 	
	</div>
                    <div class="unfoldable"><a href="#" class="javalink" id="tropical__Cycle__WEIGHT_SYSTEM" onclick="swap_content('span:tropical__Cycle__WEIGHT_SYSTEM'); return false;"><b>WEIGHT_SYSTEM</b>:
		</a>common::Matrix&lt;Rational, NonSymmetric&gt;<br/><div style="display: inline" id="span:tropical__Cycle__WEIGHT_SYSTEM" class="foldit"><div class="descr"><p> The dual of <a href="&#10;&#9;&#9;&#9;#tropical__Cycle__WEIGHT_SPACE" onclick="unfold('span:tropical__Cycle__WEIGHT_SPACE');">WEIGHT_SPACE</a>.</p></div></div></div>
                  </div>
                </li>
              </ul>
            </div>
          </div>
        </div>
      </li>
      <li>
        <div class="li">
          <div class="icon" id="icon:tropical__LinesInCubic" onclick="swap_content('span:tropical__LinesInCubic'); return false;">
		 	
	</div>
          <div class="unfoldable">
            <div class="object">
              <h3>
                <a href="#" class="javalink" id="tropical__LinesInCubic" onclick="swap_content('span:tropical__LinesInCubic'); return false;">
                  <span class="object">LinesInCubic</span>
                </a>
              </h3>
            </div>
            <div style="display: inline" id="span:tropical__LinesInCubic" class="foldit">
              <div class="descr">
                <p> This represents the result of the method <a href="tropical.html#lines_in_cubic" onclick="unfold('span:lines_in_cubic');">lines_in_cubic</a>. It contains:
 The tropical polynomial representing the surface, the surface itself as a Cycle and 
 lists of lines and families of different types, each starting with LIST_... </p>
                <p> The object also has methods, starting with array_... that return the corresponding LIST_... as a perl array.
 The different (lists of) lines can be visualized nicely with <a href="tropical.html#visualize_in_surface" onclick="unfold('span:visualize_in_surface');">visualize_in_surface</a>.</p>
              </div>
              <div class="level3"/>
              <h4>Properties of LinesInCubic</h4>
              <ul class="unfoldable">
                <li>
                  <div class="li">
                    <div class="icon" id="icon:tropical__LinesInCubic__CUBIC" onclick="swap_content('span:tropical__LinesInCubic__CUBIC'); return false;">
		 	
	</div>
                    <div class="unfoldable">
                      <a href="#" class="javalink" id="tropical__LinesInCubic__CUBIC" onclick="swap_content('span:tropical__LinesInCubic__CUBIC'); return false;"><b>CUBIC</b>:
		</a>
                      <a href="tropical.html#Cycle" onclick="unfold('span:Cycle');">Cycle</a>
                      <br/>
                      <div style="display: inline" id="span:tropical__LinesInCubic__CUBIC" class="foldit">
                        <div class="descr">
                          <p> The surface containing the lines</p>
                        </div>
                      </div>
                    </div>
                  </div>
                </li>
                <li>
                  <div class="li">
                    <div class="icon" id="icon:tropical__LinesInCubic__LIST_FAMILY_FIXED_EDGE" onclick="swap_content('span:tropical__LinesInCubic__LIST_FAMILY_FIXED_EDGE'); return false;">
		 	
	</div>
                    <div class="unfoldable">
                      <a href="#" class="javalink" id="tropical__LinesInCubic__LIST_FAMILY_FIXED_EDGE" onclick="swap_content('span:tropical__LinesInCubic__LIST_FAMILY_FIXED_EDGE'); return false;"><b>LIST_FAMILY_FIXED_EDGE</b>:
		</a>
                      <a href="tropical.html#Cycle" onclick="unfold('span:Cycle');">Cycle</a>
                      <br/>
                      <div style="display: inline" id="span:tropical__LinesInCubic__LIST_FAMILY_FIXED_EDGE" class="foldit">
                        <div class="descr">
                          <p> A list of all families of lines with (part of) the bounded edge fixed</p>
                        </div>
                      </div>
                    </div>
                  </div>
                </li>
                <li>
                  <div class="li">
                    <div class="icon" id="icon:tropical__LinesInCubic__LIST_FAMILY_FIXED_VERTEX" onclick="swap_content('span:tropical__LinesInCubic__LIST_FAMILY_FIXED_VERTEX'); return false;">
		 	
	</div>
                    <div class="unfoldable">
                      <a href="#" class="javalink" id="tropical__LinesInCubic__LIST_FAMILY_FIXED_VERTEX" onclick="swap_content('span:tropical__LinesInCubic__LIST_FAMILY_FIXED_VERTEX'); return false;"><b>LIST_FAMILY_FIXED_VERTEX</b>:
		</a>
                      <a href="tropical.html#Cycle" onclick="unfold('span:Cycle');">Cycle</a>
                      <br/>
                      <div style="display: inline" id="span:tropical__LinesInCubic__LIST_FAMILY_FIXED_VERTEX" class="foldit">
                        <div class="descr">
                          <p> A list of all families of lines without a bounded edge at a fixed vertex</p>
                        </div>
                      </div>
                    </div>
                  </div>
                </li>
                <li>
                  <div class="li">
                    <div class="icon" id="icon:tropical__LinesInCubic__LIST_FAMILY_MOVING_EDGE" onclick="swap_content('span:tropical__LinesInCubic__LIST_FAMILY_MOVING_EDGE'); return false;">
		 	
	</div>
                    <div class="unfoldable">
                      <a href="#" class="javalink" id="tropical__LinesInCubic__LIST_FAMILY_MOVING_EDGE" onclick="swap_content('span:tropical__LinesInCubic__LIST_FAMILY_MOVING_EDGE'); return false;"><b>LIST_FAMILY_MOVING_EDGE</b>:
		</a>
                      <a href="tropical.html#Cycle" onclick="unfold('span:Cycle');">Cycle</a>
                      <br/>
                      <div style="display: inline" id="span:tropical__LinesInCubic__LIST_FAMILY_MOVING_EDGE" class="foldit">
                        <div class="descr">
                          <p> A list of all families of lines with the bounded edge moving in transversal direction</p>
                        </div>
                      </div>
                    </div>
                  </div>
                </li>
                <li>
                  <div class="li">
                    <div class="icon" id="icon:tropical__LinesInCubic__LIST_FAMILY_MOVING_VERTEX" onclick="swap_content('span:tropical__LinesInCubic__LIST_FAMILY_MOVING_VERTEX'); return false;">
		 	
	</div>
                    <div class="unfoldable">
                      <a href="#" class="javalink" id="tropical__LinesInCubic__LIST_FAMILY_MOVING_VERTEX" onclick="swap_content('span:tropical__LinesInCubic__LIST_FAMILY_MOVING_VERTEX'); return false;"><b>LIST_FAMILY_MOVING_VERTEX</b>:
		</a>
                      <a href="tropical.html#Cycle" onclick="unfold('span:Cycle');">Cycle</a>
                      <br/>
                      <div style="display: inline" id="span:tropical__LinesInCubic__LIST_FAMILY_MOVING_VERTEX" class="foldit">
                        <div class="descr">
                          <p> A list of all families of lines without a bounded  edge and a moving vertex.</p>
                        </div>
                      </div>
                    </div>
                  </div>
                </li>
                <li>
                  <div class="li">
                    <div class="icon" id="icon:tropical__LinesInCubic__LIST_ISOLATED_EDGE" onclick="swap_content('span:tropical__LinesInCubic__LIST_ISOLATED_EDGE'); return false;">
		 	
	</div>
                    <div class="unfoldable">
                      <a href="#" class="javalink" id="tropical__LinesInCubic__LIST_ISOLATED_EDGE" onclick="swap_content('span:tropical__LinesInCubic__LIST_ISOLATED_EDGE'); return false;"><b>LIST_ISOLATED_EDGE</b>:
		</a>
                      <a href="tropical.html#Cycle" onclick="unfold('span:Cycle');">Cycle</a>
                      <br/>
                      <div style="display: inline" id="span:tropical__LinesInCubic__LIST_ISOLATED_EDGE" class="foldit">
                        <div class="descr">
                          <p> A list of all isolated lines with bounded edge</p>
                        </div>
                      </div>
                    </div>
                  </div>
                </li>
                <li>
                  <div class="li">
                    <div class="icon" id="icon:tropical__LinesInCubic__LIST_ISOLATED_NO_EDGE" onclick="swap_content('span:tropical__LinesInCubic__LIST_ISOLATED_NO_EDGE'); return false;">
		 	
	</div>
                    <div class="unfoldable">
                      <a href="#" class="javalink" id="tropical__LinesInCubic__LIST_ISOLATED_NO_EDGE" onclick="swap_content('span:tropical__LinesInCubic__LIST_ISOLATED_NO_EDGE'); return false;"><b>LIST_ISOLATED_NO_EDGE</b>:
		</a>
                      <a href="tropical.html#Cycle" onclick="unfold('span:Cycle');">Cycle</a>
                      <br/>
                      <div style="display: inline" id="span:tropical__LinesInCubic__LIST_ISOLATED_NO_EDGE" class="foldit">
                        <div class="descr">
                          <p> A list of all isolated lines without a bounded edge</p>
                        </div>
                      </div>
                    </div>
                  </div>
                </li>
                <li>
                  <div class="li">
                    <div class="icon" id="icon:tropical__LinesInCubic__N_FAMILIES" onclick="swap_content('span:tropical__LinesInCubic__N_FAMILIES'); return false;">
		 	
	</div>
                    <div class="unfoldable"><a href="#" class="javalink" id="tropical__LinesInCubic__N_FAMILIES" onclick="swap_content('span:tropical__LinesInCubic__N_FAMILIES'); return false;"><b>N_FAMILIES</b>:
		</a>common::Int<br/><div style="display: inline" id="span:tropical__LinesInCubic__N_FAMILIES" class="foldit"><div class="descr"><p> The total number of families in <a href="&#10;&#9;&#9;&#9;#tropical__LinesInCubic__LIST_FAMILY_FIXED_EDGE" onclick="unfold('span:tropical__LinesInCubic__LIST_FAMILY_FIXED_EDGE');">LIST_FAMILY_FIXED_EDGE</a>, <a href="&#10;&#9;&#9;&#9;#tropical__LinesInCubic__LIST_FAMILY_FIXED_VERTEX" onclick="unfold('span:tropical__LinesInCubic__LIST_FAMILY_FIXED_VERTEX');">LIST_FAMILY_FIXED_VERTEX</a>, 
 <a href="&#10;&#9;&#9;&#9;#tropical__LinesInCubic__LIST_FAMILY_MOVING_EDGE" onclick="unfold('span:tropical__LinesInCubic__LIST_FAMILY_MOVING_EDGE');">LIST_FAMILY_MOVING_EDGE</a> and <a href="&#10;&#9;&#9;&#9;#tropical__LinesInCubic__LIST_FAMILY_MOVING_VERTEX" onclick="unfold('span:tropical__LinesInCubic__LIST_FAMILY_MOVING_VERTEX');">LIST_FAMILY_MOVING_VERTEX</a></p></div></div></div>
                  </div>
                </li>
                <li>
                  <div class="li">
                    <div class="icon" id="icon:tropical__LinesInCubic__N_ISOLATED" onclick="swap_content('span:tropical__LinesInCubic__N_ISOLATED'); return false;">
		 	
	</div>
                    <div class="unfoldable"><a href="#" class="javalink" id="tropical__LinesInCubic__N_ISOLATED" onclick="swap_content('span:tropical__LinesInCubic__N_ISOLATED'); return false;"><b>N_ISOLATED</b>:
		</a>common::Int<br/><div style="display: inline" id="span:tropical__LinesInCubic__N_ISOLATED" class="foldit"><div class="descr"><p> The total number of elements in <a href="&#10;&#9;&#9;&#9;#tropical__LinesInCubic__LIST_ISOLATED_EDGE" onclick="unfold('span:tropical__LinesInCubic__LIST_ISOLATED_EDGE');">LIST_ISOLATED_EDGE</a> and <a href="&#10;&#9;&#9;&#9;#tropical__LinesInCubic__LIST_ISOLATED_NO_EDGE" onclick="unfold('span:tropical__LinesInCubic__LIST_ISOLATED_NO_EDGE');">LIST_ISOLATED_NO_EDGE</a></p></div></div></div>
                  </div>
                </li>
                <li>
                  <div class="li">
                    <div class="icon" id="icon:tropical__LinesInCubic__POLYNOMIAL" onclick="swap_content('span:tropical__LinesInCubic__POLYNOMIAL'); return false;">
		 	
	</div>
                    <div class="unfoldable"><a href="#" class="javalink" id="tropical__LinesInCubic__POLYNOMIAL" onclick="swap_content('span:tropical__LinesInCubic__POLYNOMIAL'); return false;"><b>POLYNOMIAL</b>:
		</a>common::Polynomial<br/><div style="display: inline" id="span:tropical__LinesInCubic__POLYNOMIAL" class="foldit"><div class="descr"><p> The homogeneous tropical polynomial defining the surface</p></div></div></div>
                  </div>
                </li>
              </ul>
            </div>
          </div>
        </div>
      </li>
      <li>
        <div class="li">
          <div class="icon" id="icon:tropical__Morphism" onclick="swap_content('span:tropical__Morphism'); return false;">
		 	
	</div>
          <div class="unfoldable">
            <div class="object">
              <h3>
                <a href="#" class="javalink" id="tropical__Morphism" onclick="swap_content('span:tropical__Morphism'); return false;">
                  <span class="object">Morphism</span>
                </a>
              </h3>
            </div>
            <div style="display: inline" id="span:tropical__Morphism" class="foldit">
              <div class="descr">
                <p> A morphism is a function between cycles which is locally affine linear and respects the lattices.
 It is defined by a <a href="&#10;&#9;&#9;&#9;#tropical__Morphism__DOMAIN" onclick="unfold('span:tropical__Morphism__DOMAIN');">DOMAIN</a>, which is a cycle, and values on this domain, <a href="&#10;&#9;&#9;&#9;#tropical__Morphism__VERTEX_VALUES" onclick="unfold('span:tropical__Morphism__VERTEX_VALUES');">VERTEX_VALUES</a> and
 <a href="&#10;&#9;&#9;&#9;#tropical__Morphism__LINEALITY_VALUES" onclick="unfold('span:tropical__Morphism__LINEALITY_VALUES');">LINEALITY_VALUES</a>, much like <a href="&#10;&#9;&#9;&#9;#tropical__RationalFunction" onclick="unfold('span:tropical__RationalFunction');">RationalFunction</a>.
 Alternatively, it can be defined as a global affine linear function by giving a matrix and a 
 translation vector.</p>
              </div>
              <div class="level3"/>
              <h4>Properties of Morphism</h4>
              <ul class="unfoldable">
                <li>
                  <div class="li">
                    <div class="icon" id="icon:tropical__Morphism__DOMAIN" onclick="swap_content('span:tropical__Morphism__DOMAIN'); return false;">
		 	
	</div>
                    <div class="unfoldable">
                      <a href="#" class="javalink" id="tropical__Morphism__DOMAIN" onclick="swap_content('span:tropical__Morphism__DOMAIN'); return false;"><b>DOMAIN</b>:
		</a>
                      <a href="tropical.html#Cycle" onclick="unfold('span:Cycle');">Cycle</a>
                      <br/>
                      <div style="display: inline" id="span:tropical__Morphism__DOMAIN" class="foldit">
                        <div class="descr">
                          <p> This property describes the domain of the morphism. I.e. the morphism is
 defined on this complex and is locally affine integral linear. </p>
                        </div>
                      </div>
                    </div>
                  </div>
                </li>
                <li>
                  <div class="li">
                    <div class="icon" id="icon:tropical__Morphism__IS_GLOBALLY_AFFINE_LINEAR" onclick="swap_content('span:tropical__Morphism__IS_GLOBALLY_AFFINE_LINEAR'); return false;">
		 	
	</div>
                    <div class="unfoldable"><a href="#" class="javalink" id="tropical__Morphism__IS_GLOBALLY_AFFINE_LINEAR" onclick="swap_content('span:tropical__Morphism__IS_GLOBALLY_AFFINE_LINEAR'); return false;"><b>IS_GLOBALLY_AFFINE_LINEAR</b>:
		</a>common::Bool<br/><div style="display: inline" id="span:tropical__Morphism__IS_GLOBALLY_AFFINE_LINEAR" class="foldit"><div class="descr"><p> This is TRUE, if the morphism is defined on the full projective torus by a <a href="&#10;&#9;&#9;&#9;#tropical__Morphism__MATRIX" onclick="unfold('span:tropical__Morphism__MATRIX');">MATRIX</a> and a <a href="&#10;&#9;&#9;&#9;#tropical__Morphism__TRANSLATE" onclick="unfold('span:tropical__Morphism__TRANSLATE');">TRANSLATE</a>
 The rules do not actually check for completeness of the <a href="&#10;&#9;&#9;&#9;#tropical__Morphism__DOMAIN" onclick="unfold('span:tropical__Morphism__DOMAIN');">DOMAIN</a>. 
 This property will be set to TRUE, if the morphism is only defined by <a href="&#10;&#9;&#9;&#9;#tropical__Morphism__MATRIX" onclick="unfold('span:tropical__Morphism__MATRIX');">MATRIX</a> and <a href="&#10;&#9;&#9;&#9;#tropical__Morphism__TRANSLATE" onclick="unfold('span:tropical__Morphism__TRANSLATE');">TRANSLATE</a>, 
 otherwise it is false (or you can set it upon creation).</p></div></div></div>
                  </div>
                </li>
                <li>
                  <div class="li">
                    <div class="icon" id="icon:tropical__Morphism__LINEALITY_VALUES" onclick="swap_content('span:tropical__Morphism__LINEALITY_VALUES'); return false;">
		 	
	</div>
                    <div class="unfoldable"><a href="#" class="javalink" id="tropical__Morphism__LINEALITY_VALUES" onclick="swap_content('span:tropical__Morphism__LINEALITY_VALUES'); return false;"><b>LINEALITY_VALUES</b>:
		</a>common::Matrix&lt;Rational, NonSymmetric&gt;<br/><div style="display: inline" id="span:tropical__Morphism__LINEALITY_VALUES" class="foldit"><div class="descr"><p> The vector in row i describes the function value (slope) of <a href="&#10;&#9;&#9;&#9;#tropical__Morphism__DOMAIN" onclick="unfold('span:tropical__Morphism__DOMAIN');">DOMAIN</a>-&gt;LINEALITY_SPACE-&gt;row(i)</p></div></div></div>
                  </div>
                </li>
                <li>
                  <div class="li">
                    <div class="icon" id="icon:tropical__Morphism__MATRIX" onclick="swap_content('span:tropical__Morphism__MATRIX'); return false;">
		 	
	</div>
                    <div class="unfoldable"><a href="#" class="javalink" id="tropical__Morphism__MATRIX" onclick="swap_content('span:tropical__Morphism__MATRIX'); return false;"><b>MATRIX</b>:
		</a>common::Matrix&lt;Rational, NonSymmetric&gt;<br/><div style="display: inline" id="span:tropical__Morphism__MATRIX" class="foldit"><div class="descr"><p> If the morphism is a global affine linear map x |-&gt; Ax+v, then this contains the matrix A.
 Note that this must be well-defined on tropical projective coordinates, so the sum of the columns of A
 must be a multiple of the (1,..,1)-vector.
 If <a href="&#10;&#9;&#9;&#9;#tropical__Morphism__TRANSLATE" onclick="unfold('span:tropical__Morphism__TRANSLATE');">TRANSLATE</a> is set, but this property is not set, then it is the identity by default.</p></div></div></div>
                  </div>
                </li>
                <li>
                  <div class="li">
                    <div class="icon" id="icon:tropical__Morphism__TRANSLATE" onclick="swap_content('span:tropical__Morphism__TRANSLATE'); return false;">
		 	
	</div>
                    <div class="unfoldable"><a href="#" class="javalink" id="tropical__Morphism__TRANSLATE" onclick="swap_content('span:tropical__Morphism__TRANSLATE'); return false;"><b>TRANSLATE</b>:
		</a>common::Vector&lt;Rational&gt;<br/><div style="display: inline" id="span:tropical__Morphism__TRANSLATE" class="foldit"><div class="descr"><p> If the morphism is a global affine linear map x |-&gt; Ax+v, then this contains the translation vector 
 v. If <a href="&#10;&#9;&#9;&#9;#tropical__Morphism__MATRIX" onclick="unfold('span:tropical__Morphism__MATRIX');">MATRIX</a> is set, but this property is not set, then it is the zero vector by default.</p></div></div></div>
                  </div>
                </li>
                <li>
                  <div class="li">
                    <div class="icon" id="icon:tropical__Morphism__VERTEX_VALUES" onclick="swap_content('span:tropical__Morphism__VERTEX_VALUES'); return false;">
		 	
	</div>
                    <div class="unfoldable"><a href="#" class="javalink" id="tropical__Morphism__VERTEX_VALUES" onclick="swap_content('span:tropical__Morphism__VERTEX_VALUES'); return false;"><b>VERTEX_VALUES</b>:
		</a>common::Matrix&lt;Rational, NonSymmetric&gt;<br/><div style="display: inline" id="span:tropical__Morphism__VERTEX_VALUES" class="foldit"><div class="descr"><p> The vector at row i describes the function value of vertex <a href="&#10;&#9;&#9;&#9;#tropical__Morphism__DOMAIN" onclick="unfold('span:tropical__Morphism__DOMAIN');">DOMAIN</a>-&gt;SEPARATED_VERTICES-&gt;row(i). 
 (In tropical homogenous coordinates, but without leading coordinate).
 More precisely, if the corresponding vertex is not a far ray, it describes its function value. 
 If it is a directional ray, it describes the slope on that ray.</p></div></div></div>
                  </div>
                </li>
              </ul>
            </div>
          </div>
        </div>
      </li>
      <li>
        <div class="li">
          <div class="icon" id="icon:tropical__RationalCurve" onclick="swap_content('span:tropical__RationalCurve'); return false;">
		 	
	</div>
          <div class="unfoldable">
            <div class="object">
              <h3>
                <a href="#" class="javalink" id="tropical__RationalCurve" onclick="swap_content('span:tropical__RationalCurve'); return false;">
                  <span class="object">RationalCurve</span>
                </a>
              </h3>
            </div>
            <div style="display: inline" id="span:tropical__RationalCurve" class="foldit">
              <div class="descr">
                <p> An n-marked rational curve, identified by its <a href="&#10;&#9;&#9;&#9;#tropical__RationalCurve__SETS" onclick="unfold('span:tropical__RationalCurve__SETS');">SETS</a>, i.e. its partitions of {1,...,n} and
 its <a href="tropical.html#Hypersurface__COEFFICIENTS" onclick="unfold('span:Hypersurface__COEFFICIENTS');">COEFFICIENTS</a>, i.e. the lengths of the corresponding edges.</p>
              </div>
              <div class="level3"/>
              <h4>Properties of RationalCurve</h4>
              <ul class="unfoldable">
                <li>
                  <div class="li">
                    <div class="icon" id="icon:tropical__RationalCurve__COEFFS" onclick="swap_content('span:tropical__RationalCurve__COEFFS'); return false;">
		 	
	</div>
                    <div class="unfoldable"><a href="#" class="javalink" id="tropical__RationalCurve__COEFFS" onclick="swap_content('span:tropical__RationalCurve__COEFFS'); return false;"><b>COEFFS</b>:
		</a>common::Vector&lt;Rational&gt;<br/><div style="display: inline" id="span:tropical__RationalCurve__COEFFS" class="foldit"><div class="descr"><p> A list of positive rational coefficients. The list should have the same length as <a href="&#10;&#9;&#9;&#9;#tropical__RationalCurve__SETS" onclick="unfold('span:tropical__RationalCurve__SETS');">SETS</a> and
 contain only entries &gt; 0. The i-th entry then gives the length of the bounded edge defined by
 the i-th partition. If you're not sure if all your coefficients are &gt; 0, use <a href="&#10;&#9;&#9;&#9;#tropical__RationalCurve__INPUT_SETS" onclick="unfold('span:tropical__RationalCurve__INPUT_SETS');">INPUT_SETS</a> and
 <a href="&#10;&#9;&#9;&#9;#tropical__RationalCurve__INPUT_COEFFS" onclick="unfold('span:tropical__RationalCurve__INPUT_COEFFS');">INPUT_COEFFS</a> instead.
 Note that the zero curve (i.e. no bounded edges, only leaves) is represented by one empty set
 with corresponding lenghth 0.</p></div></div></div>
                  </div>
                </li>
                <li>
                  <div class="li">
                    <div class="icon" id="icon:tropical__RationalCurve__GRAPH" onclick="swap_content('span:tropical__RationalCurve__GRAPH'); return false;">
		 	
	</div>
                    <div class="unfoldable"><a href="#" class="javalink" id="tropical__RationalCurve__GRAPH" onclick="swap_content('span:tropical__RationalCurve__GRAPH'); return false;"><b>GRAPH</b>:
		</a>graph::Graph&lt;Undirected&gt;<br/><div style="display: inline" id="span:tropical__RationalCurve__GRAPH" class="foldit"><div class="descr"><p> Contains the abstract graph (non-metric) corresponding to the curve. All unbounded leaves are modelled
 as bounded edges.
 The vertices at the ends of the "leaves" are always the first <a href="&#10;&#9;&#9;&#9;#tropical__RationalCurve__N_LEAVES" onclick="unfold('span:tropical__RationalCurve__N_LEAVES');">N_LEAVES</a> vertices.</p></div></div></div>
                  </div>
                </li>
                <li>
                  <div class="li">
                    <div class="icon" id="icon:tropical__RationalCurve__GRAPH_EDGE_LENGTHS" onclick="swap_content('span:tropical__RationalCurve__GRAPH_EDGE_LENGTHS'); return false;">
		 	
	</div>
                    <div class="unfoldable"><a href="#" class="javalink" id="tropical__RationalCurve__GRAPH_EDGE_LENGTHS" onclick="swap_content('span:tropical__RationalCurve__GRAPH_EDGE_LENGTHS'); return false;"><b>GRAPH_EDGE_LENGTHS</b>:
		</a>common::Vector&lt;Rational&gt;<br/><div style="display: inline" id="span:tropical__RationalCurve__GRAPH_EDGE_LENGTHS" class="foldit"><div class="descr"><p> Contains the lengths of the edges of <a href="&#10;&#9;&#9;&#9;#tropical__RationalCurve__GRAPH" onclick="unfold('span:tropical__RationalCurve__GRAPH');">GRAPH</a> that represent bounded edges of the curve.
 The coefficients appear in the order that the corr. edges appear in <a href="tropical.html#" onclick="unfold('span:');">EDGES</a>.</p></div></div></div>
                  </div>
                </li>
                <li>
                  <div class="li">
                    <div class="icon" id="icon:tropical__RationalCurve__INPUT_COEFFS" onclick="swap_content('span:tropical__RationalCurve__INPUT_COEFFS'); return false;">
		 	
	</div>
                    <div class="unfoldable"><a href="#" class="javalink" id="tropical__RationalCurve__INPUT_COEFFS" onclick="swap_content('span:tropical__RationalCurve__INPUT_COEFFS'); return false;"><b>INPUT_COEFFS</b>:
		</a>common::Vector&lt;Rational&gt;<br/><div style="display: inline" id="span:tropical__RationalCurve__INPUT_COEFFS" class="foldit"><div class="descr"><p> Same as <a href="&#10;&#9;&#9;&#9;#tropical__RationalCurve__COEFFS" onclick="unfold('span:tropical__RationalCurve__COEFFS');">COEFFS</a>, except that entries may be &lt;=0. This should have the same length as <a href="&#10;&#9;&#9;&#9;#tropical__RationalCurve__INPUT_SETS" onclick="unfold('span:tropical__RationalCurve__INPUT_SETS');">INPUT_SETS</a>.</p></div></div></div>
                  </div>
                </li>
                <li>
                  <div class="li">
                    <div class="icon" id="icon:tropical__RationalCurve__INPUT_SETS" onclick="swap_content('span:tropical__RationalCurve__INPUT_SETS'); return false;">
		 	
	</div>
                    <div class="unfoldable"><a href="#" class="javalink" id="tropical__RationalCurve__INPUT_SETS" onclick="swap_content('span:tropical__RationalCurve__INPUT_SETS'); return false;"><b>INPUT_SETS</b>:
		</a>common::IncidenceMatrix&lt;NonSymmetric&gt;<br/><div style="display: inline" id="span:tropical__RationalCurve__INPUT_SETS" class="foldit"><div class="descr"><p> Same as <a href="&#10;&#9;&#9;&#9;#tropical__RationalCurve__SETS" onclick="unfold('span:tropical__RationalCurve__SETS');">SETS</a>, except that sets may appear several times.</p></div></div></div>
                  </div>
                </li>
                <li>
                  <div class="li">
                    <div class="icon" id="icon:tropical__RationalCurve__INPUT_STRING" onclick="swap_content('span:tropical__RationalCurve__INPUT_STRING'); return false;">
		 	
	</div>
                    <div class="unfoldable"><a href="#" class="javalink" id="tropical__RationalCurve__INPUT_STRING" onclick="swap_content('span:tropical__RationalCurve__INPUT_STRING'); return false;"><b>INPUT_STRING</b>:
		</a>common::String<br/><div style="display: inline" id="span:tropical__RationalCurve__INPUT_STRING" class="foldit"><div class="descr"><p> This property can also be used to define a rational curve: A linear combination of partitions is
 given as a string, using the following syntax:
 A partition is given as a subset of {1,..,n} and written as a comma-separated list of leaf
 indices in round brackets, e.g. "(1,2,5)"
 A linear combination can be created using rational numbers, "+","+" and "-" in the obvious way, 
 e.g. "2*(1,2,5) + 1*(3,4,7) - 2(1,2) (The "*" is optional)
 Of course, each set should contain at least two elements. If you don't specify N_LEAVES, it is
 set to be the largest leaf index occuring in the sets.
 Partitions needn't be irredundant and coefficients can be any rational number. If the resulting 
 element is not in the moduli space, an error is thrown.</p></div></div></div>
                  </div>
                </li>
                <li>
                  <div class="li">
                    <div class="icon" id="icon:tropical__RationalCurve__NODES_BY_LEAVES" onclick="swap_content('span:tropical__RationalCurve__NODES_BY_LEAVES'); return false;">
		 	
	</div>
                    <div class="unfoldable"><a href="#" class="javalink" id="tropical__RationalCurve__NODES_BY_LEAVES" onclick="swap_content('span:tropical__RationalCurve__NODES_BY_LEAVES'); return false;"><b>NODES_BY_LEAVES</b>:
		</a>common::IncidenceMatrix&lt;NonSymmetric&gt;<br/><div style="display: inline" id="span:tropical__RationalCurve__NODES_BY_LEAVES" class="foldit"><div class="descr"><p> This incidence matrix gives a list of the vertices of the curve
 Each row corresponds to a vertex and contains as a set the [[LEAVES] that are
 attached to that vertex (again, counting from 1!)</p></div></div></div>
                  </div>
                </li>
                <li>
                  <div class="li">
                    <div class="icon" id="icon:tropical__RationalCurve__NODES_BY_SETS" onclick="swap_content('span:tropical__RationalCurve__NODES_BY_SETS'); return false;">
		 	
	</div>
                    <div class="unfoldable"><a href="#" class="javalink" id="tropical__RationalCurve__NODES_BY_SETS" onclick="swap_content('span:tropical__RationalCurve__NODES_BY_SETS'); return false;"><b>NODES_BY_SETS</b>:
		</a>common::IncidenceMatrix&lt;NonSymmetric&gt;<br/><div style="display: inline" id="span:tropical__RationalCurve__NODES_BY_SETS" class="foldit"><div class="descr"><p> This incidence matrix gives a list of the vertices of the curve
 Each row corresponds to a vertex and contains as a set the row indices of the <a href="&#10;&#9;&#9;&#9;#tropical__RationalCurve__SETS" onclick="unfold('span:tropical__RationalCurve__SETS');">SETS</a> that
 correspond to edges attached to that vertex</p></div></div></div>
                  </div>
                </li>
                <li>
                  <div class="li">
                    <div class="icon" id="icon:tropical__RationalCurve__NODE_DEGREES" onclick="swap_content('span:tropical__RationalCurve__NODE_DEGREES'); return false;">
		 	
	</div>
                    <div class="unfoldable"><a href="#" class="javalink" id="tropical__RationalCurve__NODE_DEGREES" onclick="swap_content('span:tropical__RationalCurve__NODE_DEGREES'); return false;"><b>NODE_DEGREES</b>:
		</a>common::Vector&lt;Int&gt;<br/><div style="display: inline" id="span:tropical__RationalCurve__NODE_DEGREES" class="foldit"><div class="descr"><p> This gives a list of the vertices of the curve in terms of their valences
 They appear in the same order as in <a href="&#10;&#9;&#9;&#9;#tropical__RationalCurve__NODES_BY_LEAVES" onclick="unfold('span:tropical__RationalCurve__NODES_BY_LEAVES');">NODES_BY_LEAVES</a> or <a href="&#10;&#9;&#9;&#9;#tropical__RationalCurve__NODES_BY_SETS" onclick="unfold('span:tropical__RationalCurve__NODES_BY_SETS');">NODES_BY_SETS</a></p></div></div></div>
                  </div>
                </li>
                <li>
                  <div class="li">
                    <div class="icon" id="icon:tropical__RationalCurve__N_LEAVES" onclick="swap_content('span:tropical__RationalCurve__N_LEAVES'); return false;">
		 	
	</div>
                    <div class="unfoldable"><a href="#" class="javalink" id="tropical__RationalCurve__N_LEAVES" onclick="swap_content('span:tropical__RationalCurve__N_LEAVES'); return false;"><b>N_LEAVES</b>:
		</a>common::Int<br/><div style="display: inline" id="span:tropical__RationalCurve__N_LEAVES" class="foldit"><div class="descr"><p> The number of leaves of the rational curve.</p></div></div></div>
                  </div>
                </li>
                <li>
                  <div class="li">
                    <div class="icon" id="icon:tropical__RationalCurve__SETS" onclick="swap_content('span:tropical__RationalCurve__SETS'); return false;">
		 	
	</div>
                    <div class="unfoldable"><a href="#" class="javalink" id="tropical__RationalCurve__SETS" onclick="swap_content('span:tropical__RationalCurve__SETS'); return false;"><b>SETS</b>:
		</a>common::IncidenceMatrix&lt;NonSymmetric&gt;<br/><div style="display: inline" id="span:tropical__RationalCurve__SETS" class="foldit"><div class="descr"><p> A list of partitions of [n] that define the tree of the curve: For each bounded edge we have 
 the corresponding partition of the n leaves. These should be irredundant. If you want to input
 a possibly redundant list, use <a href="&#10;&#9;&#9;&#9;#tropical__RationalCurve__INPUT_SETS" onclick="unfold('span:tropical__RationalCurve__INPUT_SETS');">INPUT_SETS</a> and <a href="&#10;&#9;&#9;&#9;#tropical__RationalCurve__INPUT_COEFFS" onclick="unfold('span:tropical__RationalCurve__INPUT_COEFFS');">INPUT_COEFFS</a> instead. The number of marked leaves
 should always be given by <a href="&#10;&#9;&#9;&#9;#tropical__RationalCurve__N_LEAVES" onclick="unfold('span:tropical__RationalCurve__N_LEAVES');">N_LEAVES</a>. The sets are subsets of {1,...,n} (NOT {0,..,n-1}!)
 Note that the zero curve (i.e. no bounded edges, only leaves) is represented by one empty set
 with corresponding lenghth 0.</p></div></div></div>
                  </div>
                </li>
              </ul>
            </div>
          </div>
        </div>
      </li>
      <li>
        <div class="li">
          <div class="icon" id="icon:tropical__RationalFunction" onclick="swap_content('span:tropical__RationalFunction'); return false;">
		 	
	</div>
          <div class="unfoldable">
            <div class="object">
              <h3>
                <a href="#" class="javalink" id="tropical__RationalFunction" onclick="swap_content('span:tropical__RationalFunction'); return false;">
                  <span class="object">RationalFunction</span>
                </a>
              </h3>
            </div>
            <div style="display: inline" id="span:tropical__RationalFunction" class="foldit">
              <div class="descr">
                <p> A rational function on a polyhedral complex. It can be described by giving its <a href="&#10;&#9;&#9;&#9;#tropical__RationalFunction__DOMAIN" onclick="unfold('span:tropical__RationalFunction__DOMAIN');">DOMAIN</a>,
 a Cycle, and values on this domain - which are encoded in the properties <a href="&#10;&#9;&#9;&#9;#tropical__RationalFunction__VERTEX_VALUES" onclick="unfold('span:tropical__RationalFunction__VERTEX_VALUES');">VERTEX_VALUES</a> and 
 <a href="&#10;&#9;&#9;&#9;#tropical__RationalFunction__LINEALITY_VALUES" onclick="unfold('span:tropical__RationalFunction__LINEALITY_VALUES');">LINEALITY_VALUES</a>.
 Alternatively, it can be defined by a tropical quotient of homogeneous tropical polynomials of the 
 same degree i.e. by giving <a href="&#10;&#9;&#9;&#9;#tropical__RationalFunction__NUMERATOR" onclick="unfold('span:tropical__RationalFunction__NUMERATOR');">NUMERATOR</a> and <a href="&#10;&#9;&#9;&#9;#tropical__RationalFunction__DENOMINATOR" onclick="unfold('span:tropical__RationalFunction__DENOMINATOR');">DENOMINATOR</a>. 
 A <a href="&#10;&#9;&#9;&#9;#tropical__RationalFunction__DOMAIN" onclick="unfold('span:tropical__RationalFunction__DOMAIN');">DOMAIN</a> can be defined additionally (though one should take
 care that both functions are actually piecewise affine linear on the cells), otherwise it will be
 computed as the common refinement of the domains of affine linearity of the two polynomials.
 Note: This has nothing to do with common's RationalFunction (which is univariate). If you want to 
 access that type or use this type from another application, be sure to prepend the appropriate 
 namespace identifier.</p>
              </div>
              <div class="level3"/>
              <h4>Properties of RationalFunction</h4>
              <ul class="unfoldable">
                <li>
                  <div class="li">
                    <div class="icon" id="icon:tropical__RationalFunction__DENOMINATOR" onclick="swap_content('span:tropical__RationalFunction__DENOMINATOR'); return false;">
		 	
	</div>
                    <div class="unfoldable"><a href="#" class="javalink" id="tropical__RationalFunction__DENOMINATOR" onclick="swap_content('span:tropical__RationalFunction__DENOMINATOR'); return false;"><b>DENOMINATOR</b>:
		</a>common::Polynomial<br/><div style="display: inline" id="span:tropical__RationalFunction__DENOMINATOR" class="foldit"><div class="descr"><p> When representing the function as a quotient of tropical polynomials, this
 is the denominator. Should be a homogeneous polynomial of the same degree as <a href="&#10;&#9;&#9;&#9;#tropical__RationalFunction__NUMERATOR" onclick="unfold('span:tropical__RationalFunction__NUMERATOR');">NUMERATOR</a>. </p></div></div></div>
                  </div>
                </li>
                <li>
                  <div class="li">
                    <div class="icon" id="icon:tropical__RationalFunction__DOMAIN" onclick="swap_content('span:tropical__RationalFunction__DOMAIN'); return false;">
		 	
	</div>
                    <div class="unfoldable">
                      <a href="#" class="javalink" id="tropical__RationalFunction__DOMAIN" onclick="swap_content('span:tropical__RationalFunction__DOMAIN'); return false;"><b>DOMAIN</b>:
		</a>
                      <a href="tropical.html#Cycle" onclick="unfold('span:Cycle');">Cycle</a>
                      <br/>
                      <div style="display: inline" id="span:tropical__RationalFunction__DOMAIN" class="foldit">
                        <div class="descr">
                          <p> This property describes the affine linearity domains of the function. I.e. the function is
 affine integral linear on each maximal polytope of <a href="&#10;&#9;&#9;&#9;#tropical__RationalFunction__DOMAIN" onclick="unfold('span:tropical__RationalFunction__DOMAIN');">DOMAIN</a>. </p>
                        </div>
                      </div>
                    </div>
                  </div>
                </li>
                <li>
                  <div class="li">
                    <div class="icon" id="icon:tropical__RationalFunction__IS_GLOBALLY_DEFINED" onclick="swap_content('span:tropical__RationalFunction__IS_GLOBALLY_DEFINED'); return false;">
		 	
	</div>
                    <div class="unfoldable"><a href="#" class="javalink" id="tropical__RationalFunction__IS_GLOBALLY_DEFINED" onclick="swap_content('span:tropical__RationalFunction__IS_GLOBALLY_DEFINED'); return false;"><b>IS_GLOBALLY_DEFINED</b>:
		</a>common::Bool<br/><div style="display: inline" id="span:tropical__RationalFunction__IS_GLOBALLY_DEFINED" class="foldit"><div class="descr"><p> This is TRUE, if the function is defined on the full projective torus by a <a href="&#10;&#9;&#9;&#9;#tropical__RationalFunction__NUMERATOR" onclick="unfold('span:tropical__RationalFunction__NUMERATOR');">NUMERATOR</a> and a 
 <a href="&#10;&#9;&#9;&#9;#tropical__RationalFunction__DENOMINATOR" onclick="unfold('span:tropical__RationalFunction__DENOMINATOR');">DENOMINATOR</a>.
 The rules do not actually check for completeness of the <a href="&#10;&#9;&#9;&#9;#tropical__RationalFunction__DOMAIN" onclick="unfold('span:tropical__RationalFunction__DOMAIN');">DOMAIN</a>. This property will be set to true, 
 if the function is created only via <a href="&#10;&#9;&#9;&#9;#tropical__RationalFunction__NUMERATOR" onclick="unfold('span:tropical__RationalFunction__NUMERATOR');">NUMERATOR</a> and <a href="&#10;&#9;&#9;&#9;#tropical__RationalFunction__DENOMINATOR" onclick="unfold('span:tropical__RationalFunction__DENOMINATOR');">DENOMINATOR</a>. Otherwise
 it will be set to FALSE (or you can set it manually upon creation).</p></div></div></div>
                  </div>
                </li>
                <li>
                  <div class="li">
                    <div class="icon" id="icon:tropical__RationalFunction__LINEALITY_VALUES" onclick="swap_content('span:tropical__RationalFunction__LINEALITY_VALUES'); return false;">
		 	
	</div>
                    <div class="unfoldable"><a href="#" class="javalink" id="tropical__RationalFunction__LINEALITY_VALUES" onclick="swap_content('span:tropical__RationalFunction__LINEALITY_VALUES'); return false;"><b>LINEALITY_VALUES</b>:
		</a>common::Vector&lt;Rational&gt;<br/><div style="display: inline" id="span:tropical__RationalFunction__LINEALITY_VALUES" class="foldit"><div class="descr"><p> The value at index i describes the function value of <a href="&#10;&#9;&#9;&#9;#tropical__RationalFunction__DOMAIN" onclick="unfold('span:tropical__RationalFunction__DOMAIN');">DOMAIN</a>-&gt;LINEALITY_SPACE-&gt;row(i)</p></div></div></div>
                  </div>
                </li>
                <li>
                  <div class="li">
                    <div class="icon" id="icon:tropical__RationalFunction__NUMERATOR" onclick="swap_content('span:tropical__RationalFunction__NUMERATOR'); return false;">
		 	
	</div>
                    <div class="unfoldable"><a href="#" class="javalink" id="tropical__RationalFunction__NUMERATOR" onclick="swap_content('span:tropical__RationalFunction__NUMERATOR'); return false;"><b>NUMERATOR</b>:
		</a>common::Polynomial<br/><div style="display: inline" id="span:tropical__RationalFunction__NUMERATOR" class="foldit"><div class="descr"><p> When representing the function as a quotient of tropical polynomials, this
 is the numerator. Should be a homogeneous polynomial of the same degree as <a href="&#10;&#9;&#9;&#9;#tropical__RationalFunction__DENOMINATOR" onclick="unfold('span:tropical__RationalFunction__DENOMINATOR');">DENOMINATOR</a>.</p></div></div></div>
                  </div>
                </li>
                <li>
                  <div class="li">
                    <div class="icon" id="icon:tropical__RationalFunction__POWER" onclick="swap_content('span:tropical__RationalFunction__POWER'); return false;">
		 	
	</div>
                    <div class="unfoldable"><a href="#" class="javalink" id="tropical__RationalFunction__POWER" onclick="swap_content('span:tropical__RationalFunction__POWER'); return false;"><b>POWER</b>:
		</a>common::Int<br/><div style="display: inline" id="span:tropical__RationalFunction__POWER" class="foldit"><div class="descr"><p> This is an internally used property that should not actually be set by the user. When creating a
 rational function with the ^-operator, this property is set to the exponent. 
 The semantics is that when computing a divisor, this function should be applied so many times
 The usual application of this is a call to divisor($X, $f^4) or something similar.
 Warning: This property is not stored if the RationalFunction object is saved. Nor should be assumed
 to be preserved during any kind of arithmetic or restricting operation.</p></div></div></div>
                  </div>
                </li>
                <li>
                  <div class="li">
                    <div class="icon" id="icon:tropical__RationalFunction__VERTEX_VALUES" onclick="swap_content('span:tropical__RationalFunction__VERTEX_VALUES'); return false;">
		 	
	</div>
                    <div class="unfoldable"><a href="#" class="javalink" id="tropical__RationalFunction__VERTEX_VALUES" onclick="swap_content('span:tropical__RationalFunction__VERTEX_VALUES'); return false;"><b>VERTEX_VALUES</b>:
		</a>common::Vector&lt;Rational&gt;<br/><div style="display: inline" id="span:tropical__RationalFunction__VERTEX_VALUES" class="foldit"><div class="descr"><p> The value at index i describes the function value at <a href="&#10;&#9;&#9;&#9;#tropical__RationalFunction__DOMAIN" onclick="unfold('span:tropical__RationalFunction__DOMAIN');">DOMAIN</a>-&gt;SEPARATED_VERTICES-&gt;row(i). 
 More precisely, if the corresponding vertex is not a far ray, it describes its function value. 
 If it is a directional ray, it describes the slope on that ray.</p></div></div></div>
                  </div>
                </li>
              </ul>
            </div>
          </div>
        </div>
      </li>
    </ul>
  </div>
</div>
