/**
This is written by *César Pairetti*, and is available [(here)](http://basilisk.fr/sandbox/pairetti/bag_mode/adapt_wavelet_limited.h). I keep it in the same folder for completeness. 
*/
#define TREE 1

struct Adapt_limited {
  scalar * slist; // list of scalars
  double * max;   // tolerance for each scalar
  int (*MLFun)(double,double,double);   // give maximum level as a field
  int minlevel;   // minimum level of refinement (default 1)
  scalar * list;  // list of fields to update (default all)
};

trace
astats adapt_wavelet_limited (struct Adapt_limited p)
{
  scalar * listcm = NULL;

  if (is_constant(cm)) {
    if (p.list == NULL)
      p.list = all;
    restriction (p.slist);
  }
  else {
    if (p.list == NULL) {
      listcm = list_concat (NULL, {cm,fm});
      for (scalar s in all)
	listcm = list_add (listcm, s);
      p.list = listcm;
    }
    scalar * listr = list_concat (p.slist, {cm});
    restriction (listr);
    free (listr);
  }

  astats st = {0, 0};
  scalar * listc = NULL;
  for (scalar s in p.list)
    if (!is_constant(s) && s.restriction != no_restriction)
      listc = list_add (listc, s);

  // refinement
  if (p.minlevel < 1)
    p.minlevel = 1;
  tree->refined.n = 0;
  static const int refined = 1 << user, too_fine = 1 << (user + 1);
  foreach_cell() {
    int cellMAX = p.MLFun(x,y,z);
    if (is_active(cell)) {
      static const int too_coarse = 1 << (user + 2);
      if (is_leaf (cell)) {
	if (cell.flags & too_coarse) {
	  cell.flags &= ~too_coarse;
	  refine_cell (point, listc, refined, &tree->refined);
	  st.nf++;
	}
	continue;
      }
      else { // !is_leaf (cell)
	if (cell.flags & refined) {
	  // cell has already been refined, skip its children
	  cell.flags &= ~too_coarse;
	  continue;
	}
	// check whether the cell or any of its children is local
	bool local = is_local(cell);
	if (!local)
	  foreach_child()
	    if (is_local(cell))
	      local = true, break;
	if (local) {
	  int i = 0;
	  static const int just_fine = 1 << (user + 3);
	  for (scalar s in p.slist) {
	    double max = p.max[i++], sc[1 << dimension];
	    int c = 0;
	    foreach_child()
	      sc[c++] = s[];
	    s.prolongation (point, s);
	    c = 0;
	    foreach_child() {
	      double e = fabs(sc[c] - s[]);
	      if (e > max && level < cellMAX) {
		cell.flags &= ~too_fine;
		cell.flags |= too_coarse;
	      }
	      else if ((e <= max/1.5 || level > cellMAX) &&
		       !(cell.flags & (too_coarse|just_fine))) {
		if (level >= p.minlevel)
		  cell.flags |= too_fine;
	      }
	      else if (!(cell.flags & too_coarse)) {
		cell.flags &= ~too_fine;
		cell.flags |= just_fine;
	      }
	      s[] = sc[c++];
	    }
	  }
	  foreach_child() {
	    cell.flags &= ~just_fine;
	    if (!is_leaf(cell)) {
	      cell.flags &= ~too_coarse;
	      if (level >= cellMAX)
		cell.flags |= too_fine;
	    }
	    else if (!is_active(cell))
	      cell.flags &= ~too_coarse;
	  }
	}
      }
    }
    else // inactive cell
      continue;
  }
  mpi_boundary_refine (listc);

  // coarsening
  // the loop below is only necessary to ensure symmetry of 2:1 constraint
  for (int l = depth(); l >= p.minlevel; l--) {
    foreach_cell()
      if (!is_boundary(cell)) {
	if (level == l) {
	  if (!is_leaf(cell)) {
	    if (cell.flags & refined)
	      // cell was refined previously, unset the flag
	      cell.flags &= ~(refined|too_fine);
	    else if (cell.flags & too_fine) {
	      if (is_local(cell) && coarsen_cell (point, listc))
		st.nc++;
	      cell.flags &= ~too_fine; // do not coarsen parent
	    }
	  }
	  if (cell.flags & too_fine)
	    cell.flags &= ~too_fine;
	  else if (aparent(0).flags & too_fine)
	    aparent(0).flags &= ~too_fine;
	  continue;
	}
	else if (is_leaf(cell))
	  continue;
      }
    mpi_boundary_coarsen (l, too_fine);
  }
  free (listc);

  mpi_all_reduce (st.nf, MPI_INT, MPI_SUM);
  mpi_all_reduce (st.nc, MPI_INT, MPI_SUM);
  if (st.nc || st.nf)
    mpi_boundary_update (p.list);
  free (listcm);

  return st;
}
